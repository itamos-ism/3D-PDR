"""Drive and read the RT-tool (RTtool) post-processor.

RT-tool reads ``paramsRT.dat`` (directory, prefix, linewidth, redshift) — all of
which we fill automatically from the model's ``params.dat`` — and writes, into
the model directory:

* ``Tr_<prefix>.<coolant>.dat``  : Ntot, N(H2), Tgas, Av, N(coolant), then the
  brightness (antenna) temperature Tr for the first ``coolev`` ladder
  transitions (C+ -> 1, C/O -> 2, others -> 10).
* ``tau_<prefix>.<coolant>.dat`` : the matching optical depths (row-aligned,
  no leading columns).
* ``Tex_<prefix>.<coolant>.dat`` : the matching excitation temperatures
  (row-aligned, no leading columns).

Transition column ``c`` (0-based) is the rotational/ladder line ``(c+1)-(c)``.
Per-transition frequencies don't vary with depth, so they aren't repeated in
those per-depth files — we get them once from the per-line summary table
(Line, Freq, Tex, tau, Tr) that RT-tool prints to stdout, which we parse from
the captured run output.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from . import fortread, paths

# Number of lines to show per coolant in the summary table (sorted by freq).
TABLE_LIMIT = {"C+": 1, "C": 2, "O": 2}
DEFAULT_LIMIT = 10


def table_limit(cname: str) -> int:
    return TABLE_LIMIT.get(cname, DEFAULT_LIMIT)


def transition_label(col: int) -> str:
    """0-based Tr/tau column -> ladder transition label, e.g. 0 -> '1-0'."""
    return f"{col + 1}-{col}"


# Preferred display order for coolants: the bright PDR diagnostics first, then
# everything else alphabetically.
COOLANT_ORDER = ["CO", "C", "C+", "O"]


def order_coolants(names: list[str]) -> list[str]:
    """Reorder coolant names as CO, C, C⁺, O, then the rest (alphabetical)."""
    head = [c for c in COOLANT_ORDER if c in names]
    tail = sorted(c for c in names if c not in COOLANT_ORDER)
    return head + tail


# Astrophysical labels for the atomic fine-structure / [CI] lines, keyed by
# (coolant, 0-based column). The full form is used in the line tables; the short
# form (just the wavelength) in the per-coolant Tr/tau panels.
_LINE_LABELS_FULL = {
    ("C+", 0): "[CII] 158μm",
    ("C", 0): "[CI] (1-0)",
    ("C", 1): "[CI] (2-1)",
    ("O", 0): "[OI] 63μm",
    ("O", 1): "[OI] 146μm",
}
_TRANS_SHORT = {
    ("C+", 0): "158μm",
    ("O", 0): "63μm",
    ("O", 1): "146μm",
}


def line_table_label(cname: str, col: int) -> str:
    """Full label for the Line-tables tab, e.g. '[CII] 158μm' or 'CO(1-0)'."""
    return _LINE_LABELS_FULL.get((cname, col), f"{cname}({transition_label(col)})")


def transition_short(cname: str, col: int) -> str:
    """Short per-panel transition label for the Tr/τ tabs, e.g. '158μm' or '1-0'."""
    return _TRANS_SHORT.get((cname, col), transition_label(col))


# ---------------------------------------------------------------------------
# paramsRT.dat — generated from params.dat (the user never edits it)
# ---------------------------------------------------------------------------
def write_paramsrt(directory: str, prefix: str, vturb: float, redshift: float) -> None:
    """Write the 4-line paramsRT.dat (directory, prefix, linewidth km/s, z)."""
    text = f"{directory}\n{prefix}\n{vturb:g}\n{redshift:g}\n"
    paths.paramsrt_dat().write_text(text)


# ---------------------------------------------------------------------------
# stdout summary-table parsing
# ---------------------------------------------------------------------------
# e.g.  "CO( 1- 0) :  1.1527120E+02  1.0602129E+01  3.9039114E+02  7.2310548E+00"
_TRANS_RE = re.compile(
    r"^(?P<c>\S+)\(\s*(?P<u>\d+)\s*-\s*(?P<l>\d+)\)\s*:\s+"
    r"(?P<freq>[-\d.Ee+]+)\s+(?P<tex>[-\d.Ee+]+)\s+(?P<tau>[-\d.Ee+]+)\s+(?P<tr>[-\d.Ee+]+)"
)


@dataclass
class Line:
    transition: str   # "1-0"
    freq_ghz: float
    tex: float
    tau: float
    tr: float


def _ladder_key(line: "Line"):
    """Sort key that prefers adjacent ladder transitions (u = l+1), ordered by
    lower level. RT-tool emits *every* level pair (1-0, 2-0, 2-1, …); the bright
    diagnostic lines are the adjacent ones, so those come first."""
    u, l = (int(t) for t in line.transition.split("-"))
    return (0 if u == l + 1 else 1, l, line.freq_ghz)


def parse_tables(stdout_text: str) -> dict[str, list[Line]]:
    """Parse RT-tool's stdout into ``{coolant: [Line, ...]}``, ordered by ladder
    transition (1-0, 2-1, …) and truncated to the per-coolant display limit.

    We prefer adjacent transitions ordered by lower level (not frequency) so the
    kept lines are the lowest rotational/fine-structure ones — e.g. C⁺ keeps the
    1-0 [CII] 158μm line, and C/O keep 1-0 and 2-1, not the 2-0 overtone.
    """
    grouped: dict[str, list[Line]] = {}
    for raw in stdout_text.splitlines():
        m = _TRANS_RE.match(raw.strip())
        if not m:
            continue
        cname = m.group("c")
        grouped.setdefault(cname, []).append(Line(
            transition=f"{int(m.group('u'))}-{int(m.group('l'))}",
            freq_ghz=fortread.to_float(m.group("freq")),
            tex=fortread.to_float(m.group("tex")),
            tau=fortread.to_float(m.group("tau")),
            tr=fortread.to_float(m.group("tr")),
        ))
    return {cn: sorted(rows, key=_ladder_key)[: table_limit(cn)]
            for cn, rows in grouped.items()}


def save_tables_log(directory: str, prefix: str, stdout_text: str) -> None:
    """Persist the run's stdout so the tables survive a page reload."""
    (paths.pdr_root() / directory / f"RTtool_{prefix}.log").write_text(stdout_text)


def load_tables_log(directory: str, prefix: str) -> str | None:
    p = paths.pdr_root() / directory / f"RTtool_{prefix}.log"
    return p.read_text() if p.exists() else None


# ---------------------------------------------------------------------------
# Tr / tau data files
# ---------------------------------------------------------------------------
def list_coolants(directory: str, prefix: str) -> list[str]:
    """Coolants that have a ``Tr_<prefix>.<coolant>.dat`` file."""
    d = paths.pdr_root() / directory
    if not d.is_dir():
        return []
    head, tail = f"Tr_{prefix}.", ".dat"
    out = [p.name[len(head):-len(tail)] for p in d.glob(f"Tr_{prefix}.*.dat")]
    return sorted(out)


@dataclass
class TrData:
    Ntot: np.ndarray
    NH2: np.ndarray
    Tgas: np.ndarray
    Av: np.ndarray
    Nmol: np.ndarray
    Tr: np.ndarray        # (npoints, ntrans)
    ntrans: int


def load_tr(directory: str, prefix: str, cname: str) -> TrData:
    path = paths.pdr_root() / directory / f"Tr_{prefix}.{cname}.dat"
    data = fortread.loadtxt(path, ndmin=2)  # header line ('# ...') skipped by comments='#'
    return TrData(
        Ntot=data[:, 0], NH2=data[:, 1], Tgas=data[:, 2], Av=data[:, 3], Nmol=data[:, 4],
        Tr=data[:, 5:], ntrans=data.shape[1] - 5,
    )


def load_tau(directory: str, prefix: str, cname: str) -> np.ndarray:
    """Optical-depth matrix (npoints, ntrans), row-aligned with the Tr file."""
    path = paths.pdr_root() / directory / f"tau_{prefix}.{cname}.dat"
    return fortread.loadtxt(path, ndmin=2)


def load_tex(directory: str, prefix: str, cname: str) -> np.ndarray:
    """Excitation-temperature matrix (npoints, ntrans), row-aligned with Tr/tau."""
    path = paths.pdr_root() / directory / f"Tex_{prefix}.{cname}.dat"
    return fortread.loadtxt(path, ndmin=2)


def transition_freqs(tables: dict[str, list[Line]], cname: str, ntrans: int) -> list[float | None]:
    """Frequency (GHz) for each Tr/tau/Tex column of ``cname`` (0-based, ``(c+1)-(c)``).

    Frequencies don't vary with depth, so they come from the stdout summary
    table (``tables`` = ``parse_tables(stdout_text)``) rather than the
    per-depth files. Missing entries (shouldn't normally happen) are ``None``.
    """
    by_col: dict[int, float] = {}
    for r in tables.get(cname, []):
        u, l = (int(t) for t in r.transition.split("-"))
        if u == l + 1:
            by_col[l] = r.freq_ghz
    return [by_col.get(c) for c in range(ntrans)]


def nearest_row(arr: np.ndarray, target: float) -> int:
    """Index of the row in ``arr`` closest to ``target`` (used to pick a depth
    point from a user-chosen column-density value)."""
    return int(np.argmin(np.abs(np.asarray(arr, dtype=float) - target)))
