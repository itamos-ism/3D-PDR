"""Browse and edit the chemical reaction network in ``chemfiles/rates_*.d``.

The rates files are CSV in the UMIST/Rate05 style. ``src/read_rates.F90`` reads
the first 14 fields of each line::

    N, R1, R2, R3, P1, P2, P3, P4, ALPHA, BETA, GAMMA, CLEM, Tmin, Tmax

Empty reactant/product slots are blank fields. Anything after ``Tmax`` (a
quality letter, a trailing comma, …) is an annotation the Fortran ignores; we
preserve it verbatim. A reaction can appear on several lines with different
temperature ranges (e.g. dielectronic recombination split into T intervals),
so the file *line* is the unique handle for editing, not the reaction index.

Editing rewrites only the affected line, keeping every other byte of the file
unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from . import network as netmod
from . import paths

# Field positions in the CSV record.
_I_N = 0
_I_REAC = (1, 2, 3)
_I_PROD = (4, 5, 6, 7)
_I_ALPHA = 8
_I_BETA = 9
_I_GAMMA = 10
_I_CLEM = 11
_I_TMIN = 12
_I_TMAX = 13
_MIN_FIELDS = 14


def rates_file(network: str) -> Path:
    """Path to the rates file for ``network`` (X-ray variant when XRAYS=1)."""
    return netmod.rates_path(netmod.suffix_for(network))


@dataclass
class Reaction:
    line_no: int               # 0-based line index in the file (the edit handle)
    index: int                 # reaction number N as written in the file
    reactants: list[str]       # non-empty reactant tokens
    products: list[str]        # non-empty product tokens
    alpha: float
    beta: float
    gamma: float
    clem: str
    tmin: float
    tmax: float
    raw: list[str] = field(default_factory=list)  # full original CSV fields

    def equation(self) -> str:
        return f"{' + '.join(self.reactants)} → {' + '.join(self.products)}"

    def label(self) -> str:
        return f"#{self.index}  {self.equation()}  [{self.tmin:g}–{self.tmax:g} K]"


def _to_float(s: str, default: float = 0.0) -> float:
    try:
        return float(s)
    except (ValueError, TypeError):
        return default


def read_reactions(network: str) -> list[Reaction]:
    """Parse the network's rates file into Reaction objects (order preserved)."""
    path = rates_file(network)
    if not path.exists():
        return []
    reactions: list[Reaction] = []
    for line_no, raw_line in enumerate(path.read_text().splitlines()):
        if not raw_line.strip():
            continue
        fields = [f.strip() for f in raw_line.split(",")]
        if len(fields) < _MIN_FIELDS:
            continue
        reactants = [fields[i] for i in _I_REAC if fields[i] and fields[i] != ""]
        products = [fields[i] for i in _I_PROD if fields[i] and fields[i] != ""]
        try:
            index = int(fields[_I_N])
        except ValueError:
            continue
        reactions.append(Reaction(
            line_no=line_no,
            index=index,
            reactants=reactants,
            products=products,
            alpha=_to_float(fields[_I_ALPHA]),
            beta=_to_float(fields[_I_BETA]),
            gamma=_to_float(fields[_I_GAMMA]),
            clem=fields[_I_CLEM],
            tmin=_to_float(fields[_I_TMIN]),
            tmax=_to_float(fields[_I_TMAX]),
            raw=fields,
        ))
    return reactions


def species_tokens(reactions: list[Reaction]) -> list[str]:
    """Sorted unique set of all species tokens appearing in the network
    (reactants and products), for populating search menus."""
    toks: set[str] = set()
    for r in reactions:
        toks.update(r.reactants)
        toks.update(r.products)
    return sorted(toks, key=str.lower)


def nonspecies_reactants(reactions: list[Reaction], species: list[str] | set[str]) -> list[str]:
    """Reactant tokens that are *not* real species — the pseudo-reactants that
    drive non-collisional processes: PHOTON, CRP, CRPHOT, FREEZE, THERM and the
    X-ray ones XRLYA, XRSEC, XRPHOT, etc. ``species`` is the network's species
    list (so the set adapts to whichever network/files are active)."""
    sp = set(species)
    toks: set[str] = set()
    for r in reactions:
        for t in r.reactants:
            if t and t not in sp:
                toks.add(t)
    return sorted(toks, key=str.lower)


def search_reactions(
    reactions: list[Reaction],
    reactants: list[str] | None = None,
    products: list[str] | None = None,
) -> list[Reaction]:
    """Return reactions containing *all* of the requested reactant tokens and
    *all* of the requested product tokens. Empty lists impose no constraint.
    Matching is exact on the token (case-insensitive)."""
    req_r = [s.strip().lower() for s in (reactants or []) if s.strip()]
    req_p = [s.strip().lower() for s in (products or []) if s.strip()]

    def has_all(have: list[str], want: list[str]) -> bool:
        low = [h.lower() for h in have]
        return all(w in low for w in want)

    out = []
    for r in reactions:
        if has_all(r.reactants, req_r) and has_all(r.products, req_p):
            out.append(r)
    return out


def _fmt_alpha(x: float) -> str:
    return f"{x:.2E}"


def _fmt_beta(x: float) -> str:
    return f"{x:.2f}"


def _fmt_gamma(x: float) -> str:
    return f"{x:.1f}"


def _fmt_temp(x: float) -> str:
    """Whole numbers as integers (matches the file), else one decimal."""
    return str(int(x)) if float(x).is_integer() else f"{x:.1f}"


def update_reaction(
    network: str,
    line_no: int,
    *,
    alpha: float,
    beta: float,
    gamma: float,
    tmin: float | None = None,
    tmax: float | None = None,
    clem: str | None = None,
) -> None:
    """Rewrite a single line of the rates file with new coefficients, preserving
    every other line and all trailing annotation fields on the edited line."""
    path = rates_file(network)
    lines = path.read_text().splitlines(keepends=True)
    if not (0 <= line_no < len(lines)):
        raise IndexError(f"line {line_no} out of range for {path.name}")

    raw_line = lines[line_no]
    eol = "\n" if raw_line.endswith("\n") else ""
    fields = raw_line.rstrip("\n").split(",")
    # Preserve leading/trailing spacing minimally by working on stripped copies;
    # the rates files are not whitespace-padded so this is faithful.
    fields[_I_ALPHA] = _fmt_alpha(alpha)
    fields[_I_BETA] = _fmt_beta(beta)
    fields[_I_GAMMA] = _fmt_gamma(gamma)
    if tmin is not None:
        fields[_I_TMIN] = _fmt_temp(tmin)
    if tmax is not None:
        fields[_I_TMAX] = _fmt_temp(tmax)
    if clem is not None:
        fields[_I_CLEM] = clem

    lines[line_no] = ",".join(fields) + eol
    path.write_text("".join(lines))
