"""Load 3D-PDR model outputs (``<outdir>/<prefix>.pdr.fin``) for analysis.

The 1D ``.pdr.fin`` layout (see ``writeoutputs.F90``) is, per grid point::

    p, x, Av, Tgas, Tdust, etype, nH, UV, abundance(1..nspec)

i.e. 8 leading columns then one abundance column per species, in the order the
species appear in ``chemfiles/species_<network>.d``. The network is *not* stored
in the file, so we infer it from the number of abundance columns
(REDUCED/MEDIUM/FULL have distinct species counts) — important because an old
output may use a different network than the one currently configured.
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import fortread, network as netmod, paths

# ---------------------------------------------------------------------------
# Heating terms (see src/heatingfunctions.F90, HEATING_RATE 1..12).
# Only the terms actually summed into TOTAL_HEATING are plotted; HEATING_RATE(1)
# (plain photoelectric) and (3) (Weingartner) are commented out of the total in
# the Fortran, so they are excluded here. Index 12 is the total.
# ---------------------------------------------------------------------------
HEATING_LABELS: dict[int, str] = {
    1: "Photoelectric (unused)",
    2: "PAH photoelectric",
    3: "Weingartner (unused)",
    4: "C ionization",
    5: "H₂ formation",
    6: "H₂ photodissociation",
    7: "FUV pumping",
    8: "Cosmic-ray",
    9: "Turbulent",
    10: "Chemical",
    11: "Gas-grain",
}
# Rate indices that contribute to TOTAL_HEATING (in plot order).
HEATING_CONTRIB = [2, 4, 5, 6, 7, 8, 9, 10, 11]
HEATING_TOTAL_IDX = 12

# number of leading (non-abundance) columns in the 1D .pdr.fin
_N_LEAD = 8
_COL = {"p": 0, "x": 1, "Av": 2, "Tgas": 3, "Tdust": 4, "etype": 5, "nH": 6, "UV": 7}


def detect_network(nspec: int) -> str | None:
    """Human label for the network whose species count equals ``nspec`` (X-ray
    aware — an X-ray build is reported as e.g. 'REDUCED+Xrays')."""
    return netmod.label_for_suffix(netmod.detect_suffix(nspec))


def list_models(outdir: str) -> list[str]:
    """Prefixes of all ``*.pdr.fin`` files in ``<root>/<outdir>`` (alphabetical)."""
    d = paths.pdr_root() / outdir
    if not d.is_dir():
        return []
    return sorted(p.name[: -len(".pdr.fin")] for p in d.glob("*.pdr.fin"))


def latest_model(outdir: str) -> str | None:
    """Prefix of the most recently written ``.pdr.fin`` (by mtime) — i.e. the
    model most likely just run. None if the directory has no outputs."""
    d = paths.pdr_root() / outdir
    if not d.is_dir():
        return None
    fins = list(d.glob("*.pdr.fin"))
    if not fins:
        return None
    newest = max(fins, key=lambda p: p.stat().st_mtime)
    return newest.name[: -len(".pdr.fin")]


@dataclass
class Model:
    prefix: str
    network: str | None
    nspec: int
    Av: np.ndarray
    x: np.ndarray
    Tgas: np.ndarray
    Tdust: np.ndarray
    nH: np.ndarray
    UV: np.ndarray
    species: list[str]          # column order
    abundances: np.ndarray      # shape (npoints, nspec)

    def abundance(self, token: str) -> np.ndarray | None:
        """Abundance array for a species token, or None if absent."""
        try:
            j = self.species.index(token)
        except ValueError:
            return None
        return self.abundances[:, j]


def load_pdr(outdir: str, prefix: str) -> Model:
    """Load a model's ``.pdr.fin`` into a :class:`Model`."""
    path = paths.pdr_root() / outdir / f"{prefix}.pdr.fin"
    data = fortread.loadtxt(path, ndmin=2)
    nspec = data.shape[1] - _N_LEAD
    suffix = netmod.detect_suffix(nspec)
    names = netmod.species_list(suffix) if suffix else []
    if len(names) != nspec:
        # fall back to generic names so the file still loads/plots
        names = [f"sp{j+1}" for j in range(nspec)]
    return Model(
        prefix=prefix,
        network=netmod.label_for_suffix(suffix),
        nspec=nspec,
        Av=data[:, _COL["Av"]],
        x=data[:, _COL["x"]],
        Tgas=data[:, _COL["Tgas"]],
        Tdust=data[:, _COL["Tdust"]],
        nH=data[:, _COL["nH"]],
        UV=data[:, _COL["UV"]],
        species=names,
        abundances=data[:, _N_LEAD:_N_LEAD + nspec],
    )


def is_uniform_density(nH: np.ndarray, rtol: float = 1e-3) -> bool:
    """True if the cloud has (essentially) constant number density, i.e. every
    nH value is within ``rtol`` of the mean. Uniform-density clouds have no
    meaningful nH profile to plot against."""
    nH = np.asarray(nH, dtype=float)
    if nH.size == 0:
        return True
    mean = float(np.mean(nH))
    if mean <= 0:
        return bool(np.ptp(nH) == 0)
    return bool(np.ptp(nH) <= rtol * mean)


_PC_CM = 3.085677581e18  # 1 parsec in cm


def column_density(x_pc: np.ndarray, nH: np.ndarray) -> np.ndarray:
    """Cumulative line-of-sight column density N(<x) [cm⁻²] = ∫ n_H dl, with the
    position x (column 1 of .pdr.fin, in pc) converted to cm and n_H (column 6).
    Same length/order as the inputs (the first point is 0)."""
    x = np.asarray(x_pc, dtype=float) * _PC_CM
    nh = np.asarray(nH, dtype=float)
    if x.size < 2:
        return np.zeros_like(nh)
    dN = 0.5 * (nh[1:] + nh[:-1]) * np.diff(x)
    return np.concatenate([[0.0], np.cumsum(dN)])


def _crossing(Av, x_HI, x_H2):
    """Locate the first HI/H2 crossing along depth (Av ascending).

    Returns ``(order, i, t)`` where ``order`` is the Av-sorting permutation,
    ``i`` the lower bracketing index and ``t`` the interpolation fraction into
    ``[i, i+1]``; or None if the curves never cross.
    """
    if Av is None or x_HI is None or x_H2 is None:
        return None
    Av = np.asarray(Av, dtype=float)
    d = np.asarray(x_HI, dtype=float) - np.asarray(x_H2, dtype=float)
    order = np.argsort(Av)
    ds = d[order]
    for i in range(len(ds) - 1):
        if ds[i] == 0.0:
            return order, i, 0.0
        if ds[i] * ds[i + 1] < 0.0:
            denom = ds[i + 1] - ds[i]
            t = 0.0 if denom == 0.0 else -ds[i] / denom
            return order, i, t
    return None


def value_at_crossing(arr: np.ndarray, Av: np.ndarray,
                      x_HI: np.ndarray, x_H2: np.ndarray) -> float | None:
    """Interpolated value of ``arr`` at the HI/H2 crossing (same depth point as
    the transition). Lets the transition be marked on any x-axis (Av or nH)."""
    c = _crossing(Av, x_HI, x_H2)
    if c is None or arr is None:
        return None
    order, i, t = c
    a = np.asarray(arr, dtype=float)[order]
    return float(a[i] + t * (a[i + 1] - a[i]))


def transition_av(Av: np.ndarray, x_HI: np.ndarray, x_H2: np.ndarray) -> float | None:
    """Av where the HI and H2 abundance curves cross (first crossing, by depth),
    linearly interpolated, or None if they never cross."""
    return value_at_crossing(Av, Av, x_HI, x_H2)


def _crossing_2h2(Av, x_HI, x_H2):
    """First depth where x(H) = 2·x(H₂) — the 50 % H₂ molecular-fraction point."""
    if Av is None or x_HI is None or x_H2 is None:
        return None
    Av = np.asarray(Av, dtype=float)
    d = np.asarray(x_HI, dtype=float) - 2.0 * np.asarray(x_H2, dtype=float)
    order = np.argsort(Av)
    ds = d[order]
    for i in range(len(ds) - 1):
        if ds[i] == 0.0:
            return order, i, 0.0
        if ds[i] * ds[i + 1] < 0.0:
            denom = ds[i + 1] - ds[i]
            t = 0.0 if denom == 0.0 else -ds[i] / denom
            return order, i, t
    return None


def value_at_crossing_2h2(arr: np.ndarray, Av: np.ndarray,
                           x_HI: np.ndarray, x_H2: np.ndarray) -> float | None:
    """Value of ``arr`` at the x(H) = 2·x(H₂) crossing (50 % H₂ fraction)."""
    c = _crossing_2h2(Av, x_HI, x_H2)
    if c is None or arr is None:
        return None
    order, i, t = c
    a = np.asarray(arr, dtype=float)[order]
    return float(a[i] + t * (a[i + 1] - a[i]))


def transition_av_2h2(Av: np.ndarray, x_HI: np.ndarray, x_H2: np.ndarray) -> float | None:
    """Av where x(H) = 2·x(H₂) — the 50 % H₂ molecular-fraction point."""
    return value_at_crossing_2h2(Av, Av, x_HI, x_H2)


# ===========================================================================
# Analysis II — heating, cooling, emissivities, level populations
# ===========================================================================
# All of these 1D output files share the leading columns  p, x, Av  then the
# payload (see writeoutputs.F90). They have the same number of rows, in the same
# grid order, as the .pdr.fin file, so the x-axis (Av/nH) and the HI/H2
# transition are taken from the Model and aligned by row index.
_N_LEAD3 = 3  # leading columns p, x, Av in heat/cool/line/spop files


def _path(outdir: str, name: str) -> Path:
    return paths.pdr_root() / outdir / name


@functools.lru_cache(maxsize=128)
def coolant_display_name(datfile: str) -> str:
    """Molecule name from a LAMDA coolant file (line 2), e.g. '12co.dat' -> 'CO'.
    Falls back to the file stem."""
    p = paths.chemfiles_dir() / datfile
    if p.exists():
        lines = p.read_text().splitlines()
        if len(lines) >= 2 and lines[1].strip():
            return lines[1].strip()
    return Path(datfile).stem


def rtspop_coolants(outdir: str, prefix: str) -> list[str]:
    """Coolant .dat basenames in the order the model used them, read from the
    ``<prefix>.RTspop.fin`` header (…coolfiles…, ENDCOOLFILES). The header may or
    may not start with a network line (custom networks omit it), so we simply keep
    every line that names a ``.dat`` coolant file and ignore any other header line."""
    p = _path(outdir, f"{prefix}.RTspop.fin")
    if not p.exists():
        return []
    names: list[str] = []
    with p.open() as fh:
        for line in fh:
            tok = line.strip()
            if tok == "ENDCOOLFILES":
                break
            if tok.endswith(".dat"):
                names.append(Path(tok).name)  # strip 'chemfiles/' prefix
    return names


def load_heating(outdir: str, prefix: str):
    """Return ``(labels, series, total)`` for the heating terms that contribute
    to TOTAL_HEATING. ``series`` maps label -> array (erg s⁻¹ cm⁻³)."""
    data = fortread.loadtxt(_path(outdir, f"{prefix}.heat.fin"), ndmin=2)
    series: dict[str, np.ndarray] = {}
    for k in HEATING_CONTRIB:
        col = _N_LEAD3 + (k - 1)
        if col < data.shape[1]:
            series[HEATING_LABELS[k]] = data[:, col]
    total = data[:, _N_LEAD3 + (HEATING_TOTAL_IDX - 1)]
    return list(series.keys()), series, total


def load_cooling(outdir: str, prefix: str):
    """Return ``(labels, series, total)`` for the per-coolant cooling rates.
    ``series`` maps coolant name -> array (erg s⁻¹ cm⁻³); total is the last
    column (totalcooling)."""
    data = fortread.loadtxt(_path(outdir, f"{prefix}.cool.fin"), ndmin=2)
    ncoo = data.shape[1] - _N_LEAD3 - 1  # last column is the total
    cool_files = rtspop_coolants(outdir, prefix)
    # Name each cooling column from its coolant file (line 2); number any extra.
    names = [coolant_display_name(f) for f in cool_files]
    labels = [names[i] if i < len(names) else f"coolant {i + 1}" for i in range(ncoo)]
    series = {labels[i]: data[:, _N_LEAD3 + i] for i in range(ncoo)}
    total = data[:, -1]
    return labels, series, total


def _coolants_with_ext(outdir: str, prefix: str, ext: str) -> list[str]:
    """Coolant cnames that have a ``<prefix>.<cname>.<ext>`` file."""
    d = paths.pdr_root() / outdir
    if not d.is_dir():
        return []
    head, tail = f"{prefix}.", f".{ext}"
    out = []
    for p in d.glob(f"{prefix}.*.{ext}"):
        name = p.name
        if name.startswith(head) and name.endswith(tail):
            out.append(name[len(head): -len(tail)])
    return sorted(out)


def list_line_coolants(outdir: str, prefix: str) -> list[str]:
    return _coolants_with_ext(outdir, prefix, "line.fin")


def list_spop_coolants(outdir: str, prefix: str) -> list[str]:
    return _coolants_with_ext(outdir, prefix, "spop.fin")


def _ncols(path: Path) -> int:
    with path.open() as fh:
        first = fh.readline()
    return len(first.split())


def line_ntrans(outdir: str, prefix: str, cname: str) -> int:
    """Number of radiative transitions in a coolant's .line.fin (= nlev-1)."""
    p = _path(outdir, f"{prefix}.{cname}.line.fin")
    return max(0, _ncols(p) - _N_LEAD3) if p.exists() else 0


def spop_nlev(outdir: str, prefix: str, cname: str) -> int:
    """Number of levels in a coolant's .spop.fin."""
    p = _path(outdir, f"{prefix}.{cname}.spop.fin")
    return max(0, _ncols(p) - _N_LEAD3) if p.exists() else 0


def load_line(outdir: str, prefix: str, cname: str) -> np.ndarray:
    """Emissivity matrix (npoints, ntrans). Column t (0-based) is the
    (t+2)->(t+1) transition, in erg s⁻¹ cm⁻³."""
    data = fortread.loadtxt(_path(outdir, f"{prefix}.{cname}.line.fin"), ndmin=2)
    return data[:, _N_LEAD3:]


def load_spop(outdir: str, prefix: str, cname: str) -> np.ndarray:
    """Level-population matrix (npoints, nlev). Column l (0-based) is level l+1."""
    data = fortread.loadtxt(_path(outdir, f"{prefix}.{cname}.spop.fin"), ndmin=2)
    return data[:, _N_LEAD3:]
