"""Read the ``<prefix>.chemanalysis.fin`` files written by ``src/analyse_chem.F90``.

The Fortran routine emits, for every grid point and every species in the
network, the dominant formation and destruction pathways (each summing to
~100%). The file is compact and tag-prefixed (see analyse_chem.F90 for the
authoritative description)::

    # 3D-PDR chemanalysis v1
    # nspec=<N> nreac=<M> time=<t> yr
    R <id> <reactant tokens ...> > <product tokens ...>   (legend, once)
    P <point> <x_pc> <Av> <nH> <Tgas> <FUV> <CR_s-1>      (per grid point)
    S <id> <species> <abundance> <Ptot_s-1> <Dtot_s-1>    (per species)
    <reaction_id> <signed_percent>   (+ = formation, - = destruction)

A reaction is quoted by its 1-based ``id`` into the legend, so the bulky
equations appear only once. Numeric fields use Fortran ``ES`` formatting, which
drops the ``E`` for 3-digit exponents (e.g. ``1.08-164``), so we route the whole
file through :func:`fortread.fix_fortran_exponents` before parsing.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

from . import fortread, paths

# header keys may be written with a space after '=' (Fortran ES adds one), so
# match "key= value" tolerantly.
_H_NSPEC = re.compile(r"nspec=\s*(\d+)")
_H_NREAC = re.compile(r"nreac=\s*(\d+)")
_H_TIME = re.compile(r"time=\s*([0-9.eEdD+\-]+)")


@dataclass
class SpeciesData:
    """Pathway analysis for one species at one grid point."""
    abundance: float
    ptot: float                       # total formation rate [s-1]
    dtot: float                       # total destruction rate [s-1]
    formation: list[tuple[int, int]] = field(default_factory=list)    # (reaction id, %)
    destruction: list[tuple[int, int]] = field(default_factory=list)  # (reaction id, %)


@dataclass
class Point:
    """One grid point: its physical state and per-species pathway analysis."""
    index: int
    x_pc: float
    Av: float
    nH: float
    Tgas: float
    FUV: float
    CR: float
    species: dict[str, SpeciesData] = field(default_factory=dict)


@dataclass
class ChemAnalysis:
    prefix: str
    nspec: int
    nreac: int
    time: float
    reactions: dict[int, str]         # reaction id -> "A + B → C + D"
    species_list: list[str]           # species in network (file) order
    points: list[Point]

    # -- convenience accessors ------------------------------------------------
    def x_axis(self, kind: str) -> np.ndarray:
        """Per-point x-values: ``kind`` is 'Av' or 'nH'."""
        attr = "nH" if kind == "nH" else "Av"
        return np.array([getattr(p, attr) for p in self.points], dtype=float)

    def abundance(self, species: str) -> np.ndarray:
        """Abundance of ``species`` over all points (NaN where absent)."""
        return np.array(
            [p.species[species].abundance if species in p.species else np.nan
             for p in self.points],
            dtype=float,
        )

    def is_uniform_density(self, rtol: float = 1e-3) -> bool:
        nH = self.x_axis("nH")
        if nH.size == 0:
            return True
        mean = float(np.mean(nH))
        if mean <= 0:
            return bool(np.ptp(nH) == 0)
        return bool(np.ptp(nH) <= rtol * mean)


# ---------------------------------------------------------------------------
# discovery
# ---------------------------------------------------------------------------
_EXT = ".chemanalysis.fin"


def list_models(outdir: str) -> list[str]:
    """Prefixes of all ``*.chemanalysis.fin`` files in ``<root>/<outdir>``."""
    d = paths.pdr_root() / outdir
    if not d.is_dir():
        return []
    return sorted(p.name[: -len(_EXT)] for p in d.glob(f"*{_EXT}"))


def latest_model(outdir: str) -> str | None:
    """Prefix of the most recently written ``.chemanalysis.fin`` (by mtime)."""
    d = paths.pdr_root() / outdir
    if not d.is_dir():
        return None
    fins = list(d.glob(f"*{_EXT}"))
    if not fins:
        return None
    return max(fins, key=lambda p: p.stat().st_mtime).name[: -len(_EXT)]


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------
def _equation(tokens: list[str]) -> str:
    """Build 'A + B → C + D' from a legend line's tokens after the id."""
    if ">" in tokens:
        gt = tokens.index(">")
        reactants, products = tokens[:gt], tokens[gt + 1:]
    else:                       # tolerate a missing separator
        reactants, products = tokens, []
    lhs = " + ".join(reactants) if reactants else "?"
    rhs = " + ".join(products) if products else "?"
    return f"{lhs} → {rhs}"


def load(outdir: str, prefix: str) -> ChemAnalysis:
    """Parse ``<outdir>/<prefix>.chemanalysis.fin`` into a :class:`ChemAnalysis`."""
    path = paths.pdr_root() / outdir / f"{prefix}{_EXT}"
    text = fortread.fix_fortran_exponents(path.read_text())

    nspec = nreac = 0
    time = 0.0
    reactions: dict[int, str] = {}
    points: list[Point] = []
    species_order: list[str] = []
    seen_species: set[str] = set()

    cur_point: Point | None = None
    cur_species: SpeciesData | None = None

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue

        if line.startswith("#"):
            # parse the metadata header: "# nspec=.. nreac=.. time=.. yr"
            if (mm := _H_NSPEC.search(line)):
                nspec = int(mm.group(1))
            if (mm := _H_NREAC.search(line)):
                nreac = int(mm.group(1))
            if (mm := _H_TIME.search(line)):
                try:
                    time = float(mm.group(1).replace("D", "E").replace("d", "e"))
                except ValueError:
                    pass
            continue

        tok = line.split()
        tag = tok[0]

        if tag == "R":
            rid = int(tok[1])
            reactions[rid] = _equation(tok[2:])

        elif tag == "P":
            cur_point = Point(
                index=int(tok[1]),
                x_pc=float(tok[2]), Av=float(tok[3]), nH=float(tok[4]),
                Tgas=float(tok[5]), FUV=float(tok[6]), CR=float(tok[7]),
            )
            points.append(cur_point)
            cur_species = None

        elif tag == "S":
            name = tok[2]
            cur_species = SpeciesData(
                abundance=float(tok[3]), ptot=float(tok[4]), dtot=float(tok[5]),
            )
            if cur_point is not None:
                cur_point.species[name] = cur_species
            if name not in seen_species:
                seen_species.add(name)
                species_order.append(name)

        else:
            # reaction contribution: "<id> <signed_percent>"
            if cur_species is None:
                continue
            rid, pct = int(tok[0]), int(tok[1])
            if pct >= 0:
                cur_species.formation.append((rid, pct))
            else:
                cur_species.destruction.append((rid, -pct))

    return ChemAnalysis(
        prefix=prefix, nspec=nspec or len(species_order), nreac=nreac or len(reactions),
        time=time, reactions=reactions, species_list=species_order, points=points,
    )


# ---------------------------------------------------------------------------
# Pathway evidence (deterministic digest for the LLM summary)
# ---------------------------------------------------------------------------
def pathway_evidence(ca: ChemAnalysis, species: str, top_n: int = 3) -> dict:
    """Condense a species' per-point pathways into the most important parts: the
    formation/destruction reactions that dominate across the column (with the Av
    range and typical contribution), plus a few depth snapshots (edge / abundance
    peak / interior). This is what the LLM summarises — the heavy lifting (which
    reactions matter and where) is done here, deterministically, so the model
    only has to phrase it."""
    pts = sorted(ca.points, key=lambda p: p.Av)
    if not pts or species not in ca.species_list:
        return {}

    def eqn(rid: int) -> str:
        return ca.reactions.get(rid, f"reaction {rid}")

    av = [p.Av for p in pts]
    nH = [p.nH for p in pts]
    Tg = [p.Tgas for p in pts]
    sds = [p.species.get(species) for p in pts]
    ab = [sd.abundance if sd else 0.0 for sd in sds]

    def dominant(kind: str):
        out = []
        for sd in sds:
            lst = getattr(sd, kind) if sd else []
            out.append(max(lst, key=lambda t: t[1]) if lst else None)
        return out

    def top_reactions(kind: str):
        seq = dominant(kind)
        present = [i for i, d in enumerate(seq) if d is not None]
        agg: dict[int, dict] = {}
        for i in present:
            rid, pct = seq[i]
            a = agg.setdefault(rid, {"count": 0, "pcts": [], "avs": []})
            a["count"] += 1
            a["pcts"].append(pct)
            a["avs"].append(av[i])
        n = max(1, len(present))
        items = sorted(agg.items(), key=lambda kv: -kv[1]["count"])[:top_n]
        return [{
            "reaction": eqn(rid),
            "dominant_fraction": round(a["count"] / n, 2),
            "typical_pct": int(sorted(a["pcts"])[len(a["pcts"]) // 2]),
            "av_lo": min(a["avs"]), "av_hi": max(a["avs"]),
        } for rid, a in items]

    def top2(lst):
        return [{"reaction": eqn(r), "pct": pc}
                for r, pc in sorted(lst, key=lambda t: -t[1])[:2]] if lst else []

    def snapshot(i: int, label: str):
        sd = sds[i]
        return {"label": label, "Av": av[i], "nH": nH[i], "Tgas": Tg[i],
                "abundance": ab[i],
                "formation": top2(sd.formation) if sd else [],
                "destruction": top2(sd.destruction) if sd else []}

    aba = np.asarray(ab, dtype=float)
    peak_i = int(np.argmax(aba)) if aba.size and float(np.max(aba)) > 0 else 0
    snaps, seen = [], set()
    for i, label in [(0, "cloud edge (low Av)"), (peak_i, "abundance peak"),
                     (len(pts) - 1, "cloud interior (high Av)")]:
        if i not in seen:
            seen.add(i)
            snaps.append(snapshot(i, label))

    nHmin, nHmax = min(nH), max(nH)
    return {
        "species": species,
        "npoints": len(pts),
        "av_range": (av[0], av[-1]),
        "nH_range": (nHmin, nHmax),
        "Tgas_range": (min(Tg), max(Tg)),
        "uniform_density": (nHmax - nHmin) <= 1e-3 * nHmax if nHmax > 0 else True,
        "abundance": {"edge": ab[0], "deep": ab[-1],
                      "max": float(aba.max()) if aba.size else 0.0,
                      "peak_Av": av[peak_i]},
        "formation_top": top_reactions("formation"),
        "destruction_top": top_reactions("destruction"),
        "snapshots": snaps,
    }
