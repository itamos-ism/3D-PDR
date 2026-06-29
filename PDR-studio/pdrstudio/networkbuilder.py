"""Build a new chemical network for 3D-PDR with MakeRates, then import it.

Workflow (see proposal/makerates.txt):

1. The user picks the elements to include (e-, H, H2, He, C+, O are always in).
2. From an imported master reaction file we collect the species whose elemental
   constituents are a subset of the selected elements and write a species file.
   (Alternatively the user imports their own species file.)
3. MakeRates (embedded, run headless) turns reaction file + species file into
   ``species_<NAME>.d`` / ``rates_<NAME>.d`` / ``odes_<NAME>.c`` (all initial
   abundances zeroed in the species file).
4. On request these are imported into the 3D-PDR tree: species/rates ->
   ``chemfiles/``, odes -> ``src/``, and ``config.mk`` / ``makefile`` /
   the suffix ``#ifdef`` blocks in ``read_species.F90`` / ``read_rates.F90`` /
   ``input_parameters.F90`` are updated to add ``NETWORK = <NAME>``.

Chemical heating needs no per-network edit thanks to the ``AUTOCHEMHEAT`` flag,
so ``heatingfunctions.F90`` is left untouched.

Note: MakeRates emits the ODEs as C (``odes_<NAME>.c``); 3D-PDR compiles them via
the makefile's ``%.o: %.c`` rule, so the file is ``.c`` (not ``.d``).
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from . import paths

# ---------------------------------------------------------------------------
# Elements
# ---------------------------------------------------------------------------
# Always included (cannot be deselected). Map UI species -> element symbol.
MANDATORY_ELEMENTS = ["H", "He", "C", "O"]      # e- handled via charge
# Optional, user-selectable. Label shown -> element symbol.
OPTIONAL_ELEMENTS: dict[str, str] = {
    "Mg+": "Mg", "N": "N", "S": "S", "Fe": "Fe", "Cl": "Cl",
    "P": "P", "Si": "Si", "Na": "Na", "F": "F",
}

# Element symbols recognised in species formulae, longest/ambiguous first so e.g.
# "He" is matched before "H" (mirrors MakeRates' find_constituents element list).
_ELEMENTS = ['PAH', 'He', 'Li', 'Na', 'Mg', 'Si', 'Cl', 'Ca', 'Fe',
             'H', 'D', 'C', 'N', 'O', 'F', 'P', 'S', '#', '+', '-']
# Tokens that are not real species (pseudo-reactants / markers).
_PSEUDO = {'', '#', 'e-', 'ELECTR', 'PHOTON', 'CRP', 'CRPHOT', 'XRAY', 'XRSEC',
           'XRLYA', 'XRPHOT', 'FREEZE', 'CRH', 'PHOTD', 'THERM'}


def _parse(species: str) -> tuple[dict[str, int], str]:
    """Return (element->count, residue). ``residue`` is whatever alphabetic text
    is left once every recognised element/digit/charge has been consumed — non-
    empty means the species contains an element we don't model (Ar, Ti, Al, …)."""
    s = species
    out: dict[str, int] = {}
    for element in _ELEMENTS:
        n = 0
        for _ in range(s.count(element)):
            index = s.index(element) + len(element)
            if s[index:index + 2].isdigit():
                n += int(s[index:index + 2]); s = s[:index - len(element)] + s[index + 2:]
            elif s[index:index + 1].isdigit():
                n += int(s[index:index + 1]); s = s[:index - len(element)] + s[index + 1:]
            else:
                n += 1; s = s[:index - len(element)] + s[index:]
        if n:
            out[element] = n
    residue = "".join(ch for ch in s if ch.isalpha())
    return out, residue


def constituents(species: str) -> dict[str, int]:
    """Element -> count for a species formula (same algorithm as MakeRates'
    ``find_constituents``). Charge ('+'/'-') and grain ('#') are returned too."""
    return _parse(species)[0]


def has_unknown_element(species: str) -> bool:
    """True if the species contains an element we don't model (Ar, Ti, Al, …)."""
    return bool(_parse(species)[1])


def total_atoms(species: str) -> int:
    """Total number of atoms (incl. H) in a species; charge/grain excluded."""
    return sum(c for el, c in constituents(species).items() if el not in ("+", "-", "#"))


def species_mass(species: str) -> float:
    mass = {'PAH': 420.0, 'He': 4.0, 'Li': 7.0, 'Na': 23.0, 'Mg': 24.0, 'Si': 28.0,
            'Cl': 35.0, 'Ca': 40.0, 'Fe': 56.0, 'H': 1.0, 'D': 2.0, 'C': 12.0,
            'N': 14.0, 'O': 16.0, 'F': 19.0, 'P': 31.0, 'S': 32.0}
    return sum(mass.get(el, 0.0) * cnt for el, cnt in constituents(species).items())


# ---------------------------------------------------------------------------
# Reaction file -> species (Rate05 CSV master networks)
# ---------------------------------------------------------------------------
def species_in_reaction_file(text: str) -> list[str]:
    """All real species appearing in a Rate05 (CSV) master reaction file, in
    order of first appearance. Reactant slots are fields 1-3, products 4-7."""
    seen: list[str] = []
    known: set[str] = set()
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith(("#", "!", "C ")):
            continue
        f = [t.strip() for t in line.split(",")]
        if len(f) < 11 or not f[0].lstrip("-").isdigit():
            continue
        for tok in f[1:8]:
            if tok and tok not in _PSEUDO and tok not in known:
                known.add(tok); seen.append(tok)
    return seen


def is_rate05_reaction_file(text: str) -> bool:
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith(("#", "!", "C ")):
            continue
        return len(s.split(",")) > 10
    return False


def heavy_atom_count(species: str) -> int:
    """Number of non-hydrogen atoms in a species (charge markers excluded)."""
    return sum(c for el, c in constituents(species).items()
               if el not in ("+", "-", "H", "D"))


def select_species_by_elements(all_species: list[str], elements: set[str],
                               max_heavy: int | None = None,
                               max_carbon: int | None = None,
                               include_anions: bool = False) -> list[str]:
    """Species whose chemical constituents are all within ``elements`` (cations
    always allowed), optionally limited in complexity. Species containing an
    element we don't model (Ar, Ti, Al, …) are always excluded.

    * ``max_carbon`` — keep only species with at most this many C atoms (the main
      lever against carbon chains, e.g. set 1 for a reduced-style network);
    * ``max_heavy``  — keep only species with at most this many heavy (non-H)
      atoms, excluding long/heavy molecules such as C4;
    * ``include_anions`` — when False, drop negatively-charged species (e.g. C-).

    'H','He','C','O' are assumed to be in ``elements``."""
    out = []
    for sp in all_species:
        cons, residue = _parse(sp)
        if residue:                          # unmodelled element (Ar, Ti, Al, …)
            continue
        if any(el not in ("+", "-") and el not in elements for el in cons):
            continue  # contains an element that wasn't selected (incl '#','PAH','D',…)
        if not include_anions and cons.get("-", 0) > 0:
            continue
        if max_carbon is not None and cons.get("C", 0) > max_carbon:
            continue
        if max_heavy is not None and heavy_atom_count(sp) > max_heavy:
            continue
        out.append(sp)
    return out


def write_species_file(path: Path, species: list[str]) -> None:
    """Write a Rate05 species file (``index,SPECIES,0.00E+00,mass``) for MakeRates."""
    lines = [f"{i + 1},{sp},0.00E+00,{species_mass(sp):.1f}"
             for i, sp in enumerate(species)]
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# MakeRates (embedded, headless)
# ---------------------------------------------------------------------------
def makerates_script() -> Path | None:
    """Locate the MakeRates script: the embedded copy first (portable), else the
    original checkout."""
    embedded = Path(__file__).resolve().parent.parent / "makerates" / "MakeRates.py"
    if embedded.exists():
        return embedded
    external = Path.home() / "data/Codes/MakeRates/OriginalMakeRates/MakeRates_python312.py"
    return external if external.exists() else None


@dataclass
class BuildResult:
    name: str
    species_file: Path     # species_<NAME>.d (abundances zeroed)
    rates_file: Path       # rates_<NAME>.d
    odes_file: Path        # odes_<NAME>.c
    n_species: int
    n_reactions: int
    log: str


def _zero_abundances(species_dat: Path, out: Path) -> int:
    """Copy a Rate05 species file forcing every abundance to 0.00E+00. Returns
    the number of species lines (excluding the trailing e-)."""
    rows = []
    n = 0
    for line in species_dat.read_text().splitlines():
        f = line.split(",")
        if len(f) >= 4:
            f[2] = "0.00E+00"
            rows.append(",".join(f))
            if f[1].strip().lower() != "e-":
                n += 1
        elif line.strip():
            rows.append(line)
    out.write_text("\n".join(rows) + "\n")
    return n


def staging_dir(name: str) -> Path:
    d = Path(__file__).resolve().parent.parent / ".network_builds" / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_makerates(reaction_file: Path, species_file: Path, name: str) -> BuildResult:
    """Run MakeRates headless and produce species_<NAME>.d / rates_<NAME>.d /
    odes_<NAME>.c in a per-network staging directory."""
    script = makerates_script()
    if script is None:
        raise RuntimeError("MakeRates script not found (expected at PDR-studio/makerates/MakeRates.py).")

    work = staging_dir(name)
    prefix = "nb_"
    cmd = [sys.executable, str(script),
           f"reactionFile={reaction_file.resolve()}",
           f"speciesFile={species_file.resolve()}",
           f"outputPrefix={prefix}", "codeFormat=C", "sortSpecies=True"]
    proc = subprocess.run(cmd, cwd=str(work), capture_output=True, text=True, timeout=1200)
    log = (proc.stdout or "") + (proc.stderr or "")

    sp_raw = work / f"{prefix}species.dat"
    rt_raw = work / f"{prefix}rates.dat"
    od_raw = work / f"{prefix}odes.c"
    if not (sp_raw.exists() and rt_raw.exists() and od_raw.exists()):
        raise RuntimeError(f"MakeRates did not produce the expected outputs.\n\n{log}")

    species_out = work / f"species_{name}.d"
    rates_out = work / f"rates_{name}.d"
    odes_out = work / f"odes_{name}.c"
    n_species = _zero_abundances(sp_raw, species_out)
    shutil.copyfile(rt_raw, rates_out)
    shutil.copyfile(od_raw, odes_out)
    n_reactions = sum(1 for ln in rates_out.read_text().splitlines() if ln.strip())

    return BuildResult(name=name, species_file=species_out, rates_file=rates_out,
                       odes_file=odes_out, n_species=n_species,
                       n_reactions=n_reactions, log=log)


# ---------------------------------------------------------------------------
# Import into the 3D-PDR tree
# ---------------------------------------------------------------------------
_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def valid_name(name: str) -> bool:
    """A network name must be a valid C-preprocessor macro / file token."""
    return bool(_NAME_RE.match(name or ""))


def _patch_suffix_fortran(path: Path, name: str) -> bool:
    """Add a ``#elif <NAME>`` branch to a file's suffix ``#ifdef`` chain, right
    after the last existing branch. Idempotent. Returns True if changed."""
    if not path.exists():
        return False
    text = path.read_text()
    if re.search(rf"#elif\s+{re.escape(name)}\b", text):
        return False
    # Match the final branch of the suffix chain: "#elif MYNETWORK\n  suffix = 'mynetwork.d'"
    m = re.search(r"(#elif\s+\w+\s*\n([ \t]*)suffix\s*=\s*'[^']*'\s*\n)(\s*#endif)", text)
    if not m:
        return False
    indent = m.group(2)
    branch = f"#elif {name}\n{indent}suffix = '{name}.d'\n"
    new = text[:m.end(1)] + branch + text[m.end(1):]
    path.write_text(new)
    return True


def _patch_makefile(path: Path, name: str) -> bool:
    """Add the ``-D<NAME>`` define and ``odes_<NAME>.o`` object to the makefile,
    after the MYNETWORK blocks. Idempotent."""
    text = path.read_text()
    changed = False
    if not re.search(rf"-D{re.escape(name)}\b", text):
        anchor = re.search(r"ifeq\s*\(\s*\$\(NETWORK\),\s*MYNETWORK\s*\)\s*\n\s*CFLAGS\s*\+=\s*-DMYNETWORK\s*\nendif\n", text)
        block = f"ifeq ($(NETWORK),{name})\n  CFLAGS += -D{name}\nendif\n"
        if anchor:
            text = text[:anchor.end()] + block + text[anchor.end():]
            changed = True
    if not re.search(rf"odes_{re.escape(name)}\.o\b", text):
        anchor = re.search(r"ifeq\s*\(\s*\$\(NETWORK\),\s*MYNETWORK\s*\)\s*\nCODE_OBJ\s*\+=\s*odes_mynetwork\.o\s*\nendif\n", text)
        block = f"ifeq ($(NETWORK),{name})\nCODE_OBJ += odes_{name}.o\nendif\n"
        if anchor:
            text = text[:anchor.end()] + block + text[anchor.end():]
            changed = True
    if changed:
        path.write_text(text)
    return changed


@dataclass
class ImportReport:
    copied: list[str]
    patched: list[str]
    skipped: list[str]
    warnings: list[str]


def import_network(result: BuildResult, set_active: bool = True) -> ImportReport:
    """Copy the built files into the 3D-PDR tree and wire up the NETWORK flag."""
    from . import config  # local import to avoid any import cycle

    name = result.name
    root = paths.pdr_root()
    chemdir = paths.chemfiles_dir()
    srcdir = paths.src_dir()
    copied, patched, skipped, warnings = [], [], [], []

    # 1) files
    shutil.copyfile(result.species_file, chemdir / f"species_{name}.d")
    shutil.copyfile(result.rates_file, chemdir / f"rates_{name}.d")
    shutil.copyfile(result.odes_file, srcdir / f"odes_{name}.c")
    copied += [f"chemfiles/species_{name}.d", f"chemfiles/rates_{name}.d", f"src/odes_{name}.c"]

    # 2) makefile
    if _patch_makefile(paths.makefile(), name):
        patched.append("src/makefile")
    else:
        skipped.append("src/makefile (already wired or anchor not found)")

    # 3) Fortran suffix blocks
    for fn in ("read_species.F90", "read_rates.F90", "input_parameters.F90"):
        if _patch_suffix_fortran(srcdir / fn, name):
            patched.append(f"src/{fn}")
        else:
            skipped.append(f"src/{fn} (already has branch or pattern not found)")

    # 4) config.mk -> NETWORK = NAME
    if set_active and paths.config_mk().exists():
        cfg = config.read_config()
        cfg["NETWORK"] = name
        config.write_config(cfg)
        patched.append("src/config.mk (NETWORK)")

    # sanity: warn if any suffix file wasn't patched (would break the build)
    for fn in ("read_species.F90", "read_rates.F90", "input_parameters.F90"):
        if f"src/{fn}" not in patched and "already has branch" not in " ".join(skipped):
            pass  # already reported in skipped
    return ImportReport(copied=copied, patched=patched, skipped=skipped, warnings=warnings)
