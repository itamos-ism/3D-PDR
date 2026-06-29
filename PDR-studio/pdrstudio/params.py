"""Read and write ``params.dat``, the coolant selection, and the per-network
initial elemental abundances in ``chemfiles/species_*.d``.

``params.dat`` (see ``src/input_parameters.F90``) is a *positional* file: header
blocks (``===`` / title / ``===``, possibly decorated with ``|`` columns as in
3D-XDR) introduce each section, then one value per line, each with an inline
``!label``; the final section is a variable-length coolant list.

Because different 3D-PDR / 3D-XDR versions reword labels and *insert* new
parameter lines (e.g. an X-ray flux entry under the cosmic-ray rate), reading is
made version-robust (:func:`_resolve_block`): a section's fields are matched by
their inline label, with a positional fallback when a section's entry count
matches and an "extras" path for genuinely new entries. Writing is
structure-preserving — section order, headers, comments, the coolant block and
any unknown ("extra") entries are kept in place, with only values substituted —
so new entries round-trip safely. A fresh file (none present) is created from
the built-in template.

The cosmic-ray field is the one curated variable-length block: its layout
depends on ``CRATTENUATION`` (0 -> ``zeta``; 1 -> an L/H letter; 2 -> norm/n0/
slope), reconstructed on write only if the file's layout no longer matches.
``Av_crit``/``v_alfv`` are always present but only matter when ``SUPRATHERMAL=1``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from . import network as netmod
from . import paths

# --- separator used for the block headers we emit -------------------------
_SEP = "================================="

# --- default values (used for new files and "reset to defaults") ----------
DEFAULTS: dict[str, object] = {
    # I/O
    "indir": "ics",
    "input": "1Dn30.dat",
    "outdir": "sims",
    "output": "model",
    # PDR
    "fuv": 1.0,            # Draine (1974) units
    "zeta": 1.0e-17,       # CR ionisation rate (CRATTENUATION=0)
    "cr_field": "L",       # CRATTENUATION=1 letter
    "cr_norm": 1.0,        # CRATTENUATION=2
    "cr_n0": 1.0e3,
    "cr_slope": 0.0,
    "d2g": 1.0,            # dust-to-gas (normalised to 1e-2); read into `metallicity`
    "vturb": 1.0,          # km/s
    "tend": 1.0e7,         # yr
    "grain_radius": 1.0e-5,  # cm
    "av_fac": 6.289e-22,
    "uv_fac": 3.02,
    "redshift": 0.0,
    "av_crit": 0.7,
    "v_alfv": 3.3,         # km/s
    # Ray-tracing
    "healpix_level": 0,
    "theta_crit": 1.3,
    # ODE solver
    "rel_tol": 1.0e-8,
    "abs_tol": 1.0e-30,
    # Chemistry iterations
    "init_iter": 8,
    "max_iter": 6000,
    # Thermal balance
    "tgas": 40.0,
    "tmin": 1.0,
    "tmax": 20000.0,
    "tdust": 20.0,
    "fcrit": 0.005,
    "tdiff": 0.01,
}

# Fields that are written in scientific notation.
_SCI_FIELDS = {"zeta", "tend", "grain_radius", "av_fac", "rel_tol", "abs_tol",
               "cr_n0"}

# Inline comments appended to each value line for readability. Fortran's
# list-directed read takes only the first whitespace-delimited token, so the
# "!comment" tail is ignored by the code but helps a human reading the file.
_COMMENTS: dict[str, str] = {
    "indir": "ICs directory",
    "input": "ICs file",
    "outdir": "Output directory",
    "output": "Output prefix",
    "fuv": "FUV strength (Draine 1974 units)",
    "zeta": "CR ionisation rate (s^-1)",
    "cr_field": "CR field choice (L/H)",
    "cr_norm": "CR attenuation norm",
    "cr_n0": "CR attenuation n0",
    "cr_slope": "CR attenuation slope",
    "d2g": "Dust-to-gas ratio (norm. to 1e-2)",
    "vturb": "Microturbulent velocity (km/s)",
    "tend": "Chemical evolution time (yr)",
    "grain_radius": "Grain radius (cm)",
    "av_fac": "N(H) -> Av conversion (mag cm^2)",
    "uv_fac": "UV factor",
    "redshift": "Redshift (sets T_CMB)",
    "av_crit": "Av_crit (suprathermal; SUPRATHERMAL=1)",
    "v_alfv": "Alfven velocity km/s (SUPRATHERMAL=1)",
    "healpix_level": "HEALPix level",
    "theta_crit": "theta_crit",
    "rel_tol": "ODE relative tolerance",
    "abs_tol": "ODE absolute tolerance",
    "init_iter": "Initial chemistry iterations",
    "max_iter": "Max iterations",
    "tgas": "Gas temperature (isothermal if THERMALBALANCE=0)",
    "tmin": "Min gas temperature",
    "tmax": "Max gas temperature",
    "tdust": "Dust temperature (used if DUST=0)",
    "fcrit": "Convergence Fcrit",
    "tdiff": "Convergence Tdiff",
}
_COMMENT_COL = 34  # column at which inline comments start


def cr_nlines(cratten: int) -> int:
    """How many lines the cosmic-ray field occupies for a given CRATTENUATION."""
    return 3 if cratten == 2 else 1


# ---------------------------------------------------------------------------
# Coolant catalog
# ---------------------------------------------------------------------------
# These four are always included (the user cannot deselect them).
FIXED_COOLANTS: list[tuple[str, str]] = [
    ("CO", "12co.dat"),
    ("C+", "12c+.dat"),
    ("C", "12c.dat"),
    ("O", "16o.dat"),
]


@dataclass(frozen=True)
class Coolant:
    label: str
    filename: str


_FIXED_FILES = {f for _, f in FIXED_COOLANTS}


def coolant_species(filename: str) -> str:
    """Coolant species name from line 2 of a LAMDA coolant file (e.g. 'HCN')."""
    p = paths.chemfiles_dir() / filename
    if not p.exists():
        return ""
    lines = p.read_text().splitlines()
    return lines[1].strip() if len(lines) >= 2 else ""


def _coolant_label(filename: str, species: str) -> str:
    return f"{species} (hfs)" if "hfs" in filename else species


def optional_for_network(network: str) -> list[Coolant]:
    """Optional coolants whose coolant species is present in ``network``'s species
    list. Discovered dynamically from chemfiles/ so new .dat files are picked up
    automatically. When two files share the same UI label (e.g. hco+.dat and
    hcop.dat both → 'HCO+'), the first in alphabetical order wins."""
    suffix = netmod.suffix_for(network)
    have = set(netmod.species_list(suffix))
    d = paths.chemfiles_dir()
    seen_labels: set[str] = set()
    out: list[Coolant] = []
    for fn in sorted(p.name for p in d.glob("*.dat") if p.name not in _FIXED_FILES):
        sp = coolant_species(fn)
        if not sp or sp not in have:
            continue
        label = _coolant_label(fn, sp)
        if label in seen_labels:
            continue
        seen_labels.add(label)
        out.append(Coolant(label, fn))
    return out


# ---------------------------------------------------------------------------
# Initial conditions — uniform 1-D grid (Python port of ics/uniform1D.f90)
# ---------------------------------------------------------------------------
_ICS_PC_CM = 3.0856e18    # parsec in cm (as in uniform1D.f90)
_ICS_AV0 = 6.289e-22      # N(H)->Av conversion (as in uniform1D.f90)


def uniform_ics_grid(nH: float, av_max: float, log_av_min: float, resolution: int) -> str:
    """Build a 1-D uniform-density initial-conditions grid (port of
    ``ics/uniform1D.f90``). Columns are ``x[pc] y z nH``: an x=0 edge point
    followed by points log-spaced in Av from ``10**log_av_min`` up to ~``av_max``
    with ``resolution`` points per Av dex (uniform density ``nH``)."""
    import math
    rows: list[tuple[float, float, float, float]] = [(0.0, 0.0, 0.0, nH)]
    lam = (math.log10(av_max) - log_av_min) * resolution
    for i in range(int(lam) + 1):
        av = 10.0 ** (log_av_min + i / resolution)
        x = av / _ICS_PC_CM / _ICS_AV0 / nH
        rows.append((x, 0.0, 0.0, nH))
    return "".join("".join(f"{v:15.7E}" for v in r) + "\n" for r in rows)


# ---------------------------------------------------------------------------
# Initial elemental abundances (per network)
# ---------------------------------------------------------------------------
# Canonical default abundances per species token (from the spec). Only those
# tokens that actually exist in a given network's species file are shown/edited,
# so a freshly-built network (Network Builder) shows exactly its element set.
ELEMENT_DEFAULTS: dict[str, float] = {
    "H": 4.00e-1,
    "H2": 3.00e-1,
    "He": 1.00e-1,
    "C+": 1.40e-4,
    "O": 3.00e-4,
    "Mg+": 2.70e-7,
    "N": 5.65e-5,
    "S": 3.51e-6,
    "Fe": 2.01e-7,
    "Cl": 3.39e-8,
    "P": 7.78e-8,
    "Si": 1.78e-6,
}
# Display order.
ELEMENT_ORDER = ["H", "H2", "He", "C+", "O", "Mg+", "N", "S", "Fe", "Cl", "P", "Si"]


def species_file(network: str) -> Path:
    """Species file for ``network`` (X-ray variant when XRAYS=1)."""
    return netmod.species_path(netmod.suffix_for(network))


def _fmt(field_name: str, value: object) -> str:
    """Format one value for params.dat, matching the existing file's style."""
    if isinstance(value, str):
        return value
    if field_name in {"healpix_level", "init_iter", "max_iter"}:
        return str(int(value))
    if field_name in _SCI_FIELDS:
        return f"{float(value):.3E}"
    return f"{float(value):.4g}"


# ---------------------------------------------------------------------------
# params.dat structure: comment-driven read, structure-preserving write
# ---------------------------------------------------------------------------
# params.dat is positional (Fortran reads one value per line, in order), but
# every value carries an inline "!label". We key on that label so a newer code
# version that *inserts* parameter lines (anywhere before the final coolant
# block) does not shift our reads, and so unknown new entries are kept verbatim
# in their original position on save.
_STR_FIELDS = {"indir", "input", "outdir", "output", "cr_field"}
_INT_FIELDS = {"healpix_level", "init_iter", "max_iter"}
_CR_KEYS = {"zeta", "cr_field", "cr_norm", "cr_n0", "cr_slope"}


def _norm(comment: str) -> str:
    """Normalise an inline comment for matching (case/space-insensitive)."""
    return re.sub(r"\s+", " ", comment).strip().lower()


_COMMENT_TO_KEY: dict[str, str] = {_norm(c): k for k, c in _COMMENTS.items()}

# Known alternative label spellings used by other 3D-PDR / 3D-XDR versions, so
# their (reworded) labels still map to the right field. Add new variants here;
# count-matching sections are also handled positionally (see _resolve_block), so
# aliases are only essential for a section that has an *inserted* new entry.
_COMMENT_ALIASES: dict[str, str] = {
    "G0 (in Draine field units) -x to +x in 1D": "fuv",
    "Cosmic-ray ionization rate (s^-1)": "zeta",
    "Dust-to-gas normalized to 1e-2": "d2g",
    "Turbulent velocity (km/s)": "vturb",
    "End time (yr)": "tend",
    "Av conversion factor": "av_fac",
    "Redshift (for CMB temperature)": "redshift",
    "Av critical (SUPRATHERMAL switch)": "av_crit",
    "Alfven velocity (km/s) (SUPRATHERMAl switch)": "v_alfv",
    "HEALPix level of refinement": "healpix_level",
    "Theta critical (0<phi<pi/2)": "theta_crit",
    "relative abundance tolerance": "rel_tol",
    "absolute abundance tolerance": "abs_tol",
    "First set of Chemical Iterations": "init_iter",
    "Total iterations": "max_iter",
    "Gas temperature for isothermal models": "tgas",
    "Floor temperature (if <Tcmb it's set to Tcmb)": "tmin",
    "Maximum allowed gas temperature": "tmax",
    "Dust temperature for isothermal dust models": "tdust",
    "Fcrit (% Accuracy)": "fcrit",
    "Tdiff (maximum temperature difference)": "tdiff",
}
for _c, _k in _COMMENT_ALIASES.items():
    _COMMENT_TO_KEY.setdefault(_norm(_c), _k)


def _cr_keys_for(cratten: int) -> list[str]:
    """Cosmic-ray field keys for a CRATTENUATION layout (0 -> ζ, 1 -> letter, 2 -> 3 nums)."""
    if cratten == 2:
        return ["cr_norm", "cr_n0", "cr_slope"]
    if cratten == 1:
        return ["cr_field"]
    return ["zeta"]


def _pdr_keys(cratten: int) -> list[str]:
    return (["fuv"] + _cr_keys_for(cratten)
            + ["d2g", "vturb", "tend", "grain_radius", "av_fac", "uv_fac",
               "redshift", "av_crit", "v_alfv"])


def _schema_for_title(title: str, cratten: int) -> list[str]:
    """Expected ordered keys for a section, matched by *title keyword* (robust to
    reworded titles and labels). Unknown sections return [] (entries -> extras)."""
    t = title.lower()
    if "input" in t or "output" in t:
        return ["indir", "input", "outdir", "output"]
    if "pdr" in t:
        return _pdr_keys(cratten)
    if "ray" in t:
        return ["healpix_level", "theta_crit"]
    if "ode" in t or "solver" in t:
        return ["rel_tol", "abs_tol"]
    if "iteration" in t:
        return ["init_iter", "max_iter"]
    if "thermal" in t:
        return ["tgas", "tmin", "tmax", "tdust", "fcrit", "tdiff"]
    return []


def _resolve_block(entries: list[dict], expected: list[str]) -> list[str | None]:
    """Resolve each entry of a section to a known key (or None = extra):

    1. every label maps to a unique known key -> use those (standard/aliased);
    2. else the entry count equals the expected count -> map positionally
       (handles fully reworded labels with no inserted/removed lines);
    3. else map the labels we recognise, others become extras (handles a
       genuinely new entry inserted among known ones, e.g. an X-ray flux line).
    """
    by_comment = [_COMMENT_TO_KEY.get(_norm(e["comment"])) for e in entries]
    if (expected and all(k is not None for k in by_comment)
            and len(set(by_comment)) == len(by_comment)):
        return by_comment
    if expected and len(entries) == len(expected):
        return list(expected)
    return by_comment


@dataclass
class Extra:
    """A params.dat entry with no curated mapping — auto-imported and preserved."""
    key: str        # stable id used as the widget/value key
    comment: str    # the inline !label (shown as the field label)
    value: str      # current raw value token
    section: str    # section title it lives under


def _extra_key(section: str, comment: str) -> str:
    return f"extra::{_norm(section)}::{_norm(comment)}"


def _coerce(key: str, raw: str):
    """Convert a raw token to the type used for a known key (fallback: default)."""
    try:
        if key in _STR_FIELDS:
            return raw.strip().upper() if key == "cr_field" else raw
        if key in _INT_FIELDS:
            return int(float(raw))
        return float(raw)
    except (ValueError, TypeError):
        return DEFAULTS.get(key, raw)


def _is_sep(s: str) -> bool:
    """A header separator line. Tolerant of decorated styles such as
    ``=====|-----`` or ``=====|`` (3D-XDR) as well as plain ``=====``: only the
    chars ``= - | space`` are allowed, with a run of at least three ``=`` or ``-``."""
    t = s.strip()
    if not t or set(t) - set("=-| "):
        return False
    return "===" in t or "---" in t


def _parse_file(text: str) -> list[dict]:
    """Parse params.dat into ordered section blocks::

        {"title": str, "raw_header": [3 lines], "is_coolant": bool,
         "entries": [{"value": str, "comment": str}], "coolant_lines": [raw str]}

    The cosmic-ray field and any new entries are ordinary entries; the coolant
    section (matched by title) keeps its lines verbatim so their comments survive.
    Section titles are taken from the text before any ``|`` column marker.
    """
    lines = text.splitlines()
    blocks: list[dict] = []
    cur: dict | None = None
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if _is_sep(line) and i + 2 < n and _is_sep(lines[i + 2]):
            title = lines[i + 1].split("|")[0].strip()
            cur = {"title": title, "raw_header": lines[i:i + 3],
                   "is_coolant": "coolant" in title.lower(),
                   "entries": [], "coolant_lines": []}
            blocks.append(cur)
            i += 3
            continue
        if line.strip() == "" or cur is None:
            i += 1
            continue
        if cur["is_coolant"]:
            cur["coolant_lines"].append(line)
        else:
            toks = line.split()
            comment = line.split("!", 1)[1].strip() if "!" in line else ""
            cur["entries"].append({"value": toks[0] if toks else "", "comment": comment})
        i += 1
    return blocks


def read_params(cratten: int = 0) -> dict:
    """Parse params.dat into a values dict keyed by curated id (matched on each
    value's inline label, so inserted lines don't shift our reads).

    Unrecognised entries are returned under ``"extras"`` as :class:`Extra`
    records; ``"coolants"`` holds the optional coolant filenames. ``cratten`` is
    accepted for backward compatibility but no longer needed (the cosmic-ray
    layout is detected from the labels).
    """
    values = dict(DEFAULTS)
    values["coolants"] = []
    values["extras"] = []
    path = paths.params_dat()
    if not path.exists():
        return values

    for blk in _parse_file(path.read_text()):
        if blk["is_coolant"]:
            cools = [ln.split()[0] for ln in blk["coolant_lines"] if ln.split()]
            values["coolants"] = [c for c in cools if c not in _FIXED_FILES]
            continue
        expected = _schema_for_title(blk["title"], cratten)
        resolved = _resolve_block(blk["entries"], expected)
        for e, key in zip(blk["entries"], resolved):
            if key:
                values[key] = _coerce(key, e["value"])
            else:
                values["extras"].append(Extra(
                    key=_extra_key(blk["title"], e["comment"]),
                    comment=e["comment"] or "(unlabelled entry)",
                    value=e["value"], section=blk["title"]))
    return values


# ---------------------------------------------------------------------------
# params.dat — write
# ---------------------------------------------------------------------------
def _emit(text: str, comment: str) -> str:
    return f"{text:<{_COMMENT_COL}}!{comment}" if comment else text


def _cr_lines(v: dict, cratten: int) -> list[str]:
    """Cosmic-ray field line(s) for the given CRATTENUATION layout."""
    if cratten == 2:
        return [_emit(_fmt(k, v[k]), _COMMENTS[k]) for k in ("cr_norm", "cr_n0", "cr_slope")]
    if cratten == 1:
        return [_emit(str(v["cr_field"]).strip().upper(), _COMMENTS["cr_field"])]
    return [_emit(_fmt("zeta", v["zeta"]), _COMMENTS["zeta"])]


def write_params(values: dict, cratten: int = 0) -> None:
    """Write params.dat from ``values``.

    If a params.dat already exists it is rewritten **structure-preserving**: the
    section order, headers, inline comments and any unknown ("extra") entries are
    kept in place and only the values are substituted — so a newer code version
    that added parameter lines round-trips safely. The cosmic-ray block is
    re-emitted for the current ``cratten`` and the final coolant block from the
    coolant selection. A fresh file (none present) uses the built-in template.
    """
    v = {**DEFAULTS, **values}
    path = paths.params_dat()
    if not path.exists():
        _write_template(v, cratten)
        return

    out: list[str] = []
    for blk in _parse_file(path.read_text()):
        out.extend(blk["raw_header"])
        if blk["is_coolant"]:
            orig = {ln.split()[0]: ln for ln in blk["coolant_lines"] if ln.split()}
            for _, fn in FIXED_COOLANTS:
                out.append(orig.get(fn, fn))  # preserve the file's line (with any comment)
            for fn in v.get("coolants", []):
                if fn not in _FIXED_FILES:
                    out.append(orig.get(fn, fn))
            continue
        expected = _schema_for_title(blk["title"], cratten)
        resolved = _resolve_block(blk["entries"], expected)
        # If the file's cosmic-ray layout already matches CRATTENUATION, keep its
        # lines (and labels) as ordinary fields; otherwise reconstruct the block.
        cr_native = [k for k in resolved if k in _CR_KEYS] == _cr_keys_for(cratten)
        cr_done = False
        for e, key in zip(blk["entries"], resolved):
            if key in _CR_KEYS and not cr_native:
                if not cr_done:
                    out.extend(_cr_lines(v, cratten))
                    cr_done = True
                continue
            if key:
                out.append(_emit(_fmt(key, v[key]), e["comment"]))
            else:
                ekey = _extra_key(blk["title"], e["comment"])
                out.append(_emit(str(v.get(ekey, e["value"])), e["comment"]))
    path.write_text("\n".join(out) + "\n")


def _write_template(v: dict, cratten: int) -> None:
    """Create a fresh params.dat from the built-in layout (no existing file)."""
    out: list[str] = []

    def block(title: str) -> None:
        out.extend([_SEP, title, _SEP])

    def emit(key: str, raw: str | None = None) -> None:
        out.append(_emit(raw if raw is not None else _fmt(key, v[key]), _COMMENTS.get(key, "")))

    block("Input/Output - Densities")
    for key in ("indir", "input", "outdir", "output"):
        emit(key, str(v[key]))
    block("PDR parameters")
    emit("fuv")
    out.extend(_cr_lines(v, cratten))
    for key in ("d2g", "vturb", "tend", "grain_radius", "av_fac", "uv_fac",
                "redshift", "av_crit", "v_alfv"):
        emit(key)
    block("Ray-tracing parameters")
    emit("healpix_level"); emit("theta_crit")
    block("ODE Solver parameters")
    emit("rel_tol"); emit("abs_tol")
    block("Chemistry iterations")
    emit("init_iter"); emit("max_iter")
    block("Thermal balance values")
    for key in ("tgas", "tmin", "tmax", "tdust", "fcrit", "tdiff"):
        emit(key)
    block("Coolant files")
    for _, fn in FIXED_COOLANTS:
        out.append(fn)
    for fn in v.get("coolants", []):
        if fn not in _FIXED_FILES:
            out.append(fn)
    paths.params_dat().write_text("\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# Initial abundances in species_*.d
# ---------------------------------------------------------------------------
def read_abundances(network: str) -> dict[str, float]:
    """Return ``{token: current_abundance}`` for the target elements that
    actually exist in this network's species file, in display order."""
    path = species_file(network)
    present: dict[str, float] = {}
    if not path.exists():
        return present
    file_vals: dict[str, float] = {}
    for line in path.read_text().splitlines():
        parts = line.strip().split(",")
        if len(parts) >= 3:
            token = parts[1].strip()
            try:
                file_vals[token] = float(parts[2])
            except ValueError:
                continue
    for token in ELEMENT_ORDER:
        if token in file_vals:
            present[token] = file_vals[token]
    return present


def write_abundances(network: str, abundances: dict[str, float]) -> list[str]:
    """Replace the abundance field of the given species tokens in the network's
    species file. Returns the list of tokens actually updated.

    The species file format is ``index,token,abundance,mass`` per line.
    """
    path = species_file(network)
    lines = path.read_text().splitlines(keepends=True)
    updated: list[str] = []
    out: list[str] = []
    for line in lines:
        raw = line.rstrip("\n")
        parts = raw.split(",")
        if len(parts) >= 4:
            token = parts[1].strip()
            if token in abundances:
                eol = "\n" if line.endswith("\n") else ""
                parts[2] = f"{abundances[token]:.2E}"
                out.append(",".join(parts) + eol)
                updated.append(token)
                continue
        out.append(line)
    path.write_text("".join(out))
    return updated
