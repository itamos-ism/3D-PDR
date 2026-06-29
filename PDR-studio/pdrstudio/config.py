"""Read/write ``src/config.mk`` and adapt automatically to different 3D-PDR versions.

PDR-studio does **not** hard-code the flag list. On every load it *discovers* the
flags from the actual code checkout, so a newer ``config.mk`` (with new compile
flags) needs no PDR-studio code change:

* ``src/config.mk``
    - the ``KEY = VALUE`` assignments: which flags exist, their current values,
      and the order to show them in;
    - the documentation comment header: a per-flag help block, the enumerated
      values it lists (e.g. ``HTT91``/``0``), and section headers
      (``#---- Compiler options ----``) that group the flags.
* ``src/makefile``
    - the ``ifeq ($(KEY),VALUE)`` branches reveal the *full* set of allowed
      values for each flag (often more than the comments list — e.g. NETWORK has
      six values in the makefile but three in the comments).

The discovered metadata is then refined by a small :data:`CURATED` overlay
(nicer labels, concise help, which flags are *advanced* or *hidden*). A flag with
no overlay entry is still shown — under its documented section, tucked into the
Advanced area — and surfaced to the user via :func:`newly_discovered` so they
know a new flag was imported. To "promote" a new flag (give it a nice label and
group), add an entry to :data:`CURATED`; nothing else changes.

Reading/writing only ever touches ``KEY = VALUE`` lines; comments, ordering and
unknown keys are preserved verbatim, so the round-trip is safe across versions.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from . import paths

# A line like  "NETWORK              = REDUCED"  (ignores commented "#..." lines)
_ASSIGN_RE = re.compile(r"^(?P<key>[A-Za-z_][A-Za-z0-9_]*)\s*=\s*(?P<val>.*?)\s*$")
# A comment section header, e.g.  "#---------- Compiler options ----------"
_SECTION_RE = re.compile(r"^#[-=]{2,}(.+?)[-=]{2,}\s*$")
# A comment flag doc block start, e.g.  "#NETWORK   : REDUCED (...)"
_DOC_RE = re.compile(r"^#\s*([A-Za-z][A-Za-z0-9_]+)\s*:\s*(.*?)\s*$")
# An enumerated value at the start of a doc line: an upper-case/numeric token
# followed by ' - ', ' (', ':' or end. (Prose words are lower-case, so skipped.)
_VALUE_RE = re.compile(r"^([A-Z0-9_+]+)\s*(?:-|\(|:|$)")
# A makefile branch  "ifeq ($(KEY),VALUE)"  (ifneq is an exclusion, not an enum)
_IFEQ_RE = re.compile(r"\bifeq\s*\(\s*\$\(\s*(\w+)\s*\)\s*,\s*([^)]*?)\s*\)")


@dataclass
class Flag:
    """UI metadata for one config.mk key (discovered, then curated)."""

    key: str
    label: str
    kind: str  # "bool" | "choice" | "int" | "str"
    help: str = ""
    choices: list[str] = field(default_factory=list)
    advanced: bool = False     # True -> tucked under the "Advanced" expander
    section: str = ""          # group/section it belongs to
    source: str = "curated"    # "curated" or "discovered" (auto-imported)


# ===========================================================================
# Curated overlay — polish for the flags we know about. Everything here is
# OPTIONAL: discovery already provides a working schema. We mainly supply nice
# labels, concise help, group placement, and which flags are advanced/hidden.
# `choices`/`kind` are intentionally left to discovery so new valid values
# (e.g. extra networks) appear automatically; override them only when needed.
# ===========================================================================
# Order in which curated groups are shown. Discovered-only sections follow.
GROUP_ORDER = [
    "Compiler", "Network", "Physics", "Geometry",
    "Convergence (recommended on)", "Run mode & output",
]

# Flags never shown in the UI (but preserved in the file). GUESS_TEMP is derived
# from THERMALBALANCE (see normalize_dependencies); CPPFLAGS is a raw build flag;
# XRAYS is intentionally not exposed.
HIDDEN = {"GUESS_TEMP", "CPPFLAGS", "XRAYS", "PYWRAP"}

CURATED: dict[str, dict] = {
    # --- Compiler ---
    "F90": dict(label="Fortran compiler", kind="str", group="Compiler",
                help="Fortran compiler command (e.g. gfortran, ifort)."),
    "CC": dict(label="C compiler", kind="str", group="Compiler",
               help="C compiler command (e.g. gcc)."),
    "OPENMP": dict(label="OpenMP parallel", group="Compiler",
                   help="1 = run in parallel (OpenMP); 0 = serial (1 CPU)."),
    "OPTIMISE": dict(label="Optimisation level", group="Compiler",
                     help="0/1 for development & debugging; 2/3 for production runs."),
    # --- Network (visible) / Geometry (advanced) ---
    "NETWORK": dict(label="Chemical network", group="Network",
                    help="Chemical network to compile (species/reactions vary per network)."),
    "DIMENSIONS": dict(label="Dimensions", group="Geometry", advanced=True,
                       help="1 = 1D models (requires RAYTHEIA=0); 3 = 3D models."),
    "RAYTHEIA": dict(label="Ray-tracing engine", group="Geometry", advanced=True,
                     help="1 = on, memory-optimised; 2 = on, no mem-opt (faster, more RAM); 0 = off."),
    "XYZ": dict(label="Coordinate ordering", group="Geometry", advanced=True,
                kind="choice", help="0 = z changes first (xyz); 1 = x changes first."),
    # --- Physics ---
    "DUST": dict(label="Dust temperature", group="Physics",
                 help="HTT91 = Hollenbach, Takahashi & Tielens (1991); 0 = isothermal dust."),
    "THERMALBALANCE": dict(label="Thermal balance", group="Physics",
                           help="1 = iterate thermal balance; 0 = isothermal run. "
                                "GUESS_TEMP follows this flag automatically."),
    "H2FORM": dict(label="H2 formation", group="Physics",
                   help="CT02 = Cazaux & Tielens; SIMPLE = 3e-18·√T·exp(−T/1e3); R07 = Roellig+07."),
    "GRAINRECOMB": dict(label="Grain recombination", group="Physics", advanced=True,
                        help="1 = electron recombination on dust grains."),
    "SUPRATHERMAL": dict(label="Suprathermal CO", group="Physics",
                         help="1 = suprathermal formation of CO via CH+."),
    "CRATTENUATION": dict(label="CR attenuation", group="Physics",
                          help="1 = model cosmic-ray attenuation (L/H of Padovani+); 0 = off (default)."),
    # --- Convergence (all advanced) ---
    "FORCECONVERGENCE": dict(label="Force convergence", group="Convergence (recommended on)", advanced=True,
                             help="1 = helps fast convergence (recommended)."),
    "NGACCEL": dict(label="Ng acceleration", group="Convergence (recommended on)", advanced=True,
                    help="1 = Ng (1974) acceleration of level-population iterations (recommended)."),
    "ILLINOIS": dict(label="Illinois regula-falsi", group="Convergence (recommended on)", advanced=True,
                     help="1 = Illinois-modified regula-falsi in ln(T) for thermal balance (recommended)."),
    "COOLANTFREEZE": dict(label="Freeze converged coolants", group="Convergence (recommended on)", advanced=True,
                          help="1 = skip re-solving coolants already converged at a point (recommended)."),
    "FREEZECOOL": dict(label="Freeze on cooling rate", group="Convergence (recommended on)", advanced=True,
                       help="1 = freeze a coolant once its cooling rate stabilises (needs COOLANTFREEZE; recommended)."),
    "EMISGUARD": dict(label="Emissivity NaN guard", group="Convergence (recommended on)", advanced=True,
                      help="1 = guard the (S-BB)/S net line-cooling factor against 0/0 NaN (recommended)."),
    "TMINRECOVERY": dict(label="Tmin recovery", group="Convergence (recommended on)", advanced=True,
                         help="1 = keep the upper T bracket at Tmin so a point can heat back up (recommended)."),
    "DETBALANCE": dict(label="Detailed balance", group="Convergence (recommended on)", advanced=True,
                       help="1 = enforce detailed balance on collisional rates at local Tgas (recommended)."),
    "SOBOLEV": dict(label="Sobolev net rates", group="Convergence (recommended on)", advanced=True,
                    help="1 = Sobolev/ALI-like net radiative rates instead of the lagged mean field (recommended)."),
    "REBRACKET": dict(label="Re-bracket on stall", group="Convergence (recommended on)", advanced=True,
                      help="1 = re-open the thermal-balance bracket when the search stalls (recommended)."),
    "NGCYCLE": dict(label="Damp flip-flop pops", group="Convergence (recommended on)", advanced=True,
                    help="1 = average the last two iterates when a level-pop update reverses (needs NGACCEL)."),
    "NGRELAX": dict(label="Adaptive under-relaxation", group="Convergence (recommended on)", advanced=True,
                    help="1 = adaptive per-point under-relaxation of oscillating pops (needs NGACCEL; recommended for 3D)."),
    "MASERFIX": dict(label="Maser fix", group="Convergence (recommended on)", advanced=True,
                     help="1 = treat inverted (masing) lines as freely escaping (beta=1) (recommended)."),
    "CHEMSTEADY": dict(label="Chemistry to steady state", group="Convergence (recommended on)", advanced=True,
                       help="1 = integrate chemistry in doubling-time blocks until <1% change (recommended)."),
    "CHEMRETRY": dict(label="CVODE retry", group="Convergence (recommended on)", advanced=True,
                      help="1 = re-integrate cells whose CVODE solve fails, with a forced shrinking step (recommended)."),
    # --- Run mode & output (all advanced) ---
    "RESTART": dict(label="Restart run", group="Run mode & output", advanced=True,
                    help="1 = checkpoint/restart an interrupted model; 0 = off (default)."),
    "OUTRAYINFO": dict(label="Per-ray output", group="Run mode & output", advanced=True,
                       help="1 = save per-ray info to separate files (3D); 0 = off."),
    "CHEMANALYSIS": dict(label="Chemical-analysis output", group="Run mode & output", advanced=True,
                         help="1 = output formation/destruction analysis per point; 0 = off."),
}


# ===========================================================================
# Low-level parsers (pure functions on file text, so they are easy to test)
# ===========================================================================
def _parse_assignments(text: str) -> list[tuple[str, str]]:
    """Ordered ``[(key, value), ...]`` for every ``KEY = VALUE`` line."""
    out: list[tuple[str, str]] = []
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            continue
        m = _ASSIGN_RE.match(line)
        if m:
            out.append((m.group("key"), m.group("val")))
    return out


def _lead_value(s: str) -> str | None:
    """Enumerated value at the start of a doc fragment, or None (prose)."""
    s = s.strip().lstrip(":").strip()
    m = _VALUE_RE.match(s)
    return m.group(1) if m else None


def _parse_comment_docs(text: str) -> dict[str, dict]:
    """Per-flag ``{key: {section, help, values}}`` parsed from the comment header.

    A doc block starts at ``#KEY : ...`` and runs over the indented continuation
    comment lines until the next flag or section header. Enumerated values are
    the upper-case/numeric tokens that start each line of the block.
    """
    docs: dict[str, dict] = {}
    section = ""
    cur: str | None = None
    for line in text.splitlines():
        if not line.startswith("#"):
            cur = None
            continue
        sec = _SECTION_RE.match(line)
        if sec:
            name = sec.group(1).strip().strip("-=").strip()
            name = re.sub(r"(?i)\s*options?$", "", name).strip()
            if name:
                section = name
            cur = None
            continue
        doc = _DOC_RE.match(line)
        if doc:
            cur = doc.group(1)
            d = docs.setdefault(cur, {"section": section, "help": "", "values": []})
            d["section"] = section
            rest = doc.group(2)
            if rest:
                d["help"] = rest
                v = _lead_value(rest)
                if v and v not in d["values"]:
                    d["values"].append(v)
            continue
        if cur is not None:  # continuation line of the current flag's block
            frag = line[1:].strip()
            d = docs[cur]
            v = _lead_value(frag)
            if v and v not in d["values"]:
                d["values"].append(v)
            if frag and not frag.startswith("-"):
                d["help"] = (d["help"] + " " + frag).strip()
    # tidy help whitespace
    for d in docs.values():
        d["help"] = re.sub(r"\s+", " ", d["help"]).strip()
    return docs


def _parse_makefile_values(text: str) -> dict[str, list[str]]:
    """``{key: [values]}`` from ``ifeq ($(KEY),VALUE)`` branches (order preserved)."""
    out: dict[str, list[str]] = {}
    for key, val in _IFEQ_RE.findall(text):
        val = val.strip()
        if not val:
            continue
        out.setdefault(key, [])
        if val not in out[key]:
            out[key].append(val)
    return out


# ===========================================================================
# Inference + merge
# ===========================================================================
def _infer(key: str, value: str, doc: dict, mk_values: list[str]) -> dict:
    """Infer kind/choices/help/section for a flag from discovery sources."""
    # union of enumerated values: comment-documented first, then makefile extras
    values: list[str] = list(doc.get("values", []))
    for v in mk_values:
        if v not in values:
            values.append(v)

    valset = set(values)
    if values and valset <= {"0", "1"}:
        kind = "bool"
    elif len(values) >= 2:
        kind = "choice"
    elif value in ("0", "1"):
        kind = "bool"
    else:
        kind = "str"

    # Always carry the value set (a curated `kind` override may turn a 0/1 flag
    # into a choice), and never drop the current value from a choice list.
    if value and value not in values:
        values = [value] + values

    return dict(kind=kind, choices=values,
                help=doc.get("help", ""), section=doc.get("section", "") or "Other")


def discover_flags(config_text: str | None = None,
                   makefile_text: str | None = None) -> list[Flag]:
    """Discover the flags from config.mk + makefile and apply the curated overlay.

    Returns an ordered list (config.mk assignment order), excluding HIDDEN flags.
    """
    if config_text is None:
        p = paths.config_mk()
        config_text = p.read_text() if p.exists() else ""
    if makefile_text is None:
        p = paths.makefile()
        makefile_text = p.read_text() if p.exists() else ""

    docs = _parse_comment_docs(config_text)
    mk_values = _parse_makefile_values(makefile_text)
    flags: list[Flag] = []
    for key, value in _parse_assignments(config_text):
        if key in HIDDEN:
            continue
        inf = _infer(key, value, docs.get(key, {}), mk_values.get(key, []))
        cur = CURATED.get(key, {})
        flags.append(Flag(
            key=key,
            label=cur.get("label", key),
            kind=cur.get("kind", inf["kind"]),
            help=cur.get("help") or inf["help"],
            choices=cur.get("choices", inf["choices"]),
            advanced=cur.get("advanced", key not in CURATED),  # new flags -> advanced
            section=cur.get("group", inf["section"]),
            source="curated" if key in CURATED else "discovered",
        ))
    return flags


def build_schema() -> dict[str, list[Flag]]:
    """Discovered+curated flags grouped by section, in display order."""
    groups: dict[str, list[Flag]] = {}
    for flag in discover_flags():
        groups.setdefault(flag.section, []).append(flag)
    # curated groups first (in GROUP_ORDER), then any discovered-only sections
    ordered = [g for g in GROUP_ORDER if g in groups]
    ordered += [g for g in groups if g not in GROUP_ORDER]
    return {g: groups[g] for g in ordered}


def visible_groups() -> dict[str, list[Flag]]:
    """Groups restricted to their always-visible (non-advanced) flags."""
    out: dict[str, list[Flag]] = {}
    for group, flags in build_schema().items():
        shown = [f for f in flags if not f.advanced]
        if shown:
            out[group] = shown
    return out


def advanced_groups() -> dict[str, list[Flag]]:
    """Groups restricted to their advanced flags."""
    out: dict[str, list[Flag]] = {}
    for group, flags in build_schema().items():
        adv = [f for f in flags if f.advanced]
        if adv:
            out[group] = adv
    return out


def newly_discovered() -> list[Flag]:
    """Flags present in config.mk but not in the curated overlay (auto-imported)."""
    return [f for f in discover_flags() if f.source == "discovered"]


# ===========================================================================
# Dependencies, read, write (KEY = VALUE only; comments/order/unknowns kept)
# ===========================================================================
def normalize_dependencies(values: dict[str, str]) -> list[str]:
    """Enforce inter-flag dependencies in place. Returns notes on any changes.

    GUESS_TEMP follows THERMALBALANCE (the makefile errors on TB=0 & GUESS_TEMP=1),
    so it is not shown in the UI and is derived here.
    """
    notes: list[str] = []
    if "THERMALBALANCE" in values:
        tb = values.get("THERMALBALANCE", "1").strip()
        want = "1" if tb == "1" else "0"
        if values.get("GUESS_TEMP", "").strip() != want:
            values["GUESS_TEMP"] = want
            notes.append(f"GUESS_TEMP set to {want} to follow THERMALBALANCE={tb}.")
    return notes


def read_config(path: Path | None = None) -> dict[str, str]:
    """Return the current ``KEY = VALUE`` assignments as a dict of strings."""
    path = path or paths.config_mk()
    return dict(_parse_assignments(path.read_text()))


def write_config(values: dict[str, str], path: Path | None = None) -> None:
    """Rewrite ``config.mk`` in place, updating only the values of existing
    assignment lines. Comments, ordering, and unknown keys are preserved.
    """
    path = path or paths.config_mk()
    lines = path.read_text().splitlines(keepends=True)
    out: list[str] = []
    for line in lines:
        stripped = line.rstrip("\n")
        if not stripped.lstrip().startswith("#"):
            m = _ASSIGN_RE.match(stripped)
            if m and m.group("key") in values:
                key = m.group("key")
                eol = "\n" if line.endswith("\n") else ""
                out.append(f"{key:<20} = {values[key]}{eol}")
                continue
        out.append(line)
    path.write_text("".join(out))
