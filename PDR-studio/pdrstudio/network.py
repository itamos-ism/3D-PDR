"""Resolve the active chemistry-network files in ``chemfiles/``.

The species and reaction files share a single suffix:
``species_<suffix>.d`` and ``rates_<suffix>.d`` always go together. The suffix
depends on the compiled ``NETWORK`` and on whether X-rays are enabled — an X-ray
build reads the ``X`` variant (``species_Xreduced.d`` / ``rates_Xreduced.d``),
as the makefile (``odes_Xreduced.o``) and ``read_species.F90`` show.

Everything is resolved against the files that actually exist and the current
``config.mk``, so PDR-studio adapts to whatever networks a given 3D-PDR / 3D-XDR
version ships (e.g. a new ``Xreduced`` network appearing when ``XRAYS=1``).
"""
from __future__ import annotations

from pathlib import Path

from . import config, paths

# NETWORK flag value -> base chemfiles suffix.
_BASE_SUFFIX: dict[str, str] = {
    "REDUCED": "reduced.d",
    "MEDIUM": "medium.d",
    "MEDIUMS": "mediumS.d",
    "SCUTUM": "scutum.d",
    "FULL": "full.d",
    "MYNETWORK": "mynetwork.d",
}


def _cfg(cfg: dict | None) -> dict:
    if cfg is not None:
        return cfg
    return config.read_config() if paths.config_mk().exists() else {}


def base_suffix(network: str) -> str:
    """Base suffix for a NETWORK value (e.g. 'REDUCED' -> 'reduced.d').

    Built-ins use their conventional lower-case suffix; a custom network
    (e.g. from the Network Builder) keeps the name's case — its files are
    ``species_<NAME>.d`` — so we prefer whichever file actually exists."""
    n = (network or "").strip()
    if n.upper() in _BASE_SUFFIX:
        return _BASE_SUFFIX[n.upper()]
    for cand in (f"{n}.d", f"{n.lower()}.d"):
        if species_path(cand).exists() or rates_path(cand).exists():
            return cand
    return f"{n}.d"


def species_path(suffix: str) -> Path:
    return paths.chemfiles_dir() / f"species_{suffix}"


def rates_path(suffix: str) -> Path:
    return paths.chemfiles_dir() / f"rates_{suffix}"


def xrays_enabled(cfg: dict | None = None) -> bool:
    return str(_cfg(cfg).get("XRAYS", "0")).strip() == "1"


def active_network(cfg: dict | None = None) -> str:
    return str(_cfg(cfg).get("NETWORK", "REDUCED")).strip()


def suffix_for(network: str, xrays: bool | None = None) -> str:
    """chemfiles suffix for a network. With X-rays on, prefer the ``X`` variant
    (``species_Xreduced.d``) when that file is present, else fall back to base."""
    if xrays is None:
        xrays = xrays_enabled()
    base = base_suffix(network)
    if xrays:
        xbase = "X" + base
        if species_path(xbase).exists() or rates_path(xbase).exists():
            return xbase
    return base


def active_suffix(cfg: dict | None = None) -> str:
    """Suffix for the currently-configured build (NETWORK + XRAYS)."""
    cfg = _cfg(cfg)
    return suffix_for(active_network(cfg), xrays_enabled(cfg))


def species_list(suffix: str) -> list[str]:
    """Species tokens (file order) from ``species_<suffix>.d`` whose lines are
    ``index,token,abundance,mass``. Empty if the file is missing."""
    p = species_path(suffix)
    out: list[str] = []
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        parts = line.split(",")
        if len(parts) >= 2 and parts[0].strip().isdigit():
            out.append(parts[1].strip())
    return out


def available_suffixes() -> list[str]:
    """All ``species_*.d`` suffixes present in ``chemfiles/`` (the real networks,
    not the ``*.d_xxx`` backups)."""
    d = paths.chemfiles_dir()
    if not d.is_dir():
        return []
    return sorted(p.name[len("species_"):] for p in d.glob("species_*.d"))


def label_for_suffix(suffix: str | None) -> str | None:
    """Human label for a suffix, e.g. 'reduced.d' -> 'REDUCED',
    'Xreduced.d' -> 'REDUCED+Xrays'."""
    if not suffix:
        return None
    for net, base in _BASE_SUFFIX.items():
        if suffix == base:
            return net
        if suffix == "X" + base:
            return f"{net}+Xrays"
    return suffix[:-2] if suffix.endswith(".d") else suffix


def detect_suffix(nspec: int, cfg: dict | None = None) -> str | None:
    """Suffix whose species count equals ``nspec``. Prefer the currently
    configured one (X-ray aware); otherwise any available species file matching."""
    act = active_suffix(cfg)
    if len(species_list(act)) == nspec:
        return act
    for suf in available_suffixes():
        if len(species_list(suf)) == nspec:
            return suf
    return None
