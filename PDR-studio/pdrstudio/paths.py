"""Locate the 3D-PDR tree and its key files/directories.

PDR-studio lives inside the 3D-PDR root (``.../3D-PDR/PDR-studio``), so the
root is simply the parent of this package's parent. Everything else is derived
from that single anchor, and can be overridden with the ``PDR_ROOT`` env var
for unusual layouts.
"""
from __future__ import annotations

import os
from pathlib import Path

# .../3D-PDR/PDR-studio/pdrstudio/paths.py  -> root is 3 parents up
_DEFAULT_ROOT = Path(__file__).resolve().parents[2]


def pdr_root() -> Path:
    """Absolute path to the 3D-PDR root directory."""
    env = os.environ.get("PDR_ROOT")
    return Path(env).resolve() if env else _DEFAULT_ROOT


def app_dir() -> Path:
    """The PDR-studio app directory (where this package lives)."""
    return Path(__file__).resolve().parents[1]


def myplots_dir() -> Path:
    """``PDR-studio/myplots`` for saved figures (created on first use)."""
    d = app_dir() / "myplots"
    d.mkdir(parents=True, exist_ok=True)
    return d


def src_dir() -> Path:
    return pdr_root() / "src"


def config_mk() -> Path:
    return src_dir() / "config.mk"


def makefile() -> Path:
    return src_dir() / "makefile"


def params_dat() -> Path:
    return pdr_root() / "params.dat"


def chemfiles_dir() -> Path:
    return pdr_root() / "chemfiles"


def ics_dir() -> Path:
    return pdr_root() / "ics"


def sims_dir() -> Path:
    """Default output directory. The actual one is set in params.dat (outdir),
    but ``sims`` is the conventional default."""
    return pdr_root() / "sims"


def executable() -> Path:
    """The compiled 3D-PDR binary (built into the root by the makefile)."""
    return pdr_root() / "3DPDR"


def rt_executable() -> Path:
    return pdr_root() / "RTtool"


def rt_tool_dir() -> Path:
    """Directory holding the RT-tool Fortran sources + makefile."""
    return pdr_root() / "RT-tool"


def paramsrt_dat() -> Path:
    return pdr_root() / "paramsRT.dat"


def studio_dir() -> Path:
    return pdr_root() / "PDR-studio"


def reports_dir() -> Path:
    """Where debug-report tarballs are written."""
    d = studio_dir() / "reports"
    d.mkdir(parents=True, exist_ok=True)
    return d
