"""Build a debug-report tarball capturing everything needed to reproduce an
issue: the source tree's config, inputs, chemistry files, and (optionally) the
outputs of the offending model.

The goal is a single file the user can email to the 3D-PDR team. We snapshot
the *current* state of the small, text inputs always, and include heavier
outputs only when asked.
"""
from __future__ import annotations

import datetime as _dt
import platform
import tarfile
from io import BytesIO
from pathlib import Path

from . import paths

# Small, always-relevant inputs (relative to root) captured in every report.
_CORE_FILES = [
    "src/config.mk",
    "src/makefile",
    "params.dat",
]
# Directories whose text inputs we snapshot wholesale.
_CORE_DIRS = ["chemfiles", "ics"]


def _manifest(note: str, model_prefix: str | None) -> str:
    now = _dt.datetime.now().isoformat(timespec="seconds")
    lines = [
        "3D-PDR / PDR-studio debug report",
        f"created: {now}",
        f"platform: {platform.platform()}",
        f"python: {platform.python_version()}",
        f"root: {paths.pdr_root()}",
        f"model_prefix: {model_prefix or '(none)'}",
        "",
        "user note:",
        note.strip() or "(none provided)",
    ]
    return "\n".join(lines) + "\n"


def build_report(
    note: str = "",
    model_prefix: str | None = None,
    include_outputs: bool = True,
) -> Path:
    """Create a ``.tar.gz`` under ``PDR-studio/reports`` and return its path.

    Always includes core inputs (config.mk, makefile, params.dat, chemfiles,
    ics). If ``include_outputs`` and ``model_prefix`` are given, also bundles
    the matching ``sims/<prefix>.*`` output files.
    """
    root = paths.pdr_root()
    stamp = _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    tag = (model_prefix or "report").replace("/", "_")
    out_path = paths.reports_dir() / f"pdrstudio-report-{tag}-{stamp}.tar.gz"

    with tarfile.open(out_path, "w:gz") as tar:
        # Manifest written from memory.
        manifest = _manifest(note, model_prefix).encode()
        info = tarfile.TarInfo(name="REPORT_MANIFEST.txt")
        info.size = len(manifest)
        tar.addfile(info, BytesIO(manifest))

        for rel in _CORE_FILES:
            p = root / rel
            if p.exists():
                tar.add(p, arcname=rel)

        for rel in _CORE_DIRS:
            d = root / rel
            if d.is_dir():
                tar.add(d, arcname=rel)

        if include_outputs and model_prefix:
            sims = paths.sims_dir()
            if sims.is_dir():
                for f in sorted(sims.glob(f"{model_prefix}.*")):
                    tar.add(f, arcname=f"sims/{f.name}")

    return out_path
