"""Shared Streamlit UI helpers used across pages."""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

import streamlit as st

from . import paths

PAGE_ICON = "🌌"  # emoji fallback if the logo image is missing
# The PDR-studio logo lives at the app root (PDR-studio/PDR-studio_logo.png);
# resolve it from this module so it's found regardless of the working directory.
LOGO = Path(__file__).resolve().parent.parent / "PDR-studio_logo.png"


def page_icon():
    return str(LOGO) if LOGO.exists() else PAGE_ICON


def configure_page(title: str) -> None:
    """Set the page config. Under ``st.navigation`` the entrypoint (app.py) has
    already configured the page, so a second call here is ignored; when a page is
    run on its own (e.g. tests) this still sets a sensible wide layout."""
    try:
        st.set_page_config(page_title=f"PDR-studio · {title}", page_icon=page_icon(),
                           layout="wide")
    except Exception:  # noqa: BLE001 — already configured by the entrypoint
        pass


def sidebar_header() -> None:
    """The PDR-studio logo + caption shown in the sidebar (called once from the
    navigation entrypoint, so it appears on every page)."""
    if LOGO.exists():
        st.logo(str(LOGO), size="large")
    with st.sidebar:
        if not LOGO.exists():
            st.markdown(f"## {PAGE_ICON} PDR-studio")
        st.caption(f"3D-PDR control panel\n\n`{paths.pdr_root()}`")
        st.divider()


def mtime_str(p: Path) -> str:
    """Human-readable modification time, or '—' if missing."""
    if not p.exists():
        return "—"
    return _dt.datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")


def status_badge(ok: bool, ok_text: str = "ready", bad_text: str = "missing") -> str:
    return f"✅ {ok_text}" if ok else f"❌ {bad_text}"


def fuv_cr_header(p: dict, cfg: dict) -> None:
    """Compact FUV + cosmic-ray summary row from params.dat / config.mk."""
    fuv = float(p.get("fuv", 1.0))
    habing = fuv * 1.71
    cratten = int(cfg.get("CRATTENUATION", "0"))
    cols = st.columns(5)
    cols[0].metric("G₀ (Draine)", f"{fuv:g}", help="External FUV field in Draine (1974) units")
    cols[1].metric("G₀ (Habing)", f"{habing:.4g}", help="1 Draine ≈ 1.71 Habing")
    if cratten == 1:
        cols[2].metric("CR field (L/H/U)", str(p.get("cr_field", "L")),
                       help="CRATTENUATION=1: L = low, H = high, U = ultra-high")
    elif cratten == 2:
        cols[2].metric("CR norm", f"{float(p.get('cr_norm', 1.0)):g}")
        cols[3].metric("CR n₀ (cm⁻³)", f"{float(p.get('cr_n0', 1e3)):.2e}")
        cols[4].metric("CR slope", f"{float(p.get('cr_slope', 0.0)):g}")
    else:
        cols[2].metric("ζ_CR (s⁻¹)", f"{float(p.get('zeta', 1e-17)):.2e}",
                       help="Cosmic-ray ionisation rate (CRATTENUATION=0)")
