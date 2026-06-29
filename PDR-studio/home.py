"""Home / welcome page (rendered by the navigation entrypoint, app.py)."""
from __future__ import annotations

import streamlit as st

from pdrstudio import config, paths, ui

ui.configure_page("Home")

if ui.LOGO.exists():
    st.image(str(ui.LOGO), width=360)
else:
    st.title("🌌 PDR-studio")
st.markdown(
    "A friendly, browser-based control panel for **3D-PDR** — configure, compile, "
    "run, and analyse photodissociation-region models without touching the terminal."
)

# --- Environment / detection ------------------------------------------------
root = paths.pdr_root()
exe = paths.executable()
cfg_path = paths.config_mk()
params_path = paths.params_dat()

cfg = config.read_config() if cfg_path.exists() else {}

col1, col2, col3 = st.columns(3)
with col1:
    st.metric("3D-PDR root", "detected" if root.exists() else "missing")
    st.caption(f"`{root}`")
with col2:
    st.metric("Executable", "built" if exe.exists() else "not built")
    st.caption(f"last built: {ui.mtime_str(exe)}")
with col3:
    net = cfg.get("NETWORK", "—")
    dims = cfg.get("DIMENSIONS", "—")
    st.metric("Current build", f"{dims}D · {net}")
    st.caption(f"config.mk: {ui.mtime_str(cfg_path)}")

st.divider()

# --- Quick file checks ------------------------------------------------------
st.subheader("Environment")
checks = {
    "config.mk": cfg_path,
    "params.dat": params_path,
    "chemfiles/": paths.chemfiles_dir(),
    "ics/": paths.ics_dir(),
    "src/": paths.src_dir(),
}
cols = st.columns(len(checks))
for (label, p), c in zip(checks.items(), cols):
    c.write(f"{ui.status_badge(p.exists())}  \n`{label}`")

st.divider()

# --- Feature guide ----------------------------------------------------------
st.subheader("Features")
st.markdown(
    """
Use the **sidebar** on the left to move between features:

| | Feature | What it does |
|---|---|---|
| ⚙️ | **Code Configuration** | Edit `config.mk` flags and compile the code. |
| 📐 | **Model Parameters** | Edit `params.dat`, set elemental abundances & reaction rates, run the model. |
| 📡 | **Line emission** | Brightness temperatures and optical depths for all modelled lines/transitions. |
| 📊 | **Abundance Profiles** | Gas-temperature & abundance profiles; plus a **Custom plot** tab to code your own matplotlib figures. |
| 🔥 | **Heating & Cooling** | Heating & cooling functions, coolant level populations, emissivities (+ custom plot). |
| 🧪 | **Chemical Analysis** | Formation & destruction pathways of any species, per depth point (optional local-LLM summary). |
| 🧬 | **Network Builder** | Build a new chemical network with MakeRates from selected elements, then import it into 3D-PDR. |
| 🐞 | **Report a bug** | Bundle the version, parameters and outputs into a tarball to send to the 3D-PDR team. |
"""
)

st.success("**Status:** all features are live.", icon="✅")
