"""Feature — Tutorial.

A hands-on, in-app walkthrough: from launching PDR-studio to a finished,
analysed 3D-PDR model. Task-oriented (the README covers the architecture); it
follows the sidebar top-to-bottom, which is the natural workflow. A small live
"where am I" status row mirrors the Home checks so the reader can see whether the
code is built / a model has been run while they read.
"""
from __future__ import annotations

import streamlit as st

from pdrstudio import config, paths, ui

ui.configure_page("Tutorial")
st.title("📖 PDR-studio — a hands-on tutorial")
st.markdown(
    "This walkthrough takes you from launching **PDR-studio** to a finished, "
    "analysed 3D-PDR model — all in the browser. Follow it top to bottom for "
    "your first model. For architecture and developer notes see `README.md`.")

# --- live status row (mirrors Home, so the reader sees their own state) ------
root = paths.pdr_root()
exe = paths.executable()
cfg = config.read_config() if paths.config_mk().exists() else {}
s1, s2, s3 = st.columns(3)
s1.metric("3D-PDR root", "detected" if root.exists() else "missing")
s2.metric("Executable", "built" if exe.exists() else "not built")
s3.metric("Current build", f"{cfg.get('DIMENSIONS', '—')}D · {cfg.get('NETWORK', '—')}")
st.caption("Tip: keep this page open in one tab and follow along in the others.")
st.divider()

tab_start, tab_model, tab_analyse, tab_further, tab_help = st.tabs(
    ["🚀 Start here", "🧭 Your first model", "📈 Analyse results",
     "🧬 Going further", "🛟 Troubleshooting"])

# ===========================================================================
# Start here
# ===========================================================================
with tab_start:
    st.header("What PDR-studio is")
    st.markdown(
        "PDR-studio is a browser-based control panel for **3D-PDR**. It lets you "
        "edit the compile-time flags and **compile** the code, set up the model "
        "parameters / initial conditions / chemistry, **run** the model with a "
        "live convergence bar, and **analyse** the outputs with interactive plots.\n\n"
        "The app lives *inside* your 3D-PDR checkout (`.../3D-PDR/PDR-studio/`) and "
        "finds the 3D-PDR root automatically as its own parent directory.")

    st.subheader("The sidebar at a glance")
    st.markdown(
        "| | Page | What you do there |\n"
        "|---|---|---|\n"
        "| 🏠 | **Home** | Check your environment: is the code built? which network? |\n"
        "| 📖 | **Tutorial** | This page. |\n"
        "| ⚙️ | **Code Configuration** | Edit `src/config.mk` flags and **compile**. |\n"
        "| 📐 | **Model Parameters** | 4-step wizard: parameters → abundances → network → **run**. |\n"
        "| 📡 | **Line emission** | Line tables, antenna temperature & optical depth, CO SLED. |\n"
        "| 📊 | **Abundance Profiles** | Tgas & abundances vs A_V / N_tot / n_H. |\n"
        "| 🔥 | **Heating & Cooling** | Heating/cooling functions, emissivities, level populations. |\n"
        "| 🧪 | **Chemical Analysis** | Formation/destruction pathways of one species per depth. |\n"
        "| 🧬 | **Network Builder** | Build a new chemical network with MakeRates and import it. |\n"
        "| 🐞 | **Report a bug** | Bundle everything into a tarball for the 3D-PDR team. |\n\n"
        "The pages are ordered as the natural workflow — top to bottom.")

    st.subheader("Launch the app")
    st.markdown("From the `PDR-studio` directory:")
    st.code("./run.sh", language="bash")
    st.markdown(
        "On the **first** run this creates a self-contained `.venv`, installs the "
        "requirements, and starts Streamlit at **http://localhost:8501**. "
        "To launch manually instead:")
    st.code(".venv/bin/streamlit run app.py", language="bash")
    st.info(
        "**Different layout?** PDR-studio assumes it sits inside the 3D-PDR tree. "
        "If yours is elsewhere, point it at the root with the `PDR_ROOT` "
        "environment variable before launching.", icon="📁")

    st.subheader("Check your environment (🏠 Home)")
    st.markdown(
        "The **Home** page is a quick health check: whether the 3D-PDR root is "
        "detected, whether the executable is **built** (and when), the compiled "
        "`DIMENSIONS` / `NETWORK`, and status badges for `config.mk`, `params.dat`, "
        "`chemfiles/`, `ics/` and `src/`. If the executable says **not built**, "
        "that's fine — you build it in the next step. The status row at the top of "
        "this page shows the same thing.")

# ===========================================================================
# Your first model
# ===========================================================================
with tab_model:
    st.markdown("The through-line is: **configure & compile → set parameters → "
                "run → analyse.**")

    st.header("1 · Configure & compile (⚙️ Code Configuration)")
    st.markdown(
        "This page edits the compile-time flags in `src/config.mk` and builds the code.\n\n"
        "1. The most-used flags are shown directly; the rest live under **⚙️ Advanced "
        "flags**. Toggles are on/off flags, dropdowns are multiple-choice, text boxes "
        "are free values, each with a help tooltip.\n"
        "2. Set the essentials — e.g. `DIMENSIONS`, `NETWORK`, and physics switches "
        "like `THERMALBALANCE`, `DUST`, `XRAYS`, `CRATTENUATION`, `SUPRATHERMAL`. "
        "**If you want to use Chemical Analysis later, enable `CHEMANALYSIS` now.**\n"
        "3. Leave **Clean build** ticked for the first build, then press "
        "**💾 Save & Compile**. The compile log streams live and collapses on success; "
        "PDR-studio also builds the **RT-tool** so line emission is ready.\n"
        "4. On success, press **Next ➡️** to go to Model Parameters.")
    st.info(
        "The flag list isn't hard-coded — PDR-studio reads it live from `config.mk` "
        "and `src/makefile`, so a newer 3D-PDR with extra flags just works. Brand-new "
        "flags are announced in a *“new flag(s) imported”* box under **Advanced**.",
        icon="🔎")

    st.header("2 · Set the model parameters (📐 Model Parameters)")
    st.markdown(
        "This page is a **4-step wizard**. Use the step control at the top, or the "
        "**💾 Save & Next ➡️** button to advance. The caption reminds you which "
        "compiled flags are active — some fields enable/disable based on them.")

    st.subheader("Step 1 — `params.dat`")
    st.markdown(
        "- **Input / Output**: the ICs directory (usually `ics`), the output "
        "directory (e.g. `sims`, created for you if missing), and a **model name** — "
        "the output prefix; every result file is `<outdir>/<prefix>.*`.\n"
        "- **Initial conditions** — the density grid 3D-PDR reads. Choose one:\n"
        "    - *Create ICs*: a uniform-density 1-D cloud, log-spaced in A_V "
        "(set n_H, A_V,max, log₁₀ A_V,min, resolution). Written as `<model>-ics.dat`.\n"
        "    - *Choose a pre-existing file*: use `ics/Aveff-nH.dat` as-is.\n"
        "    - *Import your own file*: upload a grid (columns `x[pc] y z nH`).\n"
        "- **PDR parameters**: FUV field strength (**Draine** or **Habing** units — "
        "Habing is converted on save), cosmic-ray ionisation rate ζ (or the Padovani "
        "L/H/U field / attenuation parameters, depending on `CRATTENUATION`), "
        "microturbulent velocity, and redshift (sets T_CMB).\n"
        "- **Thermal balance values**: min gas temperature, dust temperature, etc. "
        "If `THERMALBALANCE` is off, the “gas temperature” field becomes the fixed "
        "isothermal temperature.\n"
        "- **Coolant files**: CO, C⁺, C and O are always included; add optional "
        "coolants (filtered to species present in your network).\n"
        "- **Advanced** (optional): ODE tolerances, HEALPix level, iteration counts.\n\n"
        "Press **💾 Save & Next ➡️**.")

    st.subheader("Step 2 — Initial abundances")
    st.markdown(
        "Set the elemental abundances for the active network and the **dust-to-gas "
        "ratio** (1 = solar). **↩️ Reset to solar metallicity** restores defaults. "
        "Save & Next.")

    st.subheader("Step 3 — Chemical network")
    st.markdown(
        "Browse the reaction network. Filter by reactant(s), product(s), or "
        "non-species reactants (PHOTON, CRP, CRPHOT, X-ray processes…). To change a "
        "rate, pick a reaction and edit its **α, β, γ, Tmin, Tmax** — read at runtime, "
        "so the change takes effect on the next run **without recompiling**. "
        "This step is optional; press **Next ➡️** to skip.")

    st.subheader("Step 4 — Run model")
    st.markdown(
        "Press **▶️ Run 3D-PDR**. A progress bar tracks *“Thermal balance: X% "
        "converged”*; **⏹️ Stop** aborts. On success, PDR-studio automatically runs "
        "the RT-tool so line emission is ready, and **Next ➡️** takes you to the "
        "analysis pages. Outputs land in `<outdir>/<prefix>.*`.")

# ===========================================================================
# Analyse
# ===========================================================================
with tab_analyse:
    st.markdown(
        "All four analysis pages default to your most recently run model (you can "
        "pick a different output directory / prefix at the top). Plots are "
        "interactive: **drag to pan, scroll to zoom, double-click to reset, hover "
        "for values, click legend entries to toggle traces.** A dashed / dot-dashed "
        "vertical line marks the **HI→H₂** (and HI→2H₂, 50 % H₂) transition.")

    st.subheader("📡 Line emission")
    st.markdown(
        "RT-tool results in four tabs plus a scratchpad:\n"
        "- **Line tables**: per coolant, Line / Freq / Tex / τ / Tr at one depth "
        "point — pick it with the Ntot / Av / N(H₂) slider (an arrow marks the "
        "HI/H₂ transition). [CII] shows one line; C and O two; others up to ten.\n"
        "- **Antenna temperature** & **Optical depth**: how T_r / τ build up with "
        "column density, one panel per coolant, with a per-panel transition selector.\n"
        "- **CO SLED**: T_r of the first CO transitions (J→J−1) at a chosen column "
        "density, with optional normalisation to a transition of your choice.")

    st.subheader("📊 Abundance Profiles")
    st.markdown(
        "A 2×2 log-log figure vs **A_V**, **N_tot** or **n_H**: gas temperature; the "
        "HI/H₂ transition; C⁺/C/CO; and a **panel of species you choose** from the "
        "network.")

    st.subheader("🔥 Heating & Cooling")
    st.markdown(
        "Heating functions, cooling functions (both with a thick line for the total), "
        "per-coolant **emissivities**, and **level populations**, all vs A_V or n_H. "
        "Use the multiselects to pick which terms / transitions / levels to show.")

    st.subheader("🧪 Chemical Analysis")
    st.warning(
        "Requires the code to be built with the **`CHEMANALYSIS`** flag "
        "(⚙️ Code Configuration), then a fresh run — it reads "
        "`<prefix>.chemanalysis.fin`.", icon="⚠️")
    st.markdown(
        "Track **one** species: its abundance profile is plotted, and as you "
        "**hover** along the curve the panel below updates to that depth point, "
        "listing the dominant **formation (green)** and **destruction (red)** "
        "reactions with their percentage contributions and totals. Optionally press "
        "**🧠 Summarize … pathways** for a short natural-language digest "
        "(see *Going further*).")

    st.subheader("🖥️ Custom-plot scratchpad (every analysis page)")
    st.markdown(
        "Each analysis page has a **Custom plot** tab: write Python/matplotlib on the "
        "left, press **▶ Run**, and your figure renders on the right. The model's data "
        "is preloaded, e.g. `Av, nH, Tgas, Tdust, x, UV`, `abundance(\"CO\")`, "
        "`species`, and `np`/`plt`/`fig`/`ax`. On Line emission you also get "
        "`coolants`, `coolant_data(\"CO\")` (`.Tr[:, t]`, `.Nmol`, …) and "
        "`tau_data(...)`; on Heating & Cooling you get `total_heating`/"
        "`total_cooling`, `heating(...)`, `cooling(...)`, `emissivity(...)`, "
        "`level_pop(...)`. The starter code in each tab is a working example — edit "
        "it and re-run. **↺ Reset code** restores it.")

# ===========================================================================
# Going further
# ===========================================================================
with tab_further:
    st.header("Build a new chemical network (🧬 Network Builder)")
    st.markdown(
        "Create a network with the embedded, offline **MakeRates**:\n\n"
        "1. Upload a **Rate05-formatted (CSV)** master reaction file.\n"
        "2. Choose species either by **elements** (`e-, H, H2, He, C+, O` are always "
        "in; add `Mg+, N, S, Fe, …`) with complexity limits (max carbon / heavy atoms "
        "per species, anions on/off), or **import your own species file**.\n"
        "3. Name the network and press **⚙️ Build network** — this produces "
        "`species_<NAME>.d`, `rates_<NAME>.d`, `odes_<NAME>.c` (initial abundances "
        "zeroed). Download them, or…\n"
        "4. **📥 Import** into 3D-PDR: it copies the files into place and wires "
        "`NETWORK = <NAME>` into `config.mk`, the makefile and the relevant sources. "
        "**Recompile** to use it; its elements then appear in Model Parameters → "
        "Initial abundances.")

    st.header("Offline AI pathway summaries (optional)")
    st.markdown(
        "The **Chemical Analysis** summary button runs a **local, offline** LLM — the "
        "numbers come from the model, the LLM only phrases the deterministic digest. "
        "It's optional: the pathway table works without it. To enable:\n\n"
        "- **Ollama** (recommended): install it and `ollama pull llama3.2` (a ~2 GB "
        "3B model is plenty). Auto-detected afterwards; set `OLLAMA_HOST` if not on "
        "the default port.\n"
        "- **llama.cpp**: `pip install llama-cpp-python` into `.venv` and set "
        "`PDRSTUDIO_LLM_GGUF=/path/to/model.gguf` (`run.sh` already exports this var).")

    st.header("Report a bug (🐞)")
    st.markdown(
        "Describe the problem (≥ 50 characters) and press **📦 Build debug report**. "
        "PDR-studio bundles the config, parameters, chemistry files, your initial "
        "conditions and the current model's outputs into a single tarball to download "
        "and email to the 3D-PDR team.")

# ===========================================================================
# Troubleshooting
# ===========================================================================
with tab_help:
    st.header("Troubleshooting & tips")
    st.markdown(
        "- **“Executable: not built”** — compile first in ⚙️ Code Configuration.\n"
        "- **Analysis page says “no `*.pdr.fin` found”** — run a model first "
        "(📐 Model Parameters → Run), and check the output directory matches.\n"
        "- **Chemical Analysis is empty / warns about `.chemanalysis.fin`** — rebuild "
        "with `CHEMANALYSIS` enabled, then re-run the model.\n"
        "- **Line emission says RT-tool failed / not ready** — make sure the model has "
        "finished; use the **🔄 Retry RT-tool** button. RT-tool runs automatically "
        "right after a successful model run.\n"
        "- **Where do outputs go?** — `<outdir>/<prefix>.*`, with `outdir` and "
        "`prefix` set on Model Parameters → `params.dat`.\n"
        "- **Edited a rate but it “didn't take”** — reaction-rate edits apply on the "
        "next **run** (no recompile). Flag changes in Code Configuration **do** need a "
        "recompile.\n"
        "- **New 3D-PDR version?** — compile flags (`config.mk`) and parameters "
        "(`params.dat`) are discovered automatically; new entries show up under "
        "*Advanced* / *Imported parameters*, so PDR-studio keeps working.")
