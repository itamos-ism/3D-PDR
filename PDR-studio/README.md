# PDR-studio

A browser-based control panel for **3D-PDR**. Configure, compile, run, and
analyse PDR models from a friendly local web app — no terminal required.

## Quick start

```bash
cd PDR-studio
./run.sh
```

`run.sh` creates a self-contained `.venv` on first run and launches the app
(opens at http://localhost:8501). To launch manually:

```bash
.venv/bin/streamlit run app.py
```

To make PDR-studio reachable from other machines on your local network, use
`./run-lan.sh` instead — it prints the address to share (e.g.
`http://<your-ip>:8501`). Anyone who can reach that address gets the full app,
which can compile/run 3D-PDR and execute arbitrary Python typed into the
**Custom plot** tabs *on this machine*; only use it on a trusted network.

## Layout

```
PDR-studio/
  app.py                  # entry point: st.navigation (Home + feature pages)
  home.py                 # Home / welcome page
  pages/                  # one file per feature (titles/icons set in app.py)
    0_📖_Tutorial.py               # in-app walkthrough, sidebar top-to-bottom
    1_⚙️_Code_Configuration.py
    2_📐_Model_Parameters.py       # step-by-step wizard (Save & Next)
    3_📡_Line_emission.py          # runs RT-tool: line tables, Tr/τ build-up, CO SLED
    4_📊_Abundances_Profiles.py    # abundance/Tgas profiles + custom-plot tab
    5_🔥_Heating_and_Cooling.py    # heating/cooling/emissivities/level pops + custom-plot tab
    6_🧪_Chemical_Analysis.py      # reads <prefix>.chemanalysis.fin
    7_🧬_Network_Builder.py
    8_🐞_Report_a_bug.py
  pdrstudio/              # backend (one module per concern)
    paths.py          # locate the 3D-PDR tree
    config.py         # discover/read/write src/config.mk (+ makefile-aware flag schema)
    params.py         # read/write params.dat + per-network initial abundances
    network.py        # resolve species_<suffix>.d / rates_<suffix>.d (X-ray aware)
    networkbuilder.py # build a new network via MakeRates + import it
    chem.py           # browse/edit the reaction network (chemfiles/rates_*.d)
    chemanalysis.py   # read <prefix>.chemanalysis.fin (formation/destruction pathways)
    runner.py         # compile & run, streaming logs
    rttool.py         # drive & read the RT-tool post-processor
    outputs.py        # load <prefix>.pdr.fin model outputs
    fortread.py       # robust parsing of Fortran E/ES-format ASCII tables
    llm.py            # optional offline local-LLM backend (Ollama / llama.cpp)
    report.py         # build debug-report tarballs
    ui.py             # shared Streamlit helpers
  makerates/    # embedded MakeRates (MakeRates.py), run headless
  requirements.txt
  run.sh
  run-lan.sh
```

PDR-studio finds the 3D-PDR root as its own parent directory. Override with the
`PDR_ROOT` environment variable if needed.

## Feature status

| Feature | Status |
|---|---|
| 🏠 Home (detects whether 3D-PDR is built, current `config.mk` build, quick links) | ✅ |
| 📖 Tutorial (hands-on in-app walkthrough, follows the sidebar top to bottom) | ✅ |
| ⚙️ Code Configuration (edit `config.mk`, **Save & Compile**, Next → Model Parameters) | ✅ |
| 📐 Model Parameters (step wizard: `params.dat` → abundances+dust/gas → network → run → analysis) | ✅ |
| 📡 Line emission (runs RT-tool; line tables, Tr & τ build-up, CO SLED + custom-plot tab) | ✅ |
| 📊 Abundance Profiles (2×2 log-log vs A_V / n_H; + a Python custom-plot tab) | ✅ |
| 🔥 Heating & Cooling (heating, cooling, emissivities, level populations; + a custom-plot tab) | ✅ |
| 🧪 Chemical Analysis (per-species formation/destruction pathways, hover-linked to depth; optional local-LLM summary) | ✅ |
| 🧬 Network Builder (build a new network with MakeRates from selected elements, import it into 3D-PDR) | ✅ |
| 🐞 Report a bug (bundle config/params/ICs/outputs into a tarball to email the 3D-PDR team) | ✅ |

## Adapting to different 3D-PDR versions

The compile-flag list is **not** hard-coded. On every load, `config.py`
*discovers* the flags from the actual checkout, so a newer `config.mk` with new
flags works without changing PDR-studio:

- **`src/config.mk`** — the `KEY = VALUE` assignments (which flags exist, their
  values, the display order) and the comment header (per-flag help, enumerated
  values, and `#---- … options ----` section grouping).
- **`src/makefile`** — the `ifeq ($(KEY),VALUE)` branches give the full set of
  allowed values for each flag (often more than the comments list).

A small `CURATED` overlay in `config.py` adds nice labels/help and marks which
flags are *advanced* or *hidden*. A flag with no overlay entry still appears
automatically — under its documented section, in the **Advanced** area — and the
Code Configuration page reports it as *newly imported*. To promote such a flag
(custom label/grouping), add one entry to `CURATED`; nothing else changes.
Saving only rewrites `KEY = VALUE` lines, so comments, ordering and unknown keys
are preserved across versions.

The same idea applies to **`params.dat`** (`params.py`). It is a *positional*
file, so each value is matched by its inline `!label` (not its line number):
a newer version that inserts extra parameter lines (anywhere before the final
coolant block) does not shift the reads, new entries are shown under **Imported
parameters** on the Model Parameters page, and saving rewrites the file
*structure-preserving* — section order, comments and unknown entries are kept in
their original positions, with only the values substituted.

The **chemistry network files** are resolved in `network.py`: `species_<suffix>.d`
and `rates_<suffix>.d` always share one suffix, and the suffix follows the
compiled flags — an X-ray build (`XRAYS=1`) reads the `X` variant
(`species_Xreduced.d` / `rates_Xreduced.d`), matching the makefile. The
**Chemical network** tab also has a *Non-species reactants* filter that lists the
pseudo-reactants present in the active network (PHOTON, CRP, CRPHOT, the X-ray
`XRLYA`/`XRSEC`/`XRPHOT`, the grain `#`, …) so you can list every reaction driven
by a given non-collisional process.

## Local LLM pathway summaries (optional, offline)

The **Chemical Analysis** tab has a *"Summarize <species> pathways"* button that
writes a short natural-language account of how the selected species forms and is
destroyed along the depth/density column. It runs entirely **offline** on a local
model and is optional — without one, the deterministic pathway table still works.

The numbers come from the model, not the LLM: `chemanalysis.pathway_evidence`
condenses the per-point `*.chemanalysis.fin` data into the dominant pathways
(with Av ranges) plus depth snapshots, and the LLM (`llm.py`) only phrases that
digest — it's told to use only the reactions/numbers given. Backends, auto-detected:

- **[Ollama](https://ollama.com)** (recommended): install it, run `ollama pull llama3.2`
  (a ~2 GB 3B model is plenty); detected automatically afterwards, no internet needed.
- **llama.cpp**: `pip install llama-cpp-python` into `.venv` and set
  `PDRSTUDIO_LLM_GGUF=/path/to/model.gguf`.

Set `OLLAMA_HOST` if your Ollama server isn't on `http://localhost:11434`.

## Network Builder

The **Network Builder** tab creates a new chemical network with an embedded,
headless copy of **MakeRates** (`makerates/MakeRates.py`) and imports it into the
3D-PDR tree:

1. Upload a master reaction file (Rate05 CSV) and pick the elements to include
   (`e-`, `H`, `H2`, `He`, `C+`, `O` are always in; add any of
   `Mg+, N, S, Fe, Cl, P, Si, Na, F`). PDR-studio keeps the species made only of
   those elements — or you can import your own species file.
2. MakeRates builds `species_<NAME>.d`, `rates_<NAME>.d` and `odes_<NAME>.c`
   (all initial abundances set to `0.00E+00`).
3. **Import** copies species/rates → `chemfiles/`, odes → `src/`, and wires
   `NETWORK = <NAME>` into `config.mk`, the `makefile`, and the suffix `#ifdef`
   blocks of `read_species.F90` / `read_rates.F90` / `input_parameters.F90`
   (idempotent). Chemical heating needs no per-network edit: `CHEMICAL_HEATING`
   in `heatingfunctions.F90` calls the network-independent
   `chemical_heating_module` (`src/chemical_heating.F90`), which tags the
   relevant exothermic reactions automatically by matching reactants/products
   against a canonical table — so it is untouched. Recompile to use the new
   network; its elements then appear in the Model Parameters → Initial
   abundances tab.

> The ODEs are emitted as C (`odes_<NAME>.c`) because 3D-PDR compiles them via
> the makefile's `%.o: %.c` rule (the proposal's `.d` would not compile).

## Adding a feature

Drop a new `pages/N_Name.py` file; it appears in the sidebar automatically.
Put reusable logic in a `pdrstudio/` module and import it. Each page should
start with `ui.configure_page("Title")`.
