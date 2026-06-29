"""Feature 8 — Network Builder.

Build a new chemical network for 3D-PDR with MakeRates (embedded, run offline):

1. pick the elements to include (e-, H, H2, He, C+, O always in);
2. PDR-studio collects the species made only of those elements (or import your
   own species file);
3. MakeRates builds species_<NAME>.d / rates_<NAME>.d / odes_<NAME>.c
   (abundances zeroed);
4. import them into the 3D-PDR tree and wire up NETWORK=<NAME>.
"""
from __future__ import annotations

from pathlib import Path

import streamlit as st

from pdrstudio import networkbuilder as nb
from pdrstudio import ui

ui.configure_page("Network Builder")
st.title("🧬 Network Builder")
st.markdown(
    "Build a new chemical network for **3D-PDR** with MakeRates, then import it. "
    "Everything runs locally.")

# Hide the file uploader's "Limit 200MB per file" hint — these files are tiny (<1MB).
# That hint is the whole stFileUploaderDropzoneInstructions element (the "Drag and
# drop" text is a separate node), so hiding it removes only the size line.
st.markdown(
    "<style>[data-testid='stFileUploaderDropzoneInstructions']{display:none !important;}</style>",
    unsafe_allow_html=True)

if nb.makerates_script() is None:
    st.error("MakeRates was not found. Embed it at `PDR-studio/makerates/MakeRates.py` "
             "(copy of `MakeRates_python312.py`).")
    st.stop()

work = nb.staging_dir("_input")

# ===========================================================================
# 1 · Master reaction file
# ===========================================================================
st.subheader("1 · Master reaction file")
up = st.file_uploader("Reaction file — a Rate05-formatted (CSV) network is needed", key="nb_rxn")
if up is not None:
    rxn_path = work / "reactionFile.dat"
    rxn_path.write_bytes(up.getvalue())
    st.session_state["nb_rxn_path"] = str(rxn_path)

if not st.session_state.get("nb_rxn_path"):
    st.info("Upload a master reaction file to begin.", icon="⬆️")
    st.stop()

rxn_path = Path(st.session_state["nb_rxn_path"])
rxn_text = rxn_path.read_text(errors="replace")
rate05 = nb.is_rate05_reaction_file(rxn_text)
all_species = nb.species_in_reaction_file(rxn_text) if rate05 else []
if rate05:
    st.caption(f"`{up.name if up else rxn_path.name}` · Rate05 ✓ · "
               f"{len(all_species)} species detected.")
else:
    st.warning("This file is not recognised as Rate05 CSV. Element-based selection "
               "needs Rate05 — you can still import your own species file below.", icon="⚠️")

# ===========================================================================
# 2 · Species selection
# ===========================================================================
st.subheader("2 · Species selection")
mode = st.radio("How to choose species", ["Build from elements", "Import my own species file"],
                horizontal=True)

species_path: Path | None = None
n_selected = 0

if mode == "Build from elements":
    st.markdown("**Always included** (cannot be deselected): `e-`, `H`, `H2`, `He`, `C+`, `O`")
    chosen = st.multiselect("Additional elements to include",
                            options=list(nb.OPTIONAL_ELEMENTS.keys()),
                            help="Species containing only the selected elements are kept.")
    elements = set(nb.MANDATORY_ELEMENTS) | {nb.OPTIONAL_ELEMENTS[c] for c in chosen}

    st.markdown("**Set the complexity of the network by limiting the number of carbon "
                "and heavy atoms per species.** "
                "Species with unmodelled elements (Ar, Ti, Al, …) are always excluded.")
    lc1, lc2, lc3 = st.columns(3)
    with lc1:
        max_carbon = st.slider("Max carbon atoms per species", 1, 8, 1,
                               help="The main lever against carbon chains. 1 ≈ reduced-style "
                                    "(keeps CO, HCO+… but no C2/C3/C4 or large organics).")
    with lc2:
        max_heavy = st.slider("Max heavy (non-H) atoms per species", 1, 10, 2,
                              help="Caps overall molecule size, e.g. 2 keeps CO/O2/HCO+ but "
                                   "drops bigger species.")
    with lc3:
        include_anions = st.toggle("Include anions (e.g. C⁻)", value=False,
                                   help="Negatively-charged species. Off by default "
                                        "(reduced-style). Cations are always included.")

    if rate05 and all_species:
        auto = nb.select_species_by_elements(
            all_species, elements, max_heavy=max_heavy, max_carbon=max_carbon,
            include_anions=include_anions)
        st.caption(f"Auto-selected **{len(auto)}** species (elements ⊆ "
                   f"{{{', '.join(sorted(elements))}}}, ≤{max_carbon} C, ≤{max_heavy} heavy"
                   f"{', anions on' if include_anions else ''}). Untick any you don't want; "
                   "`e-` is added automatically.")

        if st.button("↩️ Re-include all", help="Restore every auto-selected species."):
            for sp in auto:
                st.session_state[f"nb_keep_{sp}"] = True
            st.rerun()

        # grouped, editable table: one column per total-atom count
        from collections import defaultdict
        groups: dict[int, list[str]] = defaultdict(list)
        for sp in auto:
            groups[nb.total_atoms(sp)].append(sp)
        ncounts = sorted(groups)
        selected: list[str] = []
        if ncounts:
            tcols = st.columns(len(ncounts))
            for col, n in zip(tcols, ncounts):
                col.markdown(f"**{n} atom{'s' if n != 1 else ''}**")
                for sp in sorted(groups[n]):
                    key = f"nb_keep_{sp}"
                    st.session_state.setdefault(key, True)
                    if col.checkbox(sp, key=key):
                        selected.append(sp)

        n_selected = len(selected)
        species_path = work / "speciesFile.dat"
        nb.write_species_file(species_path, selected)
        st.success(f"**{n_selected}** of {len(auto)} species selected for the network "
                   "(all initial abundances will be 0).")
    else:
        st.info("Provide a Rate05 reaction file, or switch to importing your own species file.")
else:
    spc_up = st.file_uploader("Species file (Rate05 / Rate95 / Rate99)", key="nb_spc")
    if spc_up is not None:
        species_path = work / "speciesFile.dat"
        species_path.write_bytes(spc_up.getvalue())
        st.caption(f"Using `{spc_up.name}`.")

# ===========================================================================
# 3 · Build
# ===========================================================================
st.subheader("3 · Build the network")
name = st.text_input("Network name (output prefix → NETWORK)", value="NEWNETWORK",
                     help="Used for species_<NAME>.d / rates_<NAME>.d / odes_<NAME>.c "
                          "and the NETWORK flag. Letters, digits, underscore.")
name_ok = nb.valid_name(name)
if name and not name_ok:
    st.warning("Name must start with a letter and contain only letters, digits and underscores.")

if st.button("⚙️ Build network", type="primary", disabled=not (species_path and name_ok)):
    try:
        with st.spinner(f"Running MakeRates to build “{name}”…"):
            st.session_state["nb_result"] = nb.run_makerates(rxn_path, species_path, name)
    except Exception as exc:  # noqa: BLE001 — surface MakeRates errors to the user
        st.error(f"MakeRates failed: {exc}")
        st.session_state.pop("nb_result", None)

res = st.session_state.get("nb_result")
if res and res.name == name:
    st.success(f"Built **{name}** — {res.n_species} species · {res.n_reactions} reactions "
               "(initial abundances set to 0.00E+00).")
    d1, d2, d3 = st.columns(3)
    d1.download_button(f"⬇️ species_{name}.d", res.species_file.read_bytes(),
                       file_name=f"species_{name}.d", key="dl_sp")
    d2.download_button(f"⬇️ rates_{name}.d", res.rates_file.read_bytes(),
                       file_name=f"rates_{name}.d", key="dl_rt")
    d3.download_button(f"⬇️ odes_{name}.c", res.odes_file.read_bytes(),
                       file_name=f"odes_{name}.c", key="dl_od")
    with st.expander("MakeRates log"):
        st.code(res.log or "(no output)", language="text")

    # =======================================================================
    # 4 · Import into 3D-PDR
    # =======================================================================
    st.subheader("4 · Import into 3D-PDR")
    st.caption(
        f"Copies `species_{name}.d` and `rates_{name}.d` → `chemfiles/`, "
        f"`odes_{name}.c` → `src/`, and wires `NETWORK = {name}` into `config.mk`, "
        "the `makefile`, the suffix `#ifdef` blocks of `read_species.F90`, "
        "`read_rates.F90` and `input_parameters.F90`, the RTspop.fin network-name "
        "block of `src/writeoutputs.F90`, and the `rt_chemsuf` case in "
        "`RT-tool/readparams.F90`. Chemical heating is handled automatically by "
        "`AUTOCHEMHEAT`, so `heatingfunctions.F90` is left untouched. "
        "Then recompile in ⚙️ Code Configuration.")
    if st.button(f"📥 Import “{name}” into 3D-PDR", type="primary"):
        try:
            rep = nb.import_network(res, set_active=True)
            st.success(f"Imported **{name}**. Recompile the code to use it "
                       "(⚙️ Code Configuration → Compile).")
            if rep.copied:
                st.markdown("**Copied:** " + ", ".join(f"`{p}`" for p in rep.copied))
            if rep.patched:
                st.markdown("**Updated:** " + ", ".join(f"`{p}`" for p in rep.patched))
            if rep.skipped:
                st.markdown("**Skipped:** " + "; ".join(rep.skipped))
            for w in rep.warnings:
                st.warning(w)
        except Exception as exc:  # noqa: BLE001
            st.error(f"Import failed: {exc}")
