"""Feature — Model Parameters (a step-by-step wizard).

Steps (a segmented control at the top; "Save & Next" / "Next" buttons advance):
  params.dat  →  Initial abundances  →  Chemical network  →  Run model  →  Analysis I

All numeric fields are plain text inputs (no +/- steppers). The dust-to-gas ratio
lives in the Initial abundances step (default 1 = solar).
"""
from __future__ import annotations

import re

import streamlit as st

import pandas as pd

from pdrstudio import chem, config, params, paths, rttool, runner, ui
from pdrstudio import network as netmod

ui.configure_page("Model Parameters")

# Hide the file uploader's "Limit 200MB per file" hint (the whole instructions
# element; the "Drag and drop" text is a separate node) — ICs files are tiny.
st.markdown(
    "<style>[data-testid='stFileUploaderDropzoneInstructions']{display:none !important;}</style>",
    unsafe_allow_html=True)

# --- compiled-mode context --------------------------------------------------
cfg = config.read_config() if paths.config_mk().exists() else {}
network = cfg.get("NETWORK", "REDUCED")
net_suffix = netmod.suffix_for(network, netmod.xrays_enabled(cfg))
cratten = int(cfg.get("CRATTENUATION", "0") or 0)
suprathermal = cfg.get("SUPRATHERMAL", "0").strip() == "1"
thermalbalance = cfg.get("THERMALBALANCE", "1").strip() == "1"
dust_flag = cfg.get("DUST", "HTT91").strip()
isothermal = not thermalbalance

cur = params.read_params(cratten)
ss = st.session_state

# --- section (wizard) state -------------------------------------------------
S_PARAMS, S_ABUND, S_NET, S_RUN = "params.dat", "Initial abundances", "Chemical network", "Run model"
SECTIONS = [S_PARAMS, S_ABUND, S_NET, S_RUN]
if "_mp_goto" in ss:                       # pending nav from a Save & Next click
    ss["mp_section"] = ss.pop("_mp_goto")
ss.setdefault("mp_section", S_PARAMS)

st.title("📐 Model Parameters")
st.segmented_control("Step", SECTIONS, key="mp_section", label_visibility="collapsed")
section = ss["mp_section"]

st.caption(
    f"Compiled mode → NETWORK **{network}** · CRATTENUATION **{cratten}** · "
    f"SUPRATHERMAL **{'on' if suprathermal else 'off'}** · "
    f"THERMALBALANCE **{'on' if thermalbalance else 'off'}** · DUST **{dust_flag}**  "
    f"(set in ⚙️ Code Configuration, then recompile)")


def _go(target: str):
    ss["_mp_goto"] = target
    st.rerun()


# scientific-notation vs integer fields (for formatting defaults)
_SCI = {"zeta", "tend", "grain_radius", "av_fac", "rel_tol", "abs_tol", "cr_n0"}
_INT = {"healpix_level", "init_iter", "max_iter"}


def _fmt(key, value) -> str:
    if key in _SCI:
        return f"{float(value):.3E}"
    if key in _INT:
        return str(int(float(value)))
    return f"{float(value):g}"


# ===========================================================================
# Step 1 — params.dat
# ===========================================================================
if section == S_PARAMS:
    new: dict = {"d2g": cur["d2g"], "input": cur["input"]}  # d2g edited later; input set by ICs choice
    extras = cur.get("extras", [])

    # --- top Save & Next button (above all content) ---------------------------
    _, _tc = st.columns([4, 1])
    top_submit = _tc.button("💾 Save & Next ➡️", key="params_save_top",
                            type="primary", width="stretch")
    msg_slot = st.container()
    st.subheader("Input / Output")
    c1, c2 = st.columns(2)
    with c1:
        indir = st.text_input("ICs directory", value=str(cur["indir"]),
                              help="Directory with initial-condition files (usually `ics`).")
        outdir = st.text_input("Output directory", value=str(cur["outdir"]),
                               help="Where 3D-PDR writes outputs (e.g. `sims`).")
        create_missing = st.checkbox("Create if it does not exist", value=True,
                                     help="Create the output directory automatically when saving.")
    with c2:
        prefix = st.text_input("Name your model (output prefix)", value=str(cur["output"]),
                               help="Prefix for all output files: `<outdir>/<prefix>.*`")
    new.update(indir=indir, outdir=outdir, output=prefix)

    model_name = (prefix or "model").strip()
    st.markdown("**Initial conditions** — the density grid 3D-PDR reads. Written when you "
                f"press *Save & Next* (as `{model_name}-ics.dat`, except the pre-existing file).")
    ics_nH = ics_avmax = ics_avmin = ics_res = ics_upload = None
    # Restore ICS mode from params.dat when entering a fresh session so the
    # radio doesn't revert to "Create ICs" after the user navigates away.
    if "ics_mode" not in ss:
        ss["ics_mode"] = ("Choose a pre-existing file"
                          if str(cur.get("input", "")) == "Aveff-nH.dat" else "Create ICs")
    ics_mode = st.radio(
        "Initial conditions", ["Create ICs", "Choose a pre-existing file", "Import your own file"],
        horizontal=True, key="ics_mode",
        help="Create a uniform 1-D grid, use ics/Aveff-nH.dat, or upload your own file.")
    if ics_mode == "Create ICs":
        st.caption("A uniform-density 1-D cloud, log-spaced in A_V (port of `ics/uniform1D.f90`).")
        g1, g2, g3, g4 = st.columns(4)
        ics_nH = g1.text_input("Density n_H [cm⁻³]", value="1e3")
        ics_avmax = g2.text_input("A_V,max (linear)", value="20")
        ics_avmin = g3.text_input("log₁₀(A_V,min)", value="-3")
        ics_res = g4.text_input("Resolution (points/dex)", value="30")
    elif ics_mode == "Choose a pre-existing file":
        st.selectbox("Pre-existing ICs file", ["Aveff-nH.dat"], key="ics_pre",
                     help="Used as-is (not renamed).")
    else:
        ics_upload = st.file_uploader("Import an ICs density file", key="ics_upload",
                                      help="Columns: x[pc] y z nH. Saved as <model>-ics.dat.")
    st.divider()
    st.subheader("PDR parameters")
    c1, c2 = st.columns(2)
    with c1:
        fuv_unit = st.radio("FUV field unit", ["Draine (1974)", "Habing (1964)"],
                            help="Habing values are converted to Draine on save (÷1.7).")
        fuv_str = st.text_input("FUV strength", value=_fmt("fuv", cur["fuv"]))
        if cratten == 0:
            zeta_str = st.text_input("Cosmic-ray ionisation rate ζ (s⁻¹)",
                                     value=_fmt("zeta", cur["zeta"]),
                                     help="CRATTENUATION=0: a number. Default 1e-17.")
        elif cratten == 1:
            opts = ["L", "H", "U"]
            cf = str(cur.get("cr_field", "L")).upper()
            new["cr_field"] = st.selectbox("Cosmic-ray field (Padovani L/H/U)", opts,
                                           index=opts.index(cf) if cf in opts else 0)
        else:
            st.markdown("**CR attenuation (CRATTENUATION=2)**")
            cr_norm_str = st.text_input("norm", value=_fmt("cr_norm", cur["cr_norm"]))
            cr_n0_str = st.text_input("n0", value=_fmt("cr_n0", cur["cr_n0"]))
            cr_slope_str = st.text_input("slope", value=_fmt("cr_slope", cur["cr_slope"]))
    with c2:
        vturb_str = st.text_input("Microturbulent velocity (km/s)", value=_fmt("vturb", cur["vturb"]))
        redshift_str = st.text_input("Redshift (sets T_CMB)", value=_fmt("redshift", cur["redshift"]))

    st.markdown("**Suprathermal CO formation** "
                + ("(active)" if suprathermal else "*(inactive — SUPRATHERMAL=0; values still stored)*"))
    c1, c2 = st.columns(2)
    with c1:
        avcrit_str = st.text_input("Av_crit", value=_fmt("av_crit", cur["av_crit"]),
                                   disabled=not suprathermal)
    with c2:
        valfv_str = st.text_input("Alfvén velocity (km/s)", value=_fmt("v_alfv", cur["v_alfv"]),
                                  disabled=not suprathermal)
    st.divider()

    st.subheader("Thermal balance values")
    if isothermal:
        st.info("THERMALBALANCE is **off**: the **gas temperature** below is the fixed "
                "isothermal temperature.", icon="🌡️")
    c1, c2, c3 = st.columns(3)
    with c1:
        tgas_str = st.text_input("Gas temperature (K)", value=_fmt("tgas", cur["tgas"]),
                                 disabled=thermalbalance,
                                 help="Isothermal temperature when THERMALBALANCE=0.")
    with c2:
        tmin_str = st.text_input("Min gas temperature (K)", value=_fmt("tmin", cur["tmin"]))
    with c3:
        tdust_str = st.text_input("Dust temperature (K)", value=_fmt("tdust", cur["tdust"]),
                                  disabled=dust_flag != "0",
                                  help="Editable only when DUST=0 (isothermal dust).")
    st.divider()

    st.subheader("Coolant files")
    st.caption("Always included: " + " · ".join(f"`{lab}`" for lab, _ in params.FIXED_COOLANTS))
    avail = params.optional_for_network(network)
    label_to_file = {c.label: c.filename for c in avail}
    file_to_label = {c.filename: c.label for c in avail}
    default_labels = [file_to_label[f] for f in cur["coolants"] if f in file_to_label]
    chosen = st.multiselect(
        f"Additional coolants available for the {network} network",
        options=list(label_to_file.keys()), default=default_labels,
        help="Filtered to coolants whose species exists in this network.")
    new["coolants"] = [label_to_file[lbl] for lbl in chosen]
    st.divider()

    with st.expander("⚙️ Advanced parameters (change only if you know what you are doing)"):
        st.markdown("**PDR parameters**")
        c1, c2 = st.columns(2)
        with c1:
            tend_str = st.text_input("Chemical evolution time (yr)", value=_fmt("tend", cur["tend"]))
            grain_str = st.text_input("Grain radius (cm)", value=_fmt("grain_radius", cur["grain_radius"]))
        with c2:
            avfac_str = st.text_input("N(H)→Av conversion (mag cm²)", value=_fmt("av_fac", cur["av_fac"]))
            uvfac_str = st.text_input("UV factor", value=_fmt("uv_fac", cur["uv_fac"]))
        st.markdown("**Thermal balance**")
        c1, c2, c3 = st.columns(3)
        tmax_str = c1.text_input("Max gas temperature (K)", value=_fmt("tmax", cur["tmax"]))
        fcrit_str = c2.text_input("Convergence Fcrit", value=_fmt("fcrit", cur["fcrit"]))
        tdiff_str = c3.text_input("Convergence Tdiff", value=_fmt("tdiff", cur["tdiff"]))
        st.markdown("**Ray-tracing · ODE solver · Chemistry iterations**")
        c1, c2, c3 = st.columns(3)
        with c1:
            hpx_str = c1.text_input("HEALPix level", value=_fmt("healpix_level", cur["healpix_level"]))
            theta_str = c1.text_input("theta_crit", value=_fmt("theta_crit", cur["theta_crit"]))
        with c2:
            reltol_str = c2.text_input("Relative tolerance", value=_fmt("rel_tol", cur["rel_tol"]))
            abstol_str = c2.text_input("Absolute tolerance", value=_fmt("abs_tol", cur["abs_tol"]))
        with c3:
            init_str = c3.text_input("Initial iterations", value=_fmt("init_iter", cur["init_iter"]))
            maxit_str = c3.text_input("Max iterations", value=_fmt("max_iter", cur["max_iter"]))

    extra_inputs: dict[str, str] = {}
    if extras:
        with st.expander(f"🆕 Imported parameters ({len(extras)})", expanded=True):
            st.caption("Entries found in `params.dat` with no curated field; saved in place.")
            ecols = st.columns(2)
            for i, ex in enumerate(extras):
                extra_inputs[ex.key] = ecols[i % 2].text_input(
                    ex.comment, value=str(ex.value), key=f"px_{ex.key}", help=f"Section: {ex.section}")

    # bottom-right Save & Next
    _, bc = st.columns([4, 1])
    bottom_submit = bc.button("💾 Save & Next ➡️", key="params_save_bottom",
                              type="primary", width="stretch")

    with st.expander("View raw params.dat"):
        if paths.params_dat().exists():
            st.code(paths.params_dat().read_text(), language="text")

    if top_submit or bottom_submit:
        errors: list[str] = []

        def pf(label, s, key, typ=float):
            try:
                new[key] = typ(float(s)) if typ is int else float(s)
            except (ValueError, TypeError):
                errors.append(f"{label}: '{s}' is not a valid number.")

        fuv_v = None
        try:
            fuv_v = float(fuv_str)
        except ValueError:
            errors.append(f"FUV strength: '{fuv_str}' is not a valid number.")
        if fuv_v is not None:
            new["fuv"] = fuv_v / 1.7 if fuv_unit.startswith("Habing") else fuv_v
        pf("Microturbulent velocity", vturb_str, "vturb")
        pf("Redshift", redshift_str, "redshift")
        pf("Av_crit", avcrit_str, "av_crit")
        pf("Alfvén velocity", valfv_str, "v_alfv")
        pf("Gas temperature", tgas_str, "tgas")
        pf("Min gas temperature", tmin_str, "tmin")
        pf("Dust temperature", tdust_str, "tdust")
        pf("Chemical evolution time", tend_str, "tend")
        pf("Grain radius", grain_str, "grain_radius")
        pf("N(H)→Av conversion", avfac_str, "av_fac")
        pf("UV factor", uvfac_str, "uv_fac")
        pf("Max gas temperature", tmax_str, "tmax")
        pf("Convergence Fcrit", fcrit_str, "fcrit")
        pf("Convergence Tdiff", tdiff_str, "tdiff")
        pf("theta_crit", theta_str, "theta_crit")
        pf("Relative tolerance", reltol_str, "rel_tol")
        pf("Absolute tolerance", abstol_str, "abs_tol")
        pf("HEALPix level", hpx_str, "healpix_level", int)
        pf("Initial iterations", init_str, "init_iter", int)
        pf("Max iterations", maxit_str, "max_iter", int)
        if cratten == 0:
            pf("Cosmic-ray ζ", zeta_str, "zeta")
        elif cratten == 2:
            pf("CR norm", cr_norm_str, "cr_norm")
            pf("CR n0", cr_n0_str, "cr_n0")
            pf("CR slope", cr_slope_str, "cr_slope")
        new.update(extra_inputs)

        # --- initial conditions (validated here, written below on success) ----
        ics_dir = paths.pdr_root() / new["indir"]
        ics_payload = None    # (data: str|bytes, filename) to write on success
        if ics_mode == "Create ICs":
            try:
                grid = params.uniform_ics_grid(float(ics_nH), float(ics_avmax),
                                               float(ics_avmin), int(float(ics_res)))
                ics_payload = (grid, f"{model_name}-ics.dat")
            except (ValueError, TypeError):
                errors.append("Create ICs: n_H, A_V,max, log(A_V,min) and resolution "
                              "must all be numbers.")
        elif ics_mode == "Choose a pre-existing file":
            new["input"] = "Aveff-nH.dat"
        else:  # Import your own file
            if ics_upload is None:
                errors.append("Import your own file: please upload an ICs density file.")
            else:
                ics_payload = (ics_upload.getvalue(), f"{model_name}-ics.dat")

        if errors:
            with msg_slot:
                for e in errors:
                    st.error(e)
        else:
            outpath = paths.pdr_root() / new["outdir"]
            if create_missing and not outpath.exists():
                outpath.mkdir(parents=True, exist_ok=True)
            if ics_payload is not None:
                data, fname = ics_payload
                ics_dir.mkdir(parents=True, exist_ok=True)
                target = ics_dir / fname
                target.write_bytes(data) if isinstance(data, bytes) else target.write_text(data)
                new["input"] = fname
            params.write_params(new, cratten)
            _go(S_ABUND)

# ===========================================================================
# Step 2 — initial abundances (+ dust-to-gas)
# ===========================================================================
elif section == S_ABUND:
    st.markdown(f"Set the initial elemental abundances for the **{network}** network "
                f"(`chemfiles/species_{net_suffix}`) and the dust-to-gas ratio.")
    current_ab = params.read_abundances(network)
    if not current_ab:
        st.error(f"Could not read species file for network {network}.")
    else:
        def _abkey(token: str) -> str:
            return f"ab_{network}_{token}"

        for token, val in current_ab.items():
            ss.setdefault(_abkey(token), f"{float(val):.2E}")
        ss.setdefault("ab_d2g", _fmt("d2g", cur["d2g"]))

        # top-right Save & Next
        _, tc = st.columns([4, 1])
        top_next = tc.button("💾 Save & Next ➡️", key="ab_next_top", type="primary", width="stretch")

        if st.button("↩️ Reset to solar metallicity"):
            for token in current_ab:
                ss[_abkey(token)] = f"{float(params.ELEMENT_DEFAULTS.get(token, current_ab[token])):.2E}"
            ss["ab_d2g"] = "1"
            st.rerun()

        d2g_str = st.text_input("Dust-to-gas ratio (norm. to 1e-2)", key="ab_d2g",
                                help="1 = solar metallicity. Lower/raise for sub-/super-solar.")

        ab_inputs: dict[str, str] = {}
        cols = st.columns(min(3, len(current_ab)))
        for i, token in enumerate(current_ab):
            ab_inputs[token] = cols[i % len(cols)].text_input(token, key=_abkey(token))

        _, bc = st.columns([4, 1])
        bottom_next = bc.button("💾 Save & Next ➡️", key="ab_next_bottom", type="primary", width="stretch")

        if top_next or bottom_next:
            parsed: dict[str, float] = {}
            ok = True
            for token, s in ab_inputs.items():
                try:
                    parsed[token] = float(s)
                except ValueError:
                    st.error(f"{token}: '{s}' is not a valid number."); ok = False
            try:
                d2g_val = float(d2g_str)
            except ValueError:
                st.error(f"Dust-to-gas ratio: '{d2g_str}' is not a valid number."); ok = False
            if ok:
                params.write_abundances(network, parsed)
                p = params.read_params(cratten)   # update only d2g in params.dat
                p["d2g"] = d2g_val
                params.write_params(p, cratten)
                _go(S_NET)

# ===========================================================================
# Step 3 — chemical network
# ===========================================================================
elif section == S_NET:
    _, tc = st.columns([4, 1])
    if tc.button("Next ➡️", key="net_next_top", type="primary", width="stretch"):
        _go(S_RUN)

    st.markdown(f"Browse the **{network}** reaction network (`chemfiles/rates_{net_suffix}`). "
                "Search by reactants/products, then edit a reaction's rate coefficients.")
    reactions = chem.read_reactions(network)
    if not reactions:
        st.error(f"Could not read the rates file for network {network}.")
    else:
        species = netmod.species_list(net_suffix)
        nonspecies = chem.nonspecies_reactants(reactions, species)
        species_only = [t for t in chem.species_tokens(reactions) if t not in set(nonspecies)]
        st.caption(f"{len(reactions)} reactions · {len(species_only)} species · "
                   f"{len(nonspecies)} non-species reactants.")
        c1, c2 = st.columns(2)
        sel_reac = c1.multiselect("Reactant(s)", species_only, max_selections=2)
        sel_prod = c2.multiselect("Product(s)", species_only, max_selections=2)
        sel_nonspec = st.multiselect("Non-species reactants", nonspecies,
                                     help="PHOTON, CRP, CRPHOT, X-ray XRLYA/XRSEC/XRPHOT, …")
        matches = chem.search_reactions(reactions, reactants=sel_reac + sel_nonspec, products=sel_prod)
        st.write(f"**{len(matches)}** matching reaction(s)"
                 + ("" if (sel_reac or sel_prod or sel_nonspec) else " (whole network)"))
        if matches:
            df = pd.DataFrame([{
                "#": r.index, "reaction": r.equation(), "α": f"{r.alpha:.2E}",
                "β": f"{r.beta:.2f}", "γ": f"{r.gamma:.1f}",
                "Tmin": f"{r.tmin:g}", "Tmax": f"{r.tmax:g}",
            } for r in matches])
            st.dataframe(df, use_container_width=True, hide_index=True, height=300)

            st.divider()
            st.subheader("Edit a reaction")
            options = {r.line_no: r.label() for r in matches}
            pick = st.selectbox("Reaction to edit", list(options.keys()),
                                format_func=lambda ln: options[ln])
            chosen = next(r for r in matches if r.line_no == pick)
            with st.form("rate_edit_form"):
                st.markdown(f"**{chosen.equation()}**  ·  reaction #{chosen.index}  ·  source `{chosen.clem}`")
                e1, e2, e3 = st.columns(3)
                a_str = e1.text_input("α (alpha)", value=f"{chosen.alpha:.2E}")
                b_str = e2.text_input("β (beta)", value=f"{chosen.beta:.2f}")
                g_str = e3.text_input("γ (gamma, K)", value=f"{chosen.gamma:.1f}")
                t1, t2 = st.columns(2)
                tmin_str = t1.text_input("Tmin (K)", value=f"{chosen.tmin:g}")
                tmax_str = t2.text_input("Tmax (K)", value=f"{chosen.tmax:g}")
                save_rate = st.form_submit_button("💾 Save reaction", type="primary")
            if save_rate:
                try:
                    chem.update_reaction(network, chosen.line_no, alpha=float(a_str),
                                         beta=float(b_str), gamma=float(g_str),
                                         tmin=float(tmin_str), tmax=float(tmax_str))
                    st.success(f"Reaction #{chosen.index} updated in rates_{net_suffix}. "
                               "Read at runtime — effective next run, no recompile.")
                except ValueError as exc:
                    st.error(f"Invalid number: {exc}")

    _, bc = st.columns([4, 1])
    if bc.button("Next ➡️", key="net_next_bottom", type="primary", width="stretch"):
        _go(S_RUN)

# ===========================================================================
# Step 4 — run model
# ===========================================================================
elif section == S_RUN:
    exe = paths.executable()

    def _go_analysis():
        st.switch_page("pages/3_📡_Line_emission.py")

    bp = ss.get("run_proc")
    running = bp is not None and bp.is_running()
    finished_ok = bp is not None and not running and bp.returncode == 0

    # top-right Next (only after a successful run)
    _, tc = st.columns([4, 1])
    if finished_ok and tc.button("Next ➡️", key="run_next_top", type="primary", width="stretch"):
        _go_analysis()

    st.markdown("Run the compiled **./3DPDR** model with the current `params.dat`.")
    st.caption(f"Executable: `{exe}` · {ui.status_badge(exe.exists(), 'built', 'not built')} "
               f"· last built {ui.mtime_str(exe)}")
    if not exe.exists():
        st.warning("The code is not compiled. Build it in ⚙️ Code Configuration first.", icon="⚠️")

    b1, b2 = st.columns(2)
    run_clicked = b1.button("▶️ Run 3D-PDR", type="primary",
                            disabled=not exe.exists() or running, width="stretch")
    stop_clicked = b2.button("⏹️ Stop 3D-PDR", disabled=not running, width="stretch")

    if run_clicked:
        proc = runner.model_process()
        if proc is None:
            st.error("Executable not found — compile the code first.")
        else:
            proc.start()
            ss["run_proc"] = proc
            # Invalidate any previous RT-tool results so they re-run after this model.
            ss.pop("rt_prefix", None)
            ss.pop("rt_stdout", None)
            ss.pop("rt_failed_prefix", None)
            st.rerun()
    if stop_clicked and bp is not None:
        bp.stop()
        st.rerun()

    _TB_RE = re.compile(r"Thermal balance is\s+([\d.]+)%\s+converged", re.IGNORECASE)

    def _parse_progress(lines: list[str]) -> float | None:
        for line in reversed(lines):
            m = _TB_RE.search(line)
            if m:
                return float(m.group(1))
        return None

    @st.fragment(run_every=1.0)
    def _live_progress():
        p = ss.get("run_proc")
        if p is None:
            return
        p.drain()
        pct = _parse_progress(p.lines)
        if pct is None:
            st.progress(0.0, text="Starting…")
        else:
            st.progress(min(pct / 100.0, 1.0),
                        text=f"Thermal balance: {pct:.1f}% converged")
        if p.is_running():
            st.caption("⏳ Running… press **Stop 3D-PDR** to abort.")
        else:
            st.rerun()

    if running:
        _live_progress()
    elif bp is not None:
        bp.drain()
        rc = bp.returncode
        with st.expander("📋 View run log", expanded=(rc not in (0, None) and rc >= 0)):
            st.code("\n".join(bp.lines) or "(empty)", language="text")
        if rc == 0:
            st.success(f"✅ Run finished. Outputs in `{paths.pdr_root() / cur['outdir']}/{cur['output']}.*` "
                       "·  press **Next ➡️** to analyse them.")
        elif rc is not None and rc < 0:
            st.warning(f"⏹️ Run stopped by user (signal {-rc}).")
        else:
            st.error(f"❌ Run failed (exit {rc}). Check the log above.")

    # Auto-run RTtool right after the model finishes successfully so that the
    # Line emission page has results ready immediately (and won't trigger a
    # premature auto-run before model output exists).
    if finished_ok:
        _rt_run_prefix = str(cur["output"])
        if ss.get("rt_prefix") != _rt_run_prefix:
            rttool.write_paramsrt(str(cur["outdir"]), _rt_run_prefix,
                                  float(cur["vturb"]), float(cur["redshift"]))
            _rt_log: list[str] = []
            _rt_exit: int | None = None
            with st.status("Running RT-tool…", expanded=True) as _rt_status:
                _rt_ph = st.empty()
                for _line in runner.run_rttool():
                    _code = runner.parse_exit(_line)
                    if _code is not None:
                        _rt_exit = _code
                        break
                    _rt_log.append(_line)
                    _rt_ph.code("\n".join(_rt_log[-300:]), language="text")
                if _rt_exit == 0:
                    _rt_status.update(label="✅ RT-tool finished — press Next to view line emission",
                                      state="complete", expanded=False)
                    _rt_text = "\n".join(_rt_log)
                    rttool.save_tables_log(str(cur["outdir"]), _rt_run_prefix, _rt_text)
                    ss["rt_stdout"] = _rt_text
                    ss["rt_prefix"] = _rt_run_prefix
                else:
                    _rt_status.update(label=f"❌ RT-tool failed (exit {_rt_exit})", state="error")
                    ss["rt_failed_prefix"] = _rt_run_prefix

    _, bc = st.columns([4, 1])
    if finished_ok and bc.button("Next ➡️", key="run_next_bottom", type="primary", width="stretch"):
        _go_analysis()
