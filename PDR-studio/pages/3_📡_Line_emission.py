"""Feature 3 — Line emission.

Runs the RTtool post-processor on the finished 3D-PDR model and shows:
  • per-coolant line tables (Line, Freq, Tex, tau, Tr) at a depth point chosen
    via a Ntot/Av/N(H2) bar, limited per coolant;
  • how the antenna temperature Tr builds vs column density (one panel/coolant);
  • the same for optical depth;
  • the CO spectral-line energy distribution (SLED) at a chosen column density.

paramsRT.dat is filled automatically from params.dat — the user only presses Run.
"""
from __future__ import annotations

import contextlib
import io
import math
import traceback

import streamlit.components.v1 as _components

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from pdrstudio import config, outputs, params, paths, rttool, runner, ui

try:
    from streamlit_ace import st_ace
    _HAVE_ACE = True
except Exception:  # noqa: BLE001
    _HAVE_ACE = False

ui.configure_page("Line emission")
st.title("📡 Line emission")

# --- parameters (auto-filled from params.dat; read-only) -------------------
p = params.read_params(0)
outdir = str(p["outdir"])
prefix = str(p["output"])
vturb = float(p["vturb"])
redshift = float(p["redshift"])

st.caption("`paramsRT.dat` is filled automatically from `params.dat` — nothing to edit here.")
c1, c2, c3, c4 = st.columns(4)
c1.metric("Directory", outdir)
c2.metric("Prefix", prefix)
c3.metric("Linewidth", f"{vturb:g} km/s")
c4.metric("Redshift", f"{redshift:g}")
ui.fuv_cr_header(p, config.read_config())

# --- Auto-run RT-tool on first visit (or when prefix changes) ---------------
_rt_ready  = st.session_state.get("rt_prefix") == prefix
_rt_failed = st.session_state.get("rt_failed_prefix") == prefix

if _rt_failed:
    st.error("RT-tool failed on the last attempt — check that the model has finished running.")
    if st.button("🔄 Retry RT-tool"):
        del st.session_state["rt_failed_prefix"]
        st.rerun()

elif not _rt_ready:
    # Guard: only auto-run RTtool if model output actually exists.  Without
    # this check, navigating here before (or during) a model run would start
    # RTtool on stale/incomplete outputs.
    _model_output = paths.pdr_root() / outdir / f"{prefix}.pdr.fin"
    if not _model_output.exists():
        st.info(
            "No model output found for this prefix. "
            "Run the model in **📐 Model Parameters** first — "
            "RT-tool will start automatically when it finishes.",
            icon="ℹ️",
        )
    else:
        rttool.write_paramsrt(outdir, prefix, vturb, redshift)
        log: list[str] = []
        exit_code = None
        with st.status("Running RT-tool…", expanded=True) as status:
            ph = st.empty()
            if runner.rttool_needs_build():
                ph.write("Building RTtool (sources changed)…")
                for line in runner.compile_rttool():
                    if runner.parse_exit(line) is not None:
                        break
                    log.append(line); ph.code("\n".join(log[-300:]), language="text")
            for line in runner.run_rttool():
                code = runner.parse_exit(line)
                if code is not None:
                    exit_code = code
                    break
                log.append(line); ph.code("\n".join(log[-400:]), language="text")
            if exit_code == 0:
                status.update(label="✅ RT-tool finished", state="complete")
            else:
                status.update(label=f"❌ RT-tool failed (exit {exit_code})", state="error")
        if exit_code == 0:
            stdout_text = "\n".join(log)
            rttool.save_tables_log(outdir, prefix, stdout_text)
            st.session_state["rt_stdout"] = stdout_text
            st.session_state["rt_prefix"] = prefix
        else:
            st.session_state["rt_failed_prefix"] = prefix

tab_tab, tab_tr, tab_tau, tab_sled, tab_custom = st.tabs(
    ["📋 Line tables", "📈 Antenna temperature", "🌫️ Optical depth", "📶 CO SLED",
     "🖥️ Custom plot"]
)

# --- shared helpers --------------------------------------------------------
# CO, C, C⁺, O first, then the rest — used everywhere coolants are listed.
# Only show coolants the model was actually run with: fixed ones (CO, C+, C, O)
# plus whatever optional coolants are listed in params.dat.  Without this filter
# the page would also pick up Tr_*.dat files left over from a previous run that
# used a different coolant selection.
_expected_coolants: set[str] = {lbl for lbl, _ in params.FIXED_COOLANTS}
for _fn in p.get("coolants", []):
    _sp = params.coolant_species(_fn)
    if _sp:
        _expected_coolants.add(_sp)
coolants = rttool.order_coolants(
    [c for c in rttool.list_coolants(outdir, prefix) if c in _expected_coolants]
)


def _pos(a):
    a = np.asarray(a, dtype=float)
    return np.where(a > 0, a, np.nan)


def _load_model():
    """pdr Model for Av and the HI/H2 transition (cached per prefix)."""
    try:
        return outputs.load_pdr(outdir, prefix)
    except Exception:  # noqa: BLE001
        return None


def _x_arrays(tr: rttool.TrData):
    return {"Ntot": tr.Ntot, "N(H₂)": tr.NH2, "N(coolant)": tr.Nmol}


def _trans_x(m, x_arr):
    """Map the HI→H₂ transition (in Av) onto a column-density x-array."""
    if m is None:
        return None
    tav = outputs.transition_av(m.Av, m.abundance("H"), m.abundance("H2"))
    if tav is None:
        return None
    av = np.asarray(m.Av, dtype=float)[: len(x_arr)]
    return float(np.interp(tav, av, np.asarray(x_arr, dtype=float)[: len(av)]))


def _trans_x_2h2(m, x_arr):
    """Map the HI→2H₂ transition (50 % H₂ fraction) onto a column-density x-array."""
    if m is None:
        return None
    tav = outputs.transition_av_2h2(m.Av, m.abundance("H"), m.abundance("H2"))
    if tav is None:
        return None
    av = np.asarray(m.Av, dtype=float)[: len(x_arr)]
    return float(np.interp(tav, av, np.asarray(x_arr, dtype=float)[: len(av)]))


def _log_bounds(arr):
    """log10 bounds of the positive entries of arr, for a log-scale slider."""
    pos = np.asarray(arr, dtype=float)
    pos = pos[pos > 0]
    return float(np.log10(pos.min())), float(np.log10(pos.max()))


def _hi_h2_bar(m, arr, label, log=True):
    """Plotly bar marking the HI→H₂ transition on a column-density/Av axis.

    ``log=True`` draws a log10 axis (Ntot, N(H₂), N(coolant)); ``log=False``
    draws ``arr`` on a linear axis (Av). Returns ``(figure, transition value
    in the axis's own units)`` or ``(None, None)`` if unavailable.
    """
    lo, hi = _log_bounds(arr) if log else (float(np.min(arr)), float(np.max(arr)))
    trans_val = _trans_x(m, arr)
    if not trans_val or (log and trans_val <= 0) or hi <= lo:
        return None, None
    pos = float(np.log10(trans_val)) if log else float(trans_val)
    pos_txt = f"{pos:.2f}" if log else f"{pos:.3g}"
    bar = go.Figure()
    bar.add_trace(go.Scatter(x=[lo, hi], y=[0, 0], mode="lines",
                             line=dict(color="#bbb", width=3), hoverinfo="skip", showlegend=False))
    bar.add_annotation(x=min(max(pos, lo), hi), y=0, ax=0, ay=-26, text=f"HI/H₂ ({pos_txt})",
                       showarrow=True, arrowhead=2, arrowwidth=2, arrowcolor="crimson",
                       font=dict(color="crimson", size=12))
    bar.update_xaxes(range=[lo, hi], title_text=label, showgrid=False, zeroline=False)
    bar.update_yaxes(visible=False, range=[-1, 0.6])
    bar.update_layout(height=110, dragmode="pan", template="plotly_white",
                      paper_bgcolor="white", plot_bgcolor="white", font=dict(color="#262730"),
                      margin=dict(t=30, l=60, r=20, b=34))
    return bar, pos


# ===========================================================================
# Tab 1 — line tables (from the captured stdout)
# ===========================================================================
@st.fragment
def render_line_tables(stdout_text: str):
    """Line tables for every coolant at a single, user-chosen depth point.

    A bar lets the user pick that depth by Ntot, Av, or N(H₂) — log-scale for
    Ntot/N(H₂), linear for Av; each coolant's table (and its column density,
    shown next to the coolant name) is read straight from the per-depth
    Tr_/tau_/Tex_ files and updates live.
    """
    tables = rttool.parse_tables(stdout_text)
    trs = {cn: rttool.load_tr(outdir, prefix, cn) for cn in coolants}
    base = trs[coolants[0]]
    m = _load_model()

    x_arrays = {"Ntot": base.Ntot, "Av": base.Av, "N(H₂)": base.NH2}
    x_units = {"Ntot": "cm⁻²", "Av": "mag", "N(H₂)": "cm⁻²"}
    x_log = {"Ntot": True, "Av": False, "N(H₂)": True}
    _tr_attr = {"Ntot": "Ntot", "Av": "Av", "N(H₂)": "NH2"}

    with st.container():
        c1, c2 = st.columns([1, 2])
        with c1:
            x_kind = st.radio("Set depth by", list(x_arrays), horizontal=True, key="tab_x_kind")
        arr = x_arrays[x_kind]
        log = x_log[x_kind]
        with c2:
            if log:
                lo, hi = _log_bounds(arr)
                val = st.slider(f"log₁₀ {x_kind} [{x_units[x_kind]}]", min_value=lo, max_value=hi,
                                value=hi, step=(hi - lo) / 200 if hi > lo else 0.1,
                                key=f"tab_depth_{x_kind}")
                target = 10 ** val
            else:
                lo, hi = float(np.min(arr)), float(np.max(arr))
                target = st.slider(f"{x_kind} [{x_units[x_kind]}]", min_value=lo, max_value=hi,
                                   value=hi, step=(hi - lo) / 200 if hi > lo else 0.1,
                                   key=f"tab_depth_{x_kind}")
        row = rttool.nearest_row(arr, target)

        bar_label = f"log₁₀ {x_kind} [{x_units[x_kind]}]" if log else f"{x_kind} [{x_units[x_kind]}]"
        bar, pos = _hi_h2_bar(m, arr, bar_label, log=log)
        if bar is not None:
            pos_txt = f"log₁₀ {x_kind} = {pos:.2f}" if log else f"{x_kind} = {pos:.3g} mag"
            st.plotly_chart(bar, use_container_width=True, theme=None, key="tab1_hi_h2_bar")
            st.caption(f"⬆️ Arrow points the HI/H₂ transition — {pos_txt}.")

    st.caption(f"Depth point: Ntot={base.Ntot[row]:.3e} cm⁻², Av={base.Av[row]:.3g} mag, "
               f"N(H₂)={base.NH2[row]:.3e} cm⁻², Tgas={base.Tgas[row]:.4g} K")
    st.caption("Lowest ladder transitions shown. C⁺ shows 1 line (the [CII] "
               "158μm 1-0); C and O show 2; all others up to 10.")

    for cname in rttool.order_coolants(list(trs)):
        tr = trs[cname]
        tau = rttool.load_tau(outdir, prefix, cname)
        tex = rttool.load_tex(outdir, prefix, cname)
        freqs = rttool.transition_freqs(tables, cname, tr.ntrans)
        row_c = rttool.nearest_row(getattr(tr, _tr_attr[x_kind]), target)
        st.markdown(f"**{cname}** — N({cname}) = {tr.Nmol[row_c]:.3e} cm⁻²")
        # sort by ascending frequency (some coolants, e.g. OH hyperfine, are
        # stored out of frequency order in the LAMDA file); columns without a
        # matched frequency (shouldn't normally happen) sort last, by column.
        cols = sorted(range(tr.ntrans), key=lambda c: (freqs[c] is None, freqs[c] if freqs[c] is not None else c))
        st.dataframe(pd.DataFrame([{
            "Line": rttool.line_table_label(cname, c),
            "Freq (GHz)": f"{freqs[c]:.8g}" if freqs[c] is not None else "—",
            "Tex (K)": f"{tex[row_c, c]:.4g}",
            "tau": f"{tau[row_c, c]:.3e}",
            "Tr (K)": f"{tr.Tr[row_c, c]:.4g}",
        } for c in cols]), hide_index=True, use_container_width=True)


with tab_tab:
    # Collapse the 0-height component iframe so it doesn't leave a gap.
    st.markdown("""<style>
iframe[height="0"]{display:none!important}
</style>""", unsafe_allow_html=True)

    # Freeze the depth-controls block (Set depth by + slider + HI/H₂ bar).
    # We use a JS injection (same-origin iframe) to directly locate the element
    # and apply position:sticky, bypassing CSS selector guesswork.  A
    # MutationObserver reapplies the style after any Streamlit rerender.
    _components.html("""
<script>
(function () {
  var HEADER_H = 60; // px — Streamlit header height (3.75rem @ 16px base)

  function applySticky() {
    try {
      var doc = window.parent.document;

      // 1. Find the radio labelled "Set depth by" in the parent document.
      var target = null;
      doc.querySelectorAll('[data-testid="stRadio"]').forEach(function (r) {
        if (!target && r.innerText && r.innerText.includes('Set depth by'))
          target = r;
      });
      if (!target) return false;

      // 2. Walk up: stRadio → stElementContainer → stVerticalBlock(col)
      //    → stColumn → stHorizontalBlock → stVerticalBlock(our container)
      var el = target;
      for (var i = 0; i < 8; i++) {
        el = el.parentElement;
        if (!el) return false;
        if (el.getAttribute('data-testid') === 'stVerticalBlock') {
          var fc = el.firstElementChild;
          if (fc && fc.getAttribute('data-testid') === 'stHorizontalBlock') {
            // Found: el is our st.container() stVerticalBlock
            el.style.setProperty('position', 'sticky', 'important');
            el.style.setProperty('top', HEADER_H + 'px', 'important');
            el.style.setProperty('z-index', '100', 'important');
            el.style.setProperty('background', 'white', 'important');
            el.style.setProperty('padding-bottom', '8px', 'important');
            el.style.setProperty('box-shadow',
              '0 2px 6px rgba(0,0,0,0.08)', 'important');
            return true;
          }
        }
      }
      return false;
    } catch (e) { return false; }
  }

  // Keep trying until the element appears; after success keep a periodic
  // check so the style survives tab-switches and fragment rerenders.
  function watch() {
    applySticky();
    try {
      var obs = new MutationObserver(function () { applySticky(); });
      obs.observe(window.parent.document.body,
        { childList: true, subtree: true });
      // Stop observing after 2 min to avoid leaks; the tab would
      // inject a fresh iframe if the user navigates away and back.
      setTimeout(function () { obs.disconnect(); }, 120000);
    } catch (e) {}
  }

  watch();
  setTimeout(watch, 500);
})();
</script>
""", height=0, scrolling=False)

    stdout_text = None
    if st.session_state.get("rt_prefix") == prefix:
        stdout_text = st.session_state.get("rt_stdout")
    if stdout_text is None:
        stdout_text = rttool.load_tables_log(outdir, prefix)

    if not stdout_text or not coolants:
        st.info("Press **Run RT-tool** above to generate the line tables.", icon="▶️")
    else:
        render_line_tables(stdout_text)


# ===========================================================================
# Tabs 2 & 3 — Tr and tau build-up (one panel per coolant)
# ===========================================================================
@st.fragment
def render_profile(kind: str, value_label: str):
    """kind = 'Tr' or 'tau'. Draws controls + a dynamic panel-per-coolant figure.

    Wrapped in a fragment so the in-tab controls (x-axis, y-scale, per-coolant
    transition) rerun only this block — Streamlit otherwise jumps back to the
    first tab on the first widget interaction inside a non-default tab.
    """
    if not coolants:
        st.info("Press **Run RT-tool** above first (no RT-tool outputs found).", icon="▶️")
        return
    m = _load_model()

    c1, c2, c3 = st.columns([1.2, 1, 1])
    with c1:
        x_kind = st.radio("X-axis (log)", ["Ntot", "N(H₂)", "N(coolant)"],
                          horizontal=True, key=f"{kind}_x")
    with c2:
        yscale = st.radio("Y-axis", ["log", "linear"], horizontal=True, key=f"{kind}_y")
    with c3:
        show_trans_2h2 = st.checkbox("HI→2H₂ transition line", value=True, key=f"{kind}_t2",
                                     disabled=m is None,
                                     help="Dot-dashed line at x(H) = 2·x(H₂) — 50 % H₂ fraction")
        show_trans = st.checkbox("HI→H₂ transition line", value=False, key=f"{kind}_t",
                                 disabled=m is None,
                                 help="Dashed line at x(H) = x(H₂)")
    ylog = yscale == "log"

    # transition selector per coolant
    st.markdown("**Transition shown in each panel:**")
    sel: dict[str, int] = {}
    nctrl = min(5, len(coolants))
    ctrl_cols = st.columns(nctrl)
    for i, cn in enumerate(coolants):
        ntr = rttool.load_tr(outdir, prefix, cn).ntrans
        opts = [rttool.transition_short(cn, c) for c in range(ntr)]
        choice = ctrl_cols[i % nctrl].selectbox(cn, opts, key=f"{kind}_sel_{cn}")
        sel[cn] = opts.index(choice)

    # build panels
    ncols = min(3, len(coolants))
    nrows = math.ceil(len(coolants) / ncols)
    fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=coolants,
                        horizontal_spacing=0.08,
                        # extra row spacing so a panel's title doesn't sit on the
                        # x-axis label of the panel above it
                        vertical_spacing=max(0.14, 0.34 / nrows))
    xlabel = {"Ntot": "N<sub>tot</sub> [cm⁻²]", "N(H₂)": "N(H₂) [cm⁻²]",
              "N(coolant)": "N(coolant) [cm⁻²]"}[x_kind]

    for i, cn in enumerate(coolants):
        r, c = i // ncols + 1, i % ncols + 1
        tr = rttool.load_tr(outdir, prefix, cn)
        xarr = _x_arrays(tr)[x_kind]
        col = sel[cn]
        yarr = tr.Tr[:, col] if kind == "Tr" else rttool.load_tau(outdir, prefix, cn)[:, col]
        lab = f"{cn} {rttool.transition_short(cn, col)}"
        fig.add_trace(go.Scatter(x=np.asarray(xarr), y=_pos(yarr), mode="lines",
                                 name=lab, line=dict(width=2), showlegend=False,
                                 hovertemplate=f"{lab}<br>%{{x:.3e}}<br>%{{y:.3e}}<extra></extra>"),
                      row=r, col=c)
        fig.update_xaxes(type="log", title_text=xlabel, row=r, col=c,
                         exponentformat="e", showexponent="all",
                         gridcolor="rgba(0,0,0,0.08)")
        if ylog:
            ymax = np.nanmax(_pos(yarr))
            upper = math.log10(max(ymax, 1e-2)) + 0.3 if np.isfinite(ymax) else -2.0
            fig.update_yaxes(type="log", range=[-3, upper], title_text=value_label,
                             row=r, col=c, exponentformat="e", showexponent="all",
                             gridcolor="rgba(0,0,0,0.08)")
        else:
            fig.update_yaxes(type="linear", title_text=value_label, row=r, col=c,
                             gridcolor="rgba(0,0,0,0.08)")
        if show_trans:
            xt = _trans_x(m, xarr)
            if xt is not None:
                fig.add_vline(x=xt, line=dict(color="gray", width=1, dash="dash"),
                              row=r, col=c, exclude_empty_subplots=False)
        if show_trans_2h2:
            xt2 = _trans_x_2h2(m, xarr)
            if xt2 is not None:
                fig.add_vline(x=xt2, line=dict(color="gray", width=1, dash="dashdot"),
                              row=r, col=c, exclude_empty_subplots=False)

    fig.update_layout(height=max(360, 340 * nrows), dragmode="pan", template="plotly_white",
                      paper_bgcolor="white", plot_bgcolor="white", font=dict(color="#262730"),
                      title=dict(text=f"{prefix} — {value_label}", x=0.5, xanchor="center"),
                      margin=dict(t=60, l=60, r=20, b=50))
    # theme=None keeps white panels + dark text even in Streamlit dark mode
    st.plotly_chart(fig, use_container_width=True, theme=None, key=f"{kind}_profile_fig")


with tab_tr:
    render_profile("Tr", "Tr [K]")
with tab_tau:
    render_profile("tau", "optical depth τ")


# ===========================================================================
# Tab 4 — CO SLED
# ===========================================================================
@st.fragment
def render_sled():
    if "CO" not in coolants:
        st.info("No CO output found — run RT-tool on a model that includes CO.", icon="ℹ️")
        return
    co = rttool.load_tr(outdir, prefix, "CO")
    njmax = min(10, co.ntrans)
    st.markdown("CO spectral-line energy distribution — Tr of the first "
                f"{njmax} transitions (J→J−1) at a chosen column density. The right "
                "panel normalises the SLED to a transition of your choice.")

    dens_kind = st.radio("Set column density by", ["Ntot", "N(H₂)", "N(CO)"],
                         horizontal=True, key="sled_dens")
    arr = {"Ntot": co.Ntot, "N(H₂)": co.NH2, "N(CO)": co.Nmol}[dens_kind]
    lo, hi = _log_bounds(arr)
    logN = st.slider(f"log₁₀ {dens_kind} [cm⁻²]", min_value=lo, max_value=hi,
                     value=hi, step=(hi - lo) / 100 if hi > lo else 0.1, key="sled_logN")
    row = rttool.nearest_row(arr, 10 ** logN)

    bar, lt = _hi_h2_bar(_load_model(), arr, f"log₁₀ {dens_kind} [cm⁻²]")
    if bar is not None:
        st.plotly_chart(bar, use_container_width=True, theme=None, key="sled_hi_h2_bar")
        st.caption(f"⬆️ Arrow points the HI/H2 transition — log₁₀ {dens_kind} = {lt:.2f}.")

    st.caption(f"Nearest model depth: {dens_kind} = {arr[row]:.3e} cm⁻² "
               f"(Ntot={co.Ntot[row]:.3e}, N(H₂)={co.NH2[row]:.3e}, N(CO)={co.Nmol[row]:.3e})")

    J = np.arange(1, njmax + 1)
    tr_vals = np.asarray(co.Tr[row, :njmax], dtype=float)
    trans_labels = [f"{j}-{j-1}" for j in J]

    # controls: y-scale per panel + the normalisation transition
    cc1, cc2, cc3 = st.columns([1, 1.4, 1])
    with cc1:
        yscale_abs = st.radio("Left Y-axis", ["log", "linear"], horizontal=True,
                              key="sled_yabs")
    with cc2:
        ref_label = st.selectbox("Normalise by transition (J→J−1)", trans_labels,
                                 index=0, key="sled_ref",
                                 help="Every line in the right panel is divided by this one.")
    with cc3:
        yscale_norm = st.radio("Right Y-axis", ["log", "linear"], horizontal=True,
                               key="sled_ynorm")
    ref_idx = trans_labels.index(ref_label)
    ref_val = tr_vals[ref_idx]
    norm_vals = tr_vals / ref_val if ref_val not in (0.0, np.nan) and ref_val > 0 else np.full_like(tr_vals, np.nan)

    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.12,
                        subplot_titles=("CO SLED", f"Normalised to CO({ref_label})"))
    fig.add_trace(go.Scatter(x=J, y=_pos(tr_vals) if yscale_abs == "log" else tr_vals,
                             mode="lines+markers", line=dict(width=2, color="#d62728"),
                             showlegend=False,
                             hovertemplate="J=%{x}<br>Tr=%{y:.3e} K<extra></extra>"),
                  row=1, col=1)
    fig.add_trace(go.Scatter(x=J, y=_pos(norm_vals) if yscale_norm == "log" else norm_vals,
                             mode="lines+markers", line=dict(width=2, color="#1f77b4"),
                             showlegend=False,
                             hovertemplate="J=%{x}<br>Tr/Tr_ref=%{y:.3e}<extra></extra>"),
                  row=1, col=2)

    fig.update_xaxes(title_text="J (upper level; J→J−1)", dtick=1, tickmode="linear")
    fig.update_yaxes(title_text="Tr [K]", type="log" if yscale_abs == "log" else "linear",
                     exponentformat="e", showexponent="all", row=1, col=1)
    fig.update_yaxes(title_text=f"Tr / Tr[CO({ref_label})]",
                     type="log" if yscale_norm == "log" else "linear",
                     exponentformat="e", showexponent="all", row=1, col=2)
    fig.update_layout(height=480, dragmode="pan", template="plotly_white",
                      paper_bgcolor="white", plot_bgcolor="white", font=dict(color="#262730"),
                      title=dict(text=f"CO SLED — {prefix}", x=0.5, xanchor="center"),
                      margin=dict(t=70, l=60, r=20, b=50))
    st.plotly_chart(fig, use_container_width=True, theme=None, key="sled_fig")


with tab_sled:
    render_sled()


# ===========================================================================
# Tab 5 — custom plot
# ===========================================================================
RT_DEFAULT_CODE = '''\
# Custom plot — write any Python/matplotlib here. Preloaded for this RT-tool run:
#
#   coolants                     list of coolants with RT output
#   coolant_data("CO")           a record with .Ntot .NH2 .Nmol .Tgas .Tr .ntrans
#                                  .Tr[:, t] is the antenna temperature of the
#                                  (t+1)-(t) transition (t=0 -> 1-0) vs depth.
#   tau_data("CO")               optical-depth matrix [:, t], same indexing
#   np, plt, fig, ax             numpy, matplotlib, and a ready-made Figure/Axes
#
#   print(coolants)              # which coolants are available
#
# This starter plots the CO(1-0) antenna temperature vs N(H2):

co = coolant_data("CO")
ax.semilogx(co.NH2, co.Tr[:, 0], lw=2)

ax.set_xlabel("N(H$_2$) [cm$^{-2}$]")
ax.set_ylabel("T$_r$ [K]")
ax.set_title("CO(1-0) antenna temperature")
ax.grid(True, which="both", alpha=0.3)
'''


def _rt_run(code: str):
    fig_u, ax_u = plt.subplots(figsize=(6.6, 5.2))

    def coolant_data(name):
        return rttool.load_tr(outdir, prefix, name)

    def tau_data(name):
        return rttool.load_tau(outdir, prefix, name)

    ns = {"np": np, "plt": plt, "fig": fig_u, "ax": ax_u, "coolants": list(coolants),
          "coolant_data": coolant_data, "tau_data": tau_data, "prefix": prefix}
    out, err, png = io.StringIO(), None, None
    try:
        with contextlib.redirect_stdout(out):
            exec(code, ns)  # noqa: S102 — local single-user scratchpad
        rfig = ns.get("fig", fig_u)
        buf = io.BytesIO()
        rfig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
        png = buf.getvalue()
    except Exception:  # noqa: BLE001
        err = traceback.format_exc()
    finally:
        plt.close("all")
    st.session_state["rt_fig_png"] = png if png is not None else st.session_state.get("rt_fig_png")
    st.session_state["rt_stdout_c"] = out.getvalue()
    st.session_state["rt_error"] = err


@st.fragment
def _custom_plot_tab():
    if not coolants:
        st.info("Run RT-tool above first (no RT-tool outputs found).", icon="▶️")
        return
    st.caption("Write Python on the left, press **Run**, and your matplotlib figure "
               "appears on the right. The RT-tool antenna temperatures and optical "
               "depths are preloaded.")
    left, right = st.columns(2)
    with left:
        st.session_state.setdefault("rt_code", RT_DEFAULT_CODE)
        st.session_state.setdefault("rt_reset", 0)
        if _HAVE_ACE:
            code = st_ace(value=st.session_state["rt_code"], language="python", theme="monokai",
                          keybinding="vscode", font_size=14, tab_size=4, show_gutter=True,
                          wrap=False, auto_update=True, min_lines=40, max_lines=40,
                          key=f"rt_ace_{st.session_state['rt_reset']}")
            st.session_state["rt_code"] = code
        else:
            st.text_area("code", key="rt_code", height=720, label_visibility="collapsed")
        b1, b2 = st.columns(2)
        run_c = b1.button("▶ Run", type="primary", width="stretch", key="rt_run_c")
        reset_c = b2.button("↺ Reset code", width="stretch", key="rt_reset_c")

    if reset_c:
        st.session_state["rt_code"] = RT_DEFAULT_CODE
        st.session_state["rt_reset"] += 1
        st.rerun(scope="fragment")
    if run_c:
        _rt_run(st.session_state["rt_code"])

    with left:
        if st.session_state.get("rt_error"):
            st.markdown("**Console — error**")
            st.code(st.session_state["rt_error"], language="text")
        elif st.session_state.get("rt_stdout_c"):
            st.markdown("**Console**")
            st.code(st.session_state["rt_stdout_c"], language="text")
    with right:
        png = st.session_state.get("rt_fig_png")
        if png:
            st.image(png, width="stretch")
        else:
            st.info("Write code on the left and press **▶ Run** to render your figure here.",
                    icon="🖥️")


with tab_custom:
    _custom_plot_tab()
