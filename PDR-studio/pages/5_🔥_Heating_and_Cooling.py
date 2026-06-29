"""Feature 5 — Heating & Cooling.

Two tabs:
  • Profiles — a 2x2 log-log figure: heating, cooling, emissivities, level
    populations vs A_V (or n_H). The heating/cooling legends sit (stacked) in the
    top margin; the total heating/cooling is a thick translucent line.
  • Custom plot — a matplotlib scratchpad with the heating/cooling/emissivity/
    level-population data preloaded (default: total heating & cooling vs A_V).
"""
from __future__ import annotations

import contextlib
import io
import traceback

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from pdrstudio import config, outputs, params, paths, ui

try:
    from streamlit_ace import st_ace
    _HAVE_ACE = True
except Exception:  # noqa: BLE001
    _HAVE_ACE = False

ui.configure_page("Heating & Cooling")
st.title("🔥 Heating & Cooling")

# --- pick a model -----------------------------------------------------------
_p = params.read_params(0)
default_outdir = str(_p.get("outdir", "sims"))
c1, c2 = st.columns([1, 2])
with c1:
    outdir = st.text_input("Output directory", value=default_outdir,
                           help="Directory containing the model output files.")
models = outputs.list_models(outdir)
with c2:
    if not models:
        st.warning(f"No `*.pdr.fin` outputs found in `{paths.pdr_root() / outdir}`. "
                   "Run a model first (📐 Model Parameters → Run).", icon="⚠️")
        st.stop()
    latest = outputs.latest_model(outdir)
    default_idx = models.index(latest) if latest in models else 0
    prefix = st.selectbox("Model", models, index=default_idx)

try:
    m = outputs.load_pdr(outdir, prefix)
    h_labels, h_series, h_total = outputs.load_heating(outdir, prefix)
    c_labels, c_series, c_total = outputs.load_cooling(outdir, prefix)
except Exception as exc:  # noqa: BLE001
    st.error(f"Could not load model/heating/cooling for `{prefix}`: {exc}")
    st.stop()

x_HI, x_H2 = m.abundance("H"), m.abundance("H2")
trans = outputs.transition_av(m.Av, x_HI, x_H2)
trans_2h2 = outputs.transition_av_2h2(m.Av, x_HI, x_H2)
uniform = outputs.is_uniform_density(m.nH)
line_coolants = outputs.list_line_coolants(outdir, prefix)
spop_coolants = outputs.list_spop_coolants(outdir, prefix)

ui.fuv_cr_header(_p, config.read_config())

TOTAL_H, TOTAL_C = "Total heating", "Total cooling"
tab_profiles, tab_custom = st.tabs(["📈 Profiles", "🖥️ Custom plot"])


def _pos(a):
    a = np.asarray(a, dtype=float)
    return np.where(a > 0, a, np.nan)


def best_corner(Xarr, y_list, n_entries):
    with np.errstate(divide="ignore", invalid="ignore"):
        lx = np.log10(np.asarray(Xarr, dtype=float))
        xs, ys = [], []
        for y in y_list:
            ly = np.log10(np.asarray(y, dtype=float))
            good = np.isfinite(lx) & np.isfinite(ly)
            xs.append(lx[good]); ys.append(ly[good])
    px = np.concatenate(xs) if xs else np.array([])
    py = np.concatenate(ys) if ys else np.array([])
    if px.size == 0:
        return "top-left"
    xmin, xmax, ymin, ymax = px.min(), px.max(), py.min(), py.max()
    xr = (xmax - xmin) or 1.0
    yr = (ymax - ymin) or 1.0
    fx, fy = 0.42, min(0.65, 0.06 * n_entries + 0.12)
    boxes = {
        "top-left": (px <= xmin + fx * xr) & (py >= ymax - fy * yr),
        "top-right": (px >= xmax - fx * xr) & (py >= ymax - fy * yr),
        "bottom-left": (px <= xmin + fx * xr) & (py <= ymin + fy * yr),
        "bottom-right": (px >= xmax - fx * xr) & (py <= ymin + fy * yr),
    }
    counts = {k: int(v.sum()) for k, v in boxes.items()}
    order = ["top-left", "top-right", "bottom-left", "bottom-right"]
    return min(order, key=lambda k: (counts[k], order.index(k)))


# ===========================================================================
# Profiles tab
# ===========================================================================
with tab_profiles:
    cc1, cc2 = st.columns([1, 3])
    with cc1:
        if uniform:
            x_kind = "A_V"
            st.radio("X-axis", ["A_V"], index=0, disabled=True,
                     help="Uniform-density cloud — only A_V is available.")
        else:
            x_kind = st.radio("X-axis", ["A_V", "n_H"], index=0, horizontal=True)
    if x_kind == "n_H":
        X = np.asarray(m.nH, dtype=float); XLABEL = "n<sub>H</sub> [cm⁻³]"; x_hover = "n<sub>H</sub>"
        x_trans = outputs.value_at_crossing(m.nH, m.Av, x_HI, x_H2)
        x_trans_2h2 = outputs.value_at_crossing_2h2(m.nH, m.Av, x_HI, x_H2)
    else:
        X = np.asarray(m.Av, dtype=float); XLABEL = "A<sub>V</sub> [mag]"; x_hover = "A<sub>V</sub>"
        x_trans = trans
        x_trans_2h2 = trans_2h2
    with cc2:
        show_trans_2h2 = st.checkbox("Show HI→2H₂ transition", value=True,
                                     disabled=x_trans_2h2 is None,
                                     help="Dot-dashed line where x(H) = 2·x(H₂) — 50 % H₂ fraction")
        show_trans = st.checkbox("Show HI→H₂ transition", value=False,
                                 disabled=x_trans is None,
                                 help="Dashed line where x(H) = x(H₂)")

    heat_choices = h_labels + [TOTAL_H]
    cool_choices = c_labels + [TOTAL_C]
    r1c1, r1c2 = st.columns(2)
    heat_sel = r1c1.multiselect("Top-left — heating functions", heat_choices, default=heat_choices)
    cool_sel = r1c2.multiselect("Top-right — cooling functions", cool_choices, default=cool_choices)
    r2c1, r2c2 = st.columns(2)
    emis_options: dict[str, tuple[str, int]] = {}
    for cn in line_coolants:
        for t in range(1, outputs.line_ntrans(outdir, prefix, cn) + 1):
            emis_options[f"{cn}: {t}→{t-1}"] = (cn, t - 1)
    pop_options: dict[str, tuple[str, int]] = {}
    for cn in spop_coolants:
        for lv in range(1, outputs.spop_nlev(outdir, prefix, cn) + 1):
            pop_options[f"{cn}: level {lv-1}"] = (cn, lv - 1)
    emis_default = [k for k in ("CO: 1→0", "CO: 2→1") if k in emis_options][:1]
    emis_sel = r2c1.multiselect("Bottom-left — emissivities (coolant: upper→lower)",
                                list(emis_options.keys()), default=emis_default)
    pop_default = [k for k in ("CO: level 0", "CO: level 1") if k in pop_options][:2]
    pop_sel = r2c2.multiselect("Bottom-right — level populations (coolant: level)",
                               list(pop_options.keys()), default=pop_default)

    n = min(len(X), len(h_total), len(c_total))
    Xc = X[:n]
    TOTAL_LINE = dict(color="rgba(0,0,0,0.4)", width=4)
    _line_cache: dict[str, np.ndarray] = {}
    _spop_cache: dict[str, np.ndarray] = {}

    def _line(y, name, legend, color=None):
        return go.Scatter(x=Xc, y=_pos(np.asarray(y)[:n]), mode="lines", name=name,
                          line=dict(width=2, color=color), legend=legend,
                          hovertemplate=f"{name}<br>{x_hover}=%{{x:.3g}}<br>%{{y:.3e}}<extra></extra>")

    def _emis(cn, col):
        if cn not in _line_cache:
            _line_cache[cn] = outputs.load_line(outdir, prefix, cn)
        return _line_cache[cn][:, col]

    def _pop(cn, col):
        if cn not in _spop_cache:
            _spop_cache[cn] = outputs.load_spop(outdir, prefix, cn)
        return _spop_cache[cn][:, col]

    # Heating and cooling are placed on their OWN full-width rows so that a
    # horizontal legend above each panel never exceeds the panel width. The
    # narrower emissivity / level-population panels share the bottom row.
    fig = make_subplots(
        rows=3, cols=2,
        specs=[[{"colspan": 2}, None], [{"colspan": 2}, None], [{}, {}]],
        row_heights=[0.34, 0.34, 0.32], vertical_spacing=0.15, horizontal_spacing=0.12,
        subplot_titles=("Heating functions", "Cooling functions",
                        "Emissivities", "Level populations"))
    panel_y: dict[str, list] = {"legend": [], "legend2": [], "legend3": [], "legend4": []}

    for lab in heat_sel:
        if lab == TOTAL_H:
            fig.add_trace(go.Scatter(x=Xc, y=_pos(h_total[:n]), mode="lines", name=lab,
                                     line=TOTAL_LINE, legend="legend",
                                     hovertemplate=f"{lab}<br>{x_hover}=%{{x:.3g}}<br>%{{y:.3e}}<extra></extra>"),
                          row=1, col=1)
            panel_y["legend"].append(h_total[:n])
        elif lab in h_series:
            fig.add_trace(_line(h_series[lab], lab, "legend"), row=1, col=1)
            panel_y["legend"].append(h_series[lab][:n])
    for lab in cool_sel:
        if lab == TOTAL_C:
            fig.add_trace(go.Scatter(x=Xc, y=_pos(c_total[:n]), mode="lines", name=lab,
                                     line=TOTAL_LINE, legend="legend2",
                                     hovertemplate=f"{lab}<br>{x_hover}=%{{x:.3g}}<br>%{{y:.3e}}<extra></extra>"),
                          row=2, col=1)
            panel_y["legend2"].append(c_total[:n])
        elif lab in c_series:
            fig.add_trace(_line(c_series[lab], lab, "legend2"), row=2, col=1)
            panel_y["legend2"].append(c_series[lab][:n])
    for key in emis_sel:
        cn, col = emis_options[key]
        y = _emis(cn, col)
        fig.add_trace(_line(y, key, "legend3"), row=3, col=1)
        panel_y["legend3"].append(y[:n])
    for key in pop_sel:
        cn, col = pop_options[key]
        y = _pop(cn, col)
        fig.add_trace(_line(y, key, "legend4"), row=3, col=2)
        panel_y["legend4"].append(y[:n])

    ERG = "erg s⁻¹ cm⁻³"
    panel_axes = [(1, 1, ERG), (2, 1, ERG),
                  (3, 1, f"emissivity [{ERG}]"), (3, 2, "level population")]

    def _logtop(arrs):
        """Upper log10 bound (data max + 0.5 dex) for a list of arrays."""
        vals = np.concatenate([np.asarray(a, float).ravel() for a in arrs]) if arrs else np.array([])
        vals = vals[np.isfinite(vals) & (vals > 0)]
        return float(np.log10(vals.max()) + 0.5) if vals.size else -20.0

    # Heating/cooling (top two) are clamped to a 1e-35 floor so tiny terms don't
    # stretch the axis; the bottom panels keep their natural autorange.
    yrange = {(1, 1): [-35, _logtop(panel_y["legend"])],
              (2, 1): [-35, _logtop(panel_y["legend2"])]}
    for r, c, yt in panel_axes:
        fig.update_xaxes(type="log", title_text=XLABEL, row=r, col=c,
                         showgrid=True, gridcolor="rgba(0,0,0,0.08)",
                         exponentformat="e", showexponent="all")
        fig.update_yaxes(type="log", title_text=yt, row=r, col=c,
                         range=yrange.get((r, c)),
                         showgrid=True, gridcolor="rgba(0,0,0,0.08)",
                         exponentformat="e", showexponent="all")
    if show_trans and x_trans is not None:
        for r, c, _yt in panel_axes:
            fig.add_vline(x=x_trans, line=dict(color="gray", width=1, dash="dash"),
                          row=r, col=c, exclude_empty_subplots=False)
    if show_trans_2h2 and x_trans_2h2 is not None:
        for r, c, _yt in panel_axes:
            fig.add_vline(x=x_trans_2h2, line=dict(color="gray", width=1, dash="dashdot"),
                          row=r, col=c, exclude_empty_subplots=False)

    # --- legends ----------------------------------------------------------
    # Heating/cooling each get a horizontal legend ABOVE their (full-width) panel;
    # the legend wraps into ~2 rows and never exceeds the panel width. The bottom
    # half-width panels use a compact inset legend in the emptiest corner.
    _INSET = 0.008

    def _entrywidth(k):
        per = max(1, -(-k // 2))   # ceil(k/2) -> about two rows
        return min(0.32, 1.0 / per)

    def _top_legend(title, k, ytop):
        return dict(orientation="h", entrywidthmode="fraction", entrywidth=_entrywidth(k),
                    xref="paper", yref="paper", x=0.5, xanchor="center",
                    y=ytop + 0.03, yanchor="bottom",
                    title=dict(text=title, side="left"), font=dict(size=12),
                    bgcolor="rgba(245,245,245,0.92)", bordercolor="rgba(0,0,0,0.2)", borderwidth=1)

    def _inset_legend(xaxis, yaxis, corner):
        xd, yd = xaxis.domain, yaxis.domain
        left, top = "left" in corner, "top" in corner
        return dict(xref="paper", yref="paper",
                    x=(xd[0] + _INSET) if left else (xd[1] - _INSET),
                    y=(yd[1] - _INSET) if top else (yd[0] + _INSET),
                    xanchor="left" if left else "right", yanchor="top" if top else "bottom",
                    bgcolor="rgba(255,255,255,0.65)", bordercolor="rgba(0,0,0,0.25)",
                    borderwidth=1, font=dict(size=9))

    legends: dict = {}
    if panel_y["legend"]:
        legends["legend"] = _top_legend("Heating", len(panel_y["legend"]),
                                        fig.layout.yaxis.domain[1])
    if panel_y["legend2"]:
        legends["legend2"] = _top_legend("Cooling", len(panel_y["legend2"]),
                                         fig.layout.yaxis2.domain[1])
    for lid, xax, yax in [("legend3", fig.layout.xaxis3, fig.layout.yaxis3),
                          ("legend4", fig.layout.xaxis4, fig.layout.yaxis4)]:
        if panel_y[lid]:
            legends[lid] = _inset_legend(xax, yax, best_corner(Xc, panel_y[lid], len(panel_y[lid])))

    fig.update_layout(height=1240, dragmode="pan", margin=dict(t=120, l=70, r=30, b=55),
                      template="plotly_white", paper_bgcolor="white", plot_bgcolor="white",
                      font=dict(color="#262730"), **legends)

    st.markdown(f"**{prefix}** · network: {m.network}")
    # theme=None keeps white panels + dark legend text even in dark mode
    st.plotly_chart(fig, use_container_width=True, theme=None)
    _tcap = ""
    if trans: _tcap += f"  Dashed: HI→H₂ (x(H)=x(H₂)), Aᵥ={trans:.4g} mag."
    if trans_2h2: _tcap += f"  Dot-dashed: HI→2H₂ (x(H)=2·x(H₂), 50 % H₂), Aᵥ={trans_2h2:.4g} mag."
    st.caption("Top row: heating functions. Middle row: cooling functions. Thick gray-transparent line: total heating / cooling. "
               "Bottom row: per-coolant `.line.fin` emissivities and `.spop.fin` level populations. "
               "Drag to pan, scroll to zoom." + _tcap)

# ===========================================================================
# Custom plot tab
# ===========================================================================
DEFAULT_CODE = '''\
# Custom plot — write any Python/matplotlib here. Preloaded for this model:
#
#   Av, nH, Tgas, Tdust, x        1-D arrays over the depth grid (edge -> interior)
#   total_heating, total_cooling  total rates [erg s^-1 cm^-3]
#   heating("PAH photoelectric")  a specific heating term  (see heating_labels)
#   cooling("CO")                 a specific coolant's cooling (see cooling_labels)
#   emissivity("CO", 0)           line emissivity, transition index 0 = (1-0)
#   level_pop("CO", 0)            level population, index 0 = level 1
#   line_coolants, spop_coolants  coolants available for emissivity / level_pop
#   np, plt, fig, ax              numpy, matplotlib, and a ready-made Figure/Axes
#
#   print(heating_labels)   # list the available heating terms
#   print(cooling_labels)   # list the available coolants
#
# This starter plots the total heating and total cooling vs Av:

ax.loglog(Av, total_heating, label="Total heating")
ax.loglog(Av, total_cooling, label="Total cooling")

ax.set_xlabel("A$_V$ [mag]")
ax.set_ylabel("erg s$^{-1}$ cm$^{-3}$")
ax.set_title("Total heating vs cooling")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
'''


def _run_user_code(code: str):
    fig_u, ax_u = plt.subplots(figsize=(6.6, 5.2))

    def heating(label):
        return h_series.get(label)

    def cooling(label):
        return c_series.get(label)

    def emissivity(coolant, t=0):
        return outputs.load_line(outdir, prefix, coolant)[:, t]

    def level_pop(coolant, lv=0):
        return outputs.load_spop(outdir, prefix, coolant)[:, lv]

    ns = {
        "np": np, "plt": plt, "fig": fig_u, "ax": ax_u, "m": m,
        "Av": m.Av, "nH": m.nH, "Tgas": m.Tgas, "Tdust": m.Tdust, "x": m.x,
        "total_heating": h_total, "total_cooling": c_total,
        "heating": heating, "cooling": cooling,
        "emissivity": emissivity, "level_pop": level_pop,
        "heating_labels": list(h_labels), "cooling_labels": list(c_labels),
        "line_coolants": list(line_coolants), "spop_coolants": list(spop_coolants),
        "prefix": prefix, "network": m.network,
    }
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
    st.session_state["a2_fig_png"] = png if png is not None else st.session_state.get("a2_fig_png")
    st.session_state["a2_stdout"] = out.getvalue()
    st.session_state["a2_error"] = err


with tab_custom:
    st.caption("Write Python on the left, press **Run**, and your matplotlib figure appears "
               "on the right. Heating/cooling/emissivity/level-population data are preloaded.")
    left, right = st.columns(2)
    with left:
        st.session_state.setdefault("a2_code", DEFAULT_CODE)
        st.session_state.setdefault("a2_reset", 0)
        if _HAVE_ACE:
            # fixed ~40-line height with a scrollbar; does not auto-grow as you type
            code = st_ace(value=st.session_state["a2_code"], language="python", theme="monokai",
                          keybinding="vscode", font_size=14, tab_size=4, show_gutter=True,
                          wrap=False, auto_update=True, height=720,
                          key=f"a2_ace_{st.session_state['a2_reset']}")
            st.session_state["a2_code"] = code
        else:
            st.text_area("code", key="a2_code", height=720, label_visibility="collapsed")
        b1, b2 = st.columns(2)
        run = b1.button("▶ Run", type="primary", width="stretch", key="a2_run")
        reset = b2.button("↺ Reset code", width="stretch", key="a2_reset_btn")

    if reset:
        st.session_state["a2_code"] = DEFAULT_CODE
        st.session_state["a2_reset"] += 1
        st.rerun()
    if run:
        _run_user_code(st.session_state["a2_code"])

    with left:
        if st.session_state.get("a2_error"):
            st.markdown("**Console — error**")
            st.code(st.session_state["a2_error"], language="text")
        elif st.session_state.get("a2_stdout"):
            st.markdown("**Console**")
            st.code(st.session_state["a2_stdout"], language="text")
    with right:
        png = st.session_state.get("a2_fig_png")
        if png:
            st.image(png, width="stretch")
        else:
            st.info("Write code on the left and press **▶ Run** to render your figure here.",
                    icon="🖥️")
