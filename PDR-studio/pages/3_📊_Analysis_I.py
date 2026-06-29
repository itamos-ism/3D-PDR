"""Feature 4 — Analysis I.

Two tabs:
  • Profiles — an interactive 2x2 log-log Plotly figure (Tgas, HI/H2, C+/C/CO,
    user-selected species) versus A_V (or n_H), with the HI->H2 transition line.
  • Custom plot — a Python "scratchpad": write matplotlib code on the left, see
    the figure on the right. The model's arrays (Av, nH, Tgas, …) and an
    abundance() helper are pre-loaded.
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

from pdrstudio import outputs, params, paths, ui

try:  # syntax-highlighting code editor (optional; falls back to a plain text area)
    from streamlit_ace import st_ace
    _HAVE_ACE = True
except Exception:  # noqa: BLE001
    _HAVE_ACE = False

ui.configure_page("Analysis I")
st.title("📊 Analysis I — temperature & abundances")

# --- pick a model ----------------------------------------------------------
default_outdir = str(params.read_params(0).get("outdir", "sims"))
c1, c2 = st.columns([1, 2])
with c1:
    outdir = st.text_input("Output directory", value=default_outdir,
                           help="Directory containing the model `.pdr.fin` files.")
models = outputs.list_models(outdir)
with c2:
    if not models:
        st.warning(f"No `*.pdr.fin` outputs found in `{paths.pdr_root() / outdir}`. "
                   "Run a model first (📐 Model Parameters → Run).", icon="⚠️")
        st.stop()
    latest = outputs.latest_model(outdir)
    default_idx = models.index(latest) if latest in models else 0
    prefix = st.selectbox("Model", models, index=default_idx,
                          help="Output prefix to analyse. Defaults to the most recently run model.")

try:
    m = outputs.load_pdr(outdir, prefix)
except Exception as exc:  # noqa: BLE001 - surface any load error to the user
    st.error(f"Could not load `{prefix}.pdr.fin`: {exc}")
    st.stop()

x_HI, x_H2 = m.abundance("H"), m.abundance("H2")
trans = outputs.transition_av(m.Av, x_HI, x_H2)
uniform = outputs.is_uniform_density(m.nH)

cols = st.columns(4)
cols[0].metric("Network", m.network or "unknown")
cols[1].metric("Grid points", len(m.Av))
cols[2].metric("Tgas range", f"{m.Tgas.min():.0f}–{m.Tgas.max():.0f} K")
cols[3].metric("HI→H₂ transition", f"Aᵥ = {trans:.3g}" if trans else "—")


# --- shared plotting helpers (pure) ----------------------------------------
def _pos(a):
    """Mask non-positive values (None) so they leave gaps on log axes."""
    a = np.asarray(a, dtype=float)
    return np.where(a > 0, a, np.nan)


def best_corner(Xarr, y_list, n_entries):
    """matplotlib-style ``loc='best'``: pick the panel corner (in log space) with
    the least data in it, so an inset legend hides as few curves as possible."""
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


_INSET = 0.008


def _panel_legend(xaxis, yaxis, corner):
    xd, yd = xaxis.domain, yaxis.domain
    left = "left" in corner
    top = "top" in corner
    return dict(
        xref="paper", yref="paper",
        x=(xd[0] + _INSET) if left else (xd[1] - _INSET),
        y=(yd[1] - _INSET) if top else (yd[0] + _INSET),
        xanchor="left" if left else "right",
        yanchor="top" if top else "bottom",
        bgcolor="rgba(255,255,255,0.65)", bordercolor="rgba(0,0,0,0.25)", borderwidth=1,
        font=dict(size=9),
    )


tab_profiles, tab_custom = st.tabs(["📈 Profiles", "🖥️ Custom plot"])

# ===========================================================================
# Tab 1 — the standard 2x2 profiles figure
# ===========================================================================
with tab_profiles:
    c1, c2, c3 = st.columns([2, 1, 1])
    with c1:
        user_species = st.multiselect(
            "Bottom-right panel — species to plot (in species-file order)",
            options=m.species,
            default=[s for s in ("HCO+", "H2O", "O2") if s in m.species][:1] or None,
            help="Choose any species from the network; order follows species_<network>.d.",
        )
    with c2:
        x_opts = ["A_V", "N_tot"] + ([] if uniform else ["n_H"])
        x_kind = st.radio(
            "X-axis", x_opts, index=0,
            help="A_V, the line-of-sight total column density N_tot = ∫ n_H dl, "
                 + ("or the local density n_H." if not uniform
                    else "(uniform-density cloud — n_H is constant)."))
        if uniform:
            st.caption(f"Uniform density: n_H = {m.nH.mean():.3g} cm⁻³")
    with c3:
        show_trans = st.checkbox("Show HI–H₂ transition line on all panels", value=True,
                                 disabled=trans is None)

    if x_kind == "n_H":
        X = np.asarray(m.nH, dtype=float)
        XLABEL = "n<sub>H</sub> [cm⁻³]"
        x_trans = outputs.value_at_crossing(m.nH, m.Av, x_HI, x_H2)
        x_hover = "n<sub>H</sub>"
    elif x_kind == "N_tot":
        X = outputs.column_density(m.x, m.nH)
        XLABEL = "N<sub>tot</sub> [cm⁻²]"
        x_trans = outputs.value_at_crossing(X, m.Av, x_HI, x_H2)
        x_hover = "N<sub>tot</sub>"
    else:
        X = np.asarray(m.Av, dtype=float)
        XLABEL = "A<sub>V</sub> [mag]"
        x_trans = trans
        x_hover = "A<sub>V</sub>"

    ABLABEL = "n(X) / n<sub>H</sub>"

    def _line(y, name, legend, color=None):
        return go.Scatter(
            x=X, y=_pos(y), mode="lines", name=name, legend=legend,
            line=dict(width=2, color=color),
            hovertemplate=f"{name}<br>{x_hover}=%{{x:.3g}}<br>%{{y:.3e}}<extra></extra>",
        )

    fig = make_subplots(
        rows=2, cols=2, horizontal_spacing=0.09, vertical_spacing=0.13,
        subplot_titles=("Gas temperature", "HI / H₂ transition",
                        "C⁺ / C / CO", "Selected species"),
    )
    panel_y: dict[str, list] = {"legend": [], "legend2": [], "legend3": [], "legend4": []}

    fig.add_trace(_line(m.Tgas, "T_gas", "legend", "#d62728"), row=1, col=1)
    panel_y["legend"].append(m.Tgas)
    if x_HI is not None:
        fig.add_trace(_line(x_HI, "HI (H)", "legend2", "#1f77b4"), row=1, col=2)
        panel_y["legend2"].append(x_HI)
    if x_H2 is not None:
        fig.add_trace(_line(x_H2, "H₂", "legend2", "#ff7f0e"), row=1, col=2)
        panel_y["legend2"].append(x_H2)
    for tok, color, lab in [("C+", "#1f77b4", "C⁺"), ("C", "#2ca02c", "C"), ("CO", "#9467bd", "CO")]:
        a = m.abundance(tok)
        if a is not None:
            fig.add_trace(_line(a, lab, "legend3", color), row=2, col=1)
            panel_y["legend3"].append(a)
    for tok in user_species:
        a = m.abundance(tok)
        if a is not None:
            fig.add_trace(_line(a, tok, "legend4"), row=2, col=2)
            panel_y["legend4"].append(a)

    for r in (1, 2):
        for c in (1, 2):
            fig.update_xaxes(type="log", title_text=XLABEL, row=r, col=c,
                             showgrid=True, gridcolor="rgba(0,0,0,0.08)",
                             exponentformat="e", showexponent="all")
            ytitle = "T_gas [K]" if (r, c) == (1, 1) else ABLABEL
            fig.update_yaxes(type="log", title_text=ytitle, row=r, col=c,
                             showgrid=True, gridcolor="rgba(0,0,0,0.08)",
                             exponentformat="e", showexponent="all")

    if show_trans and x_trans is not None:
        for r in (1, 2):
            for c in (1, 2):
                fig.add_vline(x=x_trans, line=dict(color="gray", width=1, dash="dash"),
                              row=r, col=c, exclude_empty_subplots=False)

    _axpair = {
        "legend": (fig.layout.xaxis, fig.layout.yaxis),
        "legend2": (fig.layout.xaxis2, fig.layout.yaxis2),
        "legend3": (fig.layout.xaxis3, fig.layout.yaxis3),
        "legend4": (fig.layout.xaxis4, fig.layout.yaxis4),
    }
    legends = {}
    for lid, (xax, yax) in _axpair.items():
        corner = best_corner(X, panel_y[lid], len(panel_y[lid]))
        legends[lid] = _panel_legend(xax, yax, corner)

    fig.update_layout(
        height=720,
        title=dict(text=f"{prefix}   (network: {m.network})", x=0.5, xanchor="center"),
        margin=dict(t=70, l=70, r=30, b=70),
        template="plotly_white",
        paper_bgcolor="white", plot_bgcolor="white", font=dict(color="#262730"),
        **legends,
    )
    # theme=None keeps the white plotly_white template (white panels, dark text)
    # even when Streamlit is in dark mode.
    st.plotly_chart(fig, use_container_width=True, theme=None)

    st.caption(
        f"All panels are log–log versus {x_kind}. Drag to zoom, double-click to reset, "
        "hover for values, and click legend entries to toggle traces. The thin dashed "
        "line marks the HI→H₂ transition (where the HI and H₂ abundance curves cross)."
        + (f"  Aᵥ = {trans:.4g} mag." if trans else "")
        + ("" if not uniform else "  (Uniform-density cloud — x-axis fixed to A_V.)")
    )

# ===========================================================================
# Tab 2 — custom Python plot (editor left, figure right)
# ===========================================================================
DEFAULT_CODE = '''\
# Custom plot — write any Python/matplotlib here. The current model is preloaded:
#
#   Av, nH, Tgas, Tdust, x, UV   1-D arrays over the depth grid (edge -> interior)
#                                 x is position [pc]; UV the FUV field; nH density.
#   species                      list of every species token in the network
#   abundance("CO")              n(X)/n_H array for a species (None if not present)
#   np, plt                      numpy and matplotlib.pyplot
#   fig, ax                      a ready-made Figure / Axes to draw on
#   network, prefix              network name and model prefix (strings)
#
# Tips:
#   - plot the gas temperature:      ax.plot(Av, Tgas)
#   - any species (log-log):         ax.loglog(Av, abundance("H2O"), label="H2O")
#   - several at once:               for s in ["CO","C","C+"]: ax.loglog(Av, abundance(s), label=s)
#   - second axis for temperature:   ax2 = ax.twinx(); ax2.plot(Av, Tgas, "k--")
#   - print to the console (left):   print(species)
#
# This starter reproduces the HI -> H2 transition (atomic vs molecular hydrogen):

ax.loglog(Av, abundance("H"),  label="HI (H)")
ax.loglog(Av, abundance("H2"), label="H$_2$")

ax.set_xlabel("A$_V$ [mag]")
ax.set_ylabel("n(X) / n$_\\\\mathrm{H}$")
ax.set_title("HI $\\\\to$ H$_2$ transition")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
'''


def _run_user_code(code: str):
    """Execute the user's code with the model data preloaded; capture the figure,
    stdout and any traceback into session state."""
    fig_u, ax_u = plt.subplots(figsize=(6.6, 5.2))
    ns = {
        "np": np, "plt": plt, "fig": fig_u, "ax": ax_u,
        "m": m, "model": m, "abundance": m.abundance, "species": list(m.species),
        "Av": m.Av, "nH": m.nH, "Tgas": m.Tgas, "Tdust": m.Tdust, "x": m.x, "UV": m.UV,
        "network": m.network, "prefix": prefix,
    }
    out, err, png = io.StringIO(), None, None
    try:
        with contextlib.redirect_stdout(out):
            exec(code, ns)  # noqa: S102 — local single-user scratchpad
        rfig = ns.get("fig", fig_u)
        buf = io.BytesIO()
        rfig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
        png = buf.getvalue()
    except Exception:  # noqa: BLE001 — show the traceback in the console
        err = traceback.format_exc()
    finally:
        plt.close("all")
    st.session_state["ai_fig_png"] = png if png is not None else st.session_state.get("ai_fig_png")
    st.session_state["ai_stdout"] = out.getvalue()
    st.session_state["ai_error"] = err


with tab_custom:
    st.caption("Write Python on the left, press **Run**, and your matplotlib figure "
               "appears on the right. The model's arrays and `abundance()` are preloaded.")
    left, right = st.columns(2)

    with left:
        st.session_state.setdefault("ai_code", DEFAULT_CODE)
        st.session_state.setdefault("ai_reset", 0)
        if _HAVE_ACE:
            # coloured editor with a FIXED height (~40 lines): a vertical scrollbar
            # appears past that, and it does NOT auto-grow as you type.
            # auto_update=True so Run always uses the live content (no Ctrl-Enter).
            code = st_ace(
                value=st.session_state["ai_code"], language="python", theme="monokai",
                keybinding="vscode", font_size=14, tab_size=4, show_gutter=True,
                wrap=False, auto_update=True, height=720,
                key=f"ace_{st.session_state['ai_reset']}",
            )
            st.session_state["ai_code"] = code
        else:
            st.text_area("code", key="ai_code", height=720, label_visibility="collapsed")
        b1, b2 = st.columns(2)
        run = b1.button("▶ Run", type="primary", width="stretch")
        reset = b2.button("↺ Reset code", width="stretch")

    if reset:
        st.session_state["ai_code"] = DEFAULT_CODE
        st.session_state["ai_reset"] += 1
        st.rerun()
    if run:
        _run_user_code(st.session_state["ai_code"])

    with left:
        if st.session_state.get("ai_error"):
            st.markdown("**Console — error**")
            st.code(st.session_state["ai_error"], language="text")
        elif st.session_state.get("ai_stdout"):
            st.markdown("**Console**")
            st.code(st.session_state["ai_stdout"], language="text")

    with right:
        png = st.session_state.get("ai_fig_png")
        if png:
            st.image(png, width="stretch")
        else:
            st.info("Write code on the left and press **▶ Run** to render your figure here.",
                    icon="🖥️")
