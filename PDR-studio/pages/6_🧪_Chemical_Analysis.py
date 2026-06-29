"""Feature 7 — Chemical Analysis.

Reads ``<prefix>.chemanalysis.fin`` (written by src/analyse_chem.F90) and lets
the user track ONE species: its abundance profile versus A_V (or n_H for a
non-uniform cloud), in log–log space, with the optional HI→H₂ transition line.

Below the plot a panel reproduces the per-point pathway report from the
original analyse_chem routine (formation/destruction reactions and their
percentage contributions, plus the totals). As the cursor moves over the
abundance curve, that panel updates automatically to the hovered depth point —
this is done client-side (a Plotly ``plotly_hover`` handler rewriting the
table) so it is instantaneous and needs no server round-trip.
"""
from __future__ import annotations

import functools
import html
import json

import numpy as np
import streamlit as st
from plotly.offline import get_plotlyjs

from pdrstudio import chemanalysis as ca
from pdrstudio import config, llm, outputs, params, paths, ui

ui.configure_page("Chemical Analysis")
st.title("🧪 Chemical Analysis — formation & destruction pathways")


@functools.lru_cache(maxsize=1)
def _plotlyjs() -> str:
    """The offline Plotly.js bundle (cached so we build the string once)."""
    return get_plotlyjs()


# --- pick a model -----------------------------------------------------------
_p = params.read_params(0)
default_outdir = str(_p.get("outdir", "sims"))
c1, c2 = st.columns([1, 2])
with c1:
    outdir = st.text_input("Output directory", value=default_outdir,
                           help="Directory containing the `*.chemanalysis.fin` files.")
models = ca.list_models(outdir)
with c2:
    if not models:
        st.warning(
            f"No `*.chemanalysis.fin` files found in `{paths.pdr_root() / outdir}`. "
            "Build with the `CHEMANALYSIS` flag enabled (⚙️ Code Configuration) and "
            "run a model (📐 Model Parameters → Run).", icon="⚠️")
        st.stop()
    latest = ca.latest_model(outdir)
    default_idx = models.index(latest) if latest in models else 0
    prefix = st.selectbox("Model", models, index=default_idx,
                          help="Defaults to the most recently produced chemical-analysis file.")

try:
    m = ca.load(outdir, prefix)
except Exception as exc:  # noqa: BLE001 — surface any load error to the user
    st.error(f"Could not read `{prefix}.chemanalysis.fin`: {exc}")
    st.stop()

if not m.points or not m.species_list:
    st.error(f"`{prefix}.chemanalysis.fin` contains no usable data.")
    st.stop()

network = outputs.detect_network(m.nspec)
uniform = m.is_uniform_density()

cols = st.columns(4)
cols[0].metric("Network", network or f"{m.nspec} species")
cols[1].metric("Grid points", len(m.points))
cols[2].metric("Reactions", m.nreac)
cols[3].metric("Chemistry time", f"{m.time:.3g} yr" if m.time else "—")
ui.fuv_cr_header(_p, config.read_config())

# --- controls ---------------------------------------------------------------
c1, c2, c3 = st.columns([2, 1, 1])
with c1:
    default_sp = next((s for s in ("CO", "C", "C+", "H2") if s in m.species_list),
                      m.species_list[0])
    species = st.selectbox(
        "Species to track", options=m.species_list,
        index=m.species_list.index(default_sp),
        help="Only one species at a time. Its abundance profile is plotted; the "
             "panel below shows that species' pathways at the hovered depth point.")
with c2:
    if uniform:
        x_kind = "Av"
        st.radio("X-axis", ["A_V"], index=0, disabled=True,
                 help="Uniform-density cloud — n_H is constant, so only A_V is available.")
        st.caption(f"Uniform density: n_H = {m.x_axis('nH').mean():.3g} cm⁻³")
    else:
        choice = st.radio("X-axis", ["A_V", "n_H"], index=0,
                          help="This cloud has a varying density profile, so you can "
                               "plot against n_H as well as A_V.")
        x_kind = "nH" if choice == "n_H" else "Av"
with c3:
    x_HI, x_H2 = m.abundance("H"), m.abundance("H2")
    trans_av = outputs.transition_av(m.x_axis("Av"), x_HI, x_H2)
    show_trans = st.checkbox("Show HI–H₂ transition line", value=True,
                             disabled=trans_av is None)

# --- AI pathway summary (optional, fully-local LLM) ------------------------
backend = llm.backend_name()
sc1, sc2 = st.columns([1.4, 3])
with sc1:
    summarize_clicked = st.button(
        f"🧠 Summarize {species} pathways", type="primary",
        help="Generate a short natural-language summary of how this species forms "
             "and is destroyed along the density/depth column, using a locally-run "
             "LLM (offline). The summary is based on the deterministic digest below.")
with sc2:
    st.caption(f"Local LLM backend: **{backend}**" if backend else
               "Local LLM backend: _not detected_ — see setup options below.")

# When no local LLM is present, show the install options up-front (not only on click).
if backend is None:
    with st.expander("🧠 Enable offline AI summaries (optional) — Llama not detected"):
        st.markdown(llm.setup_hint())

_sum_key = f"chemsum::{prefix}::{species}"
if summarize_clicked:
    if backend is None:
        st.info(llm.setup_hint(), icon="🧠")
    else:
        evidence = ca.pathway_evidence(m, species)
        if not evidence:
            st.warning(f"No pathway data available for {species} in this model.")
        else:
            try:
                with st.spinner(f"Summarizing {species} pathways with {backend}…"):
                    st.session_state[_sum_key] = llm.summarize(evidence)
            except Exception as exc:  # noqa: BLE001 — surface backend errors
                st.error(f"LLM summary failed: {exc}")

if st.session_state.get(_sum_key):
    st.markdown(f"#### 🧠 {species} — pathway summary")
    st.markdown(st.session_state[_sum_key])
    st.caption("Generated locally from the model's chemical-analysis output. "
               "The deterministic pathway table below remains the source of truth.")
    with st.expander("What the summary is based on (deterministic digest)"):
        st.code(llm.build_prompt(ca.pathway_evidence(m, species)), language="text")
st.divider()

# --- assemble plot data (sorted by the chosen x so the line reads left→right) -
X = m.x_axis(x_kind)
Y = m.abundance(species)
order = np.argsort(X)
pts_sorted = [m.points[i] for i in order]
Xs, Ys = X[order], Y[order]

XLABEL = "n<sub>H</sub> [cm⁻³]" if x_kind == "nH" else "A<sub>V</sub> [mag]"
x_short = "n_H" if x_kind == "nH" else "A_V"
YLABEL = f"n({species}) / n<sub>H</sub>"

# transition position on the chosen x-axis
x_trans = (outputs.value_at_crossing(m.x_axis(x_kind), m.x_axis("Av"), x_HI, x_H2)
           if show_trans else None)


def _logsafe(a: np.ndarray) -> list[float | None]:
    """Non-positive values -> None (gaps on a log axis), and JSON-serialisable."""
    return [float(v) if (np.isfinite(v) and v > 0) else None for v in a]


# --- per-point pathway report (the analyse_chem table), one HTML block each --
def _react_rows(pairs: list[tuple[int, int]], sign: str) -> str:
    if not pairs:
        return '<tr><td colspan="2" class="none">— none —</td></tr>'
    rows = []
    for rid, pct in pairs:
        eq = html.escape(m.reactions.get(rid, f"reaction {rid}"))
        rows.append(f'<tr><td class="eq">{eq}</td>'
                    f'<td class="pct {sign}">{pct}%</td></tr>')
    return "".join(rows)


def _point_html(p: ca.Point) -> str:
    sd = p.species.get(species)
    if sd is None:
        return '<p class="none">This species is not present at this point.</p>'
    head = (
        f'<div class="meta">'
        f'<b>Grid point {p.index}</b> &nbsp;·&nbsp; '
        f'<i>x</i> = {p.x_pc:.4g} pc &nbsp;·&nbsp; '
        f'A<sub>V</sub> = {p.Av:.4g} mag &nbsp;·&nbsp; '
        f'n<sub>H</sub> = {p.nH:.4g} cm⁻³ &nbsp;·&nbsp; '
        f'T<sub>gas</sub> = {p.Tgas:.4g} K<br>'
        f'FUV = {p.FUV:.4g} &nbsp;·&nbsp; CR = {p.CR:.4g} s⁻¹'
        f'</div>'
    )
    summary = (
        f'<div class="summary">'
        f'<span>Abundance <b>{sd.abundance:.4e}</b></span>'
        f'<span>Formation rate <b>{sd.ptot:.3e}</b> s⁻¹</span>'
        f'<span>Destruction rate <b>{sd.dtot:.3e}</b> s⁻¹</span>'
        f'</div>'
    )
    tables = (
        '<div class="cols">'
        '<table><thead><tr><th>Formation</th><th>%</th></tr></thead>'
        f'<tbody>{_react_rows(sd.formation, "form")}</tbody></table>'
        '<table><thead><tr><th>Destruction</th><th>%</th></tr></thead>'
        f'<tbody>{_react_rows(sd.destruction, "dest")}</tbody></table>'
        '</div>'
    )
    return f'<h4>{html.escape(species)} at grid point {p.index}</h4>{head}{summary}{tables}'


info = [_point_html(p) for p in pts_sorted]
default_pt = next((i for i, (xv, yv) in enumerate(zip(Xs, Ys))
                   if np.isfinite(xv) and xv > 0 and np.isfinite(yv) and yv > 0), 0)

# --- Plotly figure (as JSON for the inline JS) ------------------------------
fig_data = [{
    "type": "scatter", "mode": "lines+markers", "name": species,
    "x": _logsafe(Xs), "y": _logsafe(Ys),
    "line": {"color": "#1f77b4", "width": 2},
    "marker": {"size": 6, "color": "#1f77b4"},
    "hovertemplate": (f"{x_short}=%{{x:.3g}}<br>n({species})/n_H=%{{y:.3e}}"
                      "<extra></extra>"),
}]
fig_layout = {
    "height": 430,
    "margin": {"t": 50, "l": 85, "r": 25, "b": 60},
    "template": "plotly_white",
    # hovermode "x" lets the cursor sit at ANY y for a given x and still update
    # the table — the user no longer has to track the curve exactly.
    "hovermode": "x",
    "title": {"text": f"{html.escape(prefix)} — n({html.escape(species)})/n<sub>H</sub>",
              "x": 0.5, "xanchor": "center"},
    "xaxis": {"type": "log", "title": {"text": XLABEL}, "exponentformat": "e",
              "showexponent": "all", "gridcolor": "rgba(0,0,0,0.08)",
              "showspikes": True, "spikemode": "across", "spikethickness": 1,
              "spikecolor": "#888", "spikedash": "dot"},
    "yaxis": {"type": "log", "title": {"text": YLABEL}, "exponentformat": "e",
              "showexponent": "all", "gridcolor": "rgba(0,0,0,0.08)"},
}
if x_trans is not None and np.isfinite(x_trans) and x_trans > 0:
    fig_layout["shapes"] = [{
        "type": "line", "x0": float(x_trans), "x1": float(x_trans),
        "yref": "paper", "y0": 0, "y1": 1,
        "line": {"color": "gray", "width": 1, "dash": "dash"},
    }]
    fig_layout["annotations"] = [{
        "x": float(x_trans), "y": 1, "yref": "paper", "yanchor": "bottom",
        "text": "HI→H₂", "showarrow": False, "font": {"size": 11, "color": "gray"},
    }]

# theme-aware colours: in dark mode use white text + brighter green/red so the
# pathway table stays readable (the Plotly figure itself stays white).
try:
    _DARK = st.context.theme.type == "dark"
except Exception:  # noqa: BLE001
    _DARK = False
_C = (dict(txt="#e6e6e6", meta="#c9c9c9", green="#4ade80", red="#f87171",
           sumbg="rgba(255,255,255,0.08)", sumb="#93c5fd", th="rgba(255,255,255,0.30)",
           td="rgba(255,255,255,0.14)", none="#9aa0a6", hint="#aab")
      if _DARK else
      dict(txt="#262730", meta="#444", green="#15803d", red="#b91c1c",
           sumbg="#f3f6fb", sumb="#0b3d91", th="#ccc", td="#eee", none="#999", hint="#888"))

# --- self-contained HTML component: plot + auto-updating pathway table -------
html_doc = """
<!DOCTYPE html><html><head><meta charset="utf-8">
<script>__PLOTLYJS__</script>
<style>
  body{margin:0;font-family:"Source Sans Pro",system-ui,sans-serif;color:__TXT__;}
  #plot{width:100%;}
  #panel{padding:6px 10px 12px;}
  #panel h4{margin:.4em 0 .3em;font-size:1.0rem;color:__TXT__;}
  #panel .meta{font-size:.84rem;color:__META__;line-height:1.5;}
  #panel .summary{display:flex;gap:18px;flex-wrap:wrap;margin:.5em 0;
     font-size:.9rem;background:__SUMBG__;border-radius:6px;padding:6px 12px;color:__TXT__;}
  #panel .summary b{color:__SUMB__;}
  #panel .cols{display:flex;gap:20px;flex-wrap:wrap;}
  #panel table{border-collapse:collapse;flex:1 1 320px;font-size:.86rem;color:__TXT__;}
  #panel th{text-align:left;border-bottom:2px solid __TH__;padding:3px 8px;color:__TXT__;}
  #panel td{padding:2px 8px;border-bottom:1px solid __TD__;vertical-align:top;}
  #panel td.eq{font-family:ui-monospace,Menlo,Consolas,monospace;white-space:nowrap;}
  #panel td.pct{text-align:right;font-variant-numeric:tabular-nums;font-weight:600;}
  #panel td.pct.form{color:__GREEN__;}
  #panel td.pct.dest{color:__RED__;}
  #panel td.none,#panel p.none{color:__NONE__;font-style:italic;}
  #hint{font-size:.8rem;color:__HINT__;padding:0 10px;}
</style></head><body>
<div id="plot"></div>
<div id="hint">Hover over the abundance curve to inspect a depth point.</div>
<div id="panel">__DEFAULT__</div>
<script>
  var INFO = __INFO__;
  var gd = document.getElementById('plot');
  Plotly.newPlot(gd, __DATA__, __LAYOUT__, {responsive:true, displaylogo:false});
  function show(i){ if(i!=null && INFO[i]!=null){ document.getElementById('panel').innerHTML = INFO[i]; } }
  gd.on('plotly_hover', function(ev){
     if(ev.points && ev.points.length){ show(ev.points[0].pointNumber); }
  });
  show(__DEFAULTIDX__);
</script></body></html>
"""
for _tok, _key in [("__TXT__", "txt"), ("__META__", "meta"), ("__GREEN__", "green"),
                   ("__RED__", "red"), ("__SUMBG__", "sumbg"), ("__SUMB__", "sumb"),
                   ("__TH__", "th"), ("__TD__", "td"), ("__NONE__", "none"), ("__HINT__", "hint")]:
    html_doc = html_doc.replace(_tok, _C[_key])
html_doc = (html_doc
            .replace("__PLOTLYJS__", _plotlyjs())
            .replace("__DATA__", json.dumps(fig_data))
            .replace("__LAYOUT__", json.dumps(fig_layout))
            .replace("__INFO__", json.dumps(info))
            .replace("__DEFAULT__", info[default_pt] if info else "")
            .replace("__DEFAULTIDX__", str(default_pt)))

st.iframe(html_doc, height=860)

st.caption(
    f"Abundance of **{species}** versus ${x_short}$ (log–log). "
    "Move the cursor along the curve to update the pathway report below the plot "
    "to that depth point. Percentages are each reaction's share of the total "
    "formation (green) / destruction (red) rate, listed until they reach 100%."
    + (f"  HI→H₂ transition at Aᵥ = {trans_av:.4g} mag." if trans_av else "")
    + ("  (Uniform-density cloud — x-axis fixed to A_V.)" if uniform else ""))
