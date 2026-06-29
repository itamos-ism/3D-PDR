"""Feature — Code Configuration.

Edit the compile-time flags in ``src/config.mk`` and compile the code. The most-
used flags are always visible; the rest live under an "Advanced" expander. A
single "Save & Compile" button writes the flags and builds; the compile log
collapses when it finishes (still expandable). "Next" goes to Model Parameters.
"""
from __future__ import annotations

import streamlit as st

from pdrstudio import config, paths, runner, ui

ui.configure_page("Code Configuration")

_MODEL_PARAMS = "pages/2_📐_Model_Parameters.py"


def _go_model_params():
    st.session_state["mp_section"] = "params.dat"
    st.switch_page(_MODEL_PARAMS)


st.title("⚙️ Code Configuration")

cfg_path = paths.config_mk()
if not cfg_path.exists():
    st.error(f"config.mk not found at `{cfg_path}`.")
    st.stop()

current = config.read_config()
exe = paths.executable()

# "Editing…" note directly under the title
st.caption(
    f"Editing `{cfg_path}` · modified {ui.mtime_str(cfg_path)}  \n"
    f"Executable: {ui.status_badge(exe.exists(), 'built', 'not built')} · {ui.mtime_str(exe)}")


def _render_flag(flag: config.Flag, cur_values: dict, col) -> str:
    cur = cur_values.get(flag.key, "")
    label = f"{flag.label}  ·  `{flag.key}`"
    if flag.kind == "bool":
        on = col.toggle(label, value=str(cur).strip() == "1", help=flag.help, key=f"cfg_{flag.key}")
        return "1" if on else "0"
    if flag.kind == "choice":
        opts = flag.choices
        idx = opts.index(cur) if cur in opts else 0
        return col.selectbox(label, opts, index=idx, help=flag.help, key=f"cfg_{flag.key}")
    return col.text_input(label, value=cur, help=flag.help, key=f"cfg_{flag.key}")


# --- action row: Save & Compile + Next (same height) -----------------------
a1, a2, a3, a4 = st.columns([1.6, 1.2, 2.2, 1.4])
save_compile = a1.button("💾 Save & Compile", type="primary", width="stretch")
clean = a2.checkbox("Clean build", value=True, help="Run `make clean` first.")
if a4.button("Next ➡️", key="cc_next_top", type="primary", width="stretch",
             help="Go to Model Parameters"):
    _go_model_params()

status_slot = st.container()  # save feedback + compile log appear here

st.caption(
    "Flags are read live from `src/config.mk` and `src/makefile`, so PDR-studio "
    "adapts to whatever flags this 3D-PDR version defines — new flags are imported "
    "automatically.")

new_flags = config.newly_discovered()
if new_flags:
    with st.expander(f"🆕 {len(new_flags)} new flag(s) imported from this version's "
                     "config.mk", expanded=True):
        st.caption("These flags exist in this code version but PDR-studio has no curated "
                   "label/grouping for them yet. They are shown below in the **Advanced** "
                   "area and are safe to edit and save.")
        for f in new_flags:
            kind = (f"choice: {', '.join(f.choices)}" if f.kind == "choice" else f.kind)
            st.markdown(f"- **`{f.key}`** *(section: {f.section}; {kind})* — "
                        f"{f.help or 'no description in config.mk'}")
st.divider()

# --- flags ------------------------------------------------------------------
new_values: dict[str, str] = {}
for group, flags in config.visible_groups().items():
    st.subheader(group)
    cols = st.columns(2)
    for i, flag in enumerate(flags):
        new_values[flag.key] = _render_flag(flag, current, cols[i % 2])
    st.divider()

with st.expander("⚙️ Advanced flags (change only if you know what you are doing)"):
    for group, flags in config.advanced_groups().items():
        st.markdown(f"**{group}**")
        cols = st.columns(2)
        for i, flag in enumerate(flags):
            new_values[flag.key] = _render_flag(flag, current, cols[i % 2])
        st.markdown("")

with st.expander("View raw config.mk"):
    st.code(cfg_path.read_text(), language="makefile")

# --- bottom-right Next ------------------------------------------------------
_, bn = st.columns([4, 1])
if bn.button("Next ➡️", key="cc_next_bottom", type="primary", width="stretch",
             help="Go to Model Parameters"):
    _go_model_params()

# --- handle Save & Compile --------------------------------------------------
if save_compile:
    merged = {**current, **new_values}
    notes = config.normalize_dependencies(merged)
    config.write_config(merged)
    with status_slot:
        st.success("config.mk saved.")
        for note in notes:
            st.info(note, icon="🔗")
        log_lines: list[str] = []
        exit_code: int | None = None
        with st.status("Compiling 3D-PDR…", expanded=True) as status:
            placeholder = st.empty()
            for line in runner.compile_code(clean=clean):
                code = runner.parse_exit(line)
                if code is not None:
                    exit_code = code
                    break
                log_lines.append(line)
                placeholder.code("\n".join(log_lines[-400:]), language="text")
            if exit_code == 0:
                # collapse the log once finished (the user can re-expand it)
                status.update(label="✅ Compilation succeeded", state="complete", expanded=False)
            else:
                status.update(label=f"❌ Compilation failed (exit {exit_code})", state="error")
        if exit_code != 0:
            st.error("Build failed — expand the log above to see why.")

        if exit_code == 0:
            # Also build the RT-tool so its binary stays in sync.
            rt_log: list[str] = []
            rt_exit: int | None = None
            with st.status("Building RT-tool…", expanded=True) as rt_status:
                rt_ph = st.empty()
                for line in runner.compile_rttool():
                    code = runner.parse_exit(line)
                    if code is not None:
                        rt_exit = code
                        break
                    rt_log.append(line)
                    rt_ph.code("\n".join(rt_log[-300:]), language="text")
                if rt_exit == 0:
                    rt_status.update(label="✅ RT-tool built", state="complete", expanded=False)
                else:
                    rt_status.update(label=f"❌ RT-tool build failed (exit {rt_exit})", state="error")

            if rt_exit == 0:
                # Invalidate cached RT-tool results so the RT-tool page re-runs.
                st.session_state.pop("rt_prefix", None)
                st.session_state.pop("rt_stdout", None)
                st.success(f"Build complete. Executable: `{exe}`  ·  press **Next ➡️** to set up the model.")
            else:
                st.error("RT-tool build failed — expand the log above to see why.")
