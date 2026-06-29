"""Feature — Report a bug.

Bundle the current code configuration, parameters, chemistry files, initial
conditions and the current model's outputs into a single tarball to email to the
3D-PDR team. The model is the one configured in params.dat (not selectable), and
its inputs/outputs are always included.
"""
from __future__ import annotations

import streamlit as st

from pdrstudio import params, report, ui

ui.configure_page("Report a bug")
st.title("🐞 Report a bug")

prefix = str(params.read_params(0).get("output", "model"))

st.markdown(
    "Describe the problem below. PDR-studio will package the current **code "
    "configuration**, **parameters**, **chemistry files**, your **initial "
    "conditions** and the **output files** of your current model "
    f"(`{prefix}`) into a single tarball — they are all included automatically.")

_MIN_CHARS = 50
note = st.text_area(
    "Describe the issue",
    placeholder="What happened? What did you expect? Any error messages? "
                "(at least 50 characters)",
    height=160,
)
n = len(note.strip())
if n < _MIN_CHARS:
    st.caption(f"✍️ Please write at least {_MIN_CHARS} characters — {n}/{_MIN_CHARS}.")

if st.button("📦 Build debug report", type="primary", disabled=n < _MIN_CHARS):
    path = report.build_report(note=note, model_prefix=prefix, include_outputs=True)
    st.success(f"Report created: `{path}`")
    with open(path, "rb") as fh:
        st.download_button("⬇️ Download tarball", fh, file_name=path.name,
                           mime="application/gzip")

st.markdown("Please email the tarball to **Thomas Bisbas (tbisbas@gmail.com)**.")
