"""PDR-studio — a browser-based control panel for 3D-PDR.

Entry point: configures the page once, shows the logo/sidebar header, and builds
the sidebar navigation (Home + one page per feature) with custom titles & icons.

Launch with:   ./run.sh        (from the PDR-studio directory)
          or:   streamlit run app.py
"""
from __future__ import annotations

import streamlit as st

from pdrstudio import ui

st.set_page_config(page_title="PDR-studio", page_icon=ui.page_icon(), layout="wide")
ui.sidebar_header()

nav = st.navigation([
    st.Page("home.py", title="Home", icon="🏠", default=True),
    st.Page("pages/1_⚙️_Code_Configuration.py", title="Code Configuration", icon="⚙️"),
    st.Page("pages/2_📐_Model_Parameters.py", title="Model Parameters", icon="📐"),
    st.Page("pages/3_📡_Line_emission.py", title="Line emission", icon="📡"),
    st.Page("pages/4_📊_Abundances_Profiles.py", title="Abundance Profiles", icon="📊"),
    st.Page("pages/5_🔥_Heating_and_Cooling.py", title="Heating & Cooling", icon="🔥"),
    st.Page("pages/6_🧪_Chemical_Analysis.py", title="Chemical Analysis", icon="🧪"),
    st.Page("pages/7_🧬_Network_Builder.py", title="Network Builder", icon="🧬"),
    st.Page("pages/8_🐞_Report_a_bug.py", title="Report a bug", icon="🐞"),
])
nav.run()
