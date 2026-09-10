from __future__ import annotations

import streamlit as st

from app.state import init_session_state


init_session_state()

st.title("Set up")
st.info(
    "Population, strategy and parameter editing have moved to the unified Set up page. "
    "Use that page before running an analysis."
)
st.page_link("pages/0_Start.py", label="Open Set up")
