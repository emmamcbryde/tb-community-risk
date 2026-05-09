from __future__ import annotations

import streamlit as st

from app.state import get_backend, init_session_state


st.set_page_config(page_title="TB Community Risk Models", layout="wide")

init_session_state()
backend = get_backend()

with st.sidebar:
    st.header("TB Community Risk")
    status = backend.status()
    st.caption(f"APY MATLAB backend: {status['name']}")
    st.caption(f"MATLAB started: {status['started']}")
    if status.get("error"):
        st.error(status["error"])

pages = {
    "APY v9 ABM": [
        st.Page("pages/1_Scenario.py", title="Scenario"),
        st.Page("pages/2_Run_Model.py", title="Run Model"),
        st.Page("pages/3_Results.py", title="Results"),
        st.Page("pages/4_Economics.py", title="Economics"),
    ],
    "Dynamic Model": [
        st.Page("pages/5_Dynamic_Model.py", title="Dynamic Workflow"),
    ],
}

navigation = st.navigation(pages)
navigation.run()
