from __future__ import annotations

import streamlit as st

from app.state import get_backend, init_session_state


st.set_page_config(page_title="APY v9 TB Screening", layout="wide")

init_session_state()
backend = get_backend()

with st.sidebar:
    st.header("APY v9")
    status = backend.status()
    st.caption(f"Backend: {status['name']}")
    st.caption(f"Started: {status['started']}")
    if status.get("error"):
        st.error(status["error"])

pages = [
    st.Page("pages/1_Scenario.py", title="Scenario"),
    st.Page("pages/2_Run_Model.py", title="Run Model"),
    st.Page("pages/3_Results.py", title="Results"),
    st.Page("pages/4_Economics.py", title="Economics"),
]

navigation = st.navigation(pages)
navigation.run()
