from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.import_bootstrap import ensure_repo_root

ensure_repo_root(REPO_ROOT)

from __future__ import annotations

from pathlib import Path
import sys
import streamlit as st

# ---------------------------------------------------------------------
# Ensure repo root is on sys.path (needed for ui/, engine/, adapters/)
# ---------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# ---------------------------------------------------------------------
# App imports
# ---------------------------------------------------------------------
from app.state import get_backend, init_session_state

# ---------------------------------------------------------------------
# Page config (must come before any UI)
# ---------------------------------------------------------------------
st.set_page_config(
    page_title="TB Community Risk Models",
    layout="wide",
)

# ---------------------------------------------------------------------
# Initialize session + backend
# ---------------------------------------------------------------------
init_session_state()
backend = get_backend()

# ---------------------------------------------------------------------
# Sidebar (light, model-aware)
# ---------------------------------------------------------------------
with st.sidebar:
    st.markdown("### APY MATLAB backend")
    status = st.session_state.get("backend_status", "unknown")
    st.write(status)

# ---------------------------------------------------------------------
# Navigation structure
# ---------------------------------------------------------------------
navigation = st.navigation(
    {
        "APY v9 ABM": [
            st.Page("pages/1_Scenario.py", title="Scenario"),
            st.Page("pages/2_Run_Model.py", title="Run Model"),
            st.Page("pages/3_Results.py", title="Results"),
            st.Page("pages/4_Economics.py", title="Economics"),
            st.Page("pages/6_Compare.py", title="Compare"),
        ],
        "Dynamic Model": [
            st.Page("pages/5_Dynamic_Model.py", title="Dynamic Workflow"),
        ],
        "Integrated workflow": [
            st.Page("pages/7_Dynamic_ABM_Compare.py", title="Dynamic + ABM Compare"),
        ],
    }
)

# ---------------------------------------------------------------------
# Run selected page
# ---------------------------------------------------------------------
navigation.run()
