from __future__ import annotations

from pathlib import Path
import sys
REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.import_bootstrap import ensure_repo_root

ensure_repo_root(REPO_ROOT)

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
from app.state import init_session_state

# ---------------------------------------------------------------------
# Page config (must come before any UI)
# ---------------------------------------------------------------------
st.set_page_config(
    page_title="LTBI Screening Decision Tool",
    layout="wide",
)

# ---------------------------------------------------------------------
# Initialize session state
# ---------------------------------------------------------------------
init_session_state()

# ---------------------------------------------------------------------
# Navigation structure
# ---------------------------------------------------------------------
navigation = st.navigation(
    {
        "LTBI Screening Tool": [
            st.Page("pages/0_Start.py", title="Start"),
            st.Page("pages/1_Scenario.py", title="Define Strategy"),
            st.Page("pages/2_Run_Model.py", title="Run Analysis"),
            st.Page("pages/3_Results.py", title="Results"),
            st.Page("pages/5_Decision_Analysis.py", title="Explore Decisions"),
            st.Page("pages/4_Economics.py", title="Evidence & Assumptions"),
        ],
        "Research & Development": [
            st.Page("pages/8_Technical_Settings.py", title="Technical settings"),
            st.Page("pages/5_Dynamic_Model.py", title="Dynamic transmission model"),
            st.Page("pages/7_Dynamic_ABM_Compare.py", title="Model comparison"),
        ],
    }
)

# ---------------------------------------------------------------------
# Run selected page
# ---------------------------------------------------------------------
navigation.run()
