"""Entry point for the general community TB screening application.

Run with ``streamlit run general_app.py``. The frozen setting-specific workflow
remains available through ``streamlit_app.py`` and the release tag.
"""

from __future__ import annotations

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.import_bootstrap import ensure_repo_root  # noqa: E402

ensure_repo_root(REPO_ROOT)

import streamlit as st  # noqa: E402

from app.general.state import init_general_state  # noqa: E402
from app.general.terminology import APP_TITLE, WORKFLOW_PAGES  # noqa: E402

st.set_page_config(page_title=APP_TITLE, layout="wide")
init_general_state()

navigation = st.navigation(
    {APP_TITLE: [st.Page(path, title=title) for path, title in WORKFLOW_PAGES]}
)
navigation.run()
