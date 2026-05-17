from __future__ import annotations

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.import_bootstrap import ensure_repo_root

ensure_repo_root(REPO_ROOT)

from app.state import init_session_state
from ui.dynamic_ui import render_dynamic_ui

init_session_state()
render_dynamic_ui()
