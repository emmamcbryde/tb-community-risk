from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path
from unittest.mock import MagicMock


def test_economics_page_import_reaches_empty_config_stop(monkeypatch) -> None:
    fake_streamlit = types.ModuleType("streamlit")
    fake_streamlit.session_state = {}
    fake_streamlit.title = MagicMock()
    fake_streamlit.caption = MagicMock()
    fake_streamlit.info = MagicMock()
    fake_streamlit.stop = MagicMock(side_effect=RuntimeError("st.stop"))
    fake_streamlit.columns = MagicMock(
        return_value=[
            types.SimpleNamespace(button=MagicMock(return_value=False)),
            types.SimpleNamespace(button=MagicMock(return_value=False)),
        ]
    )

    fake_app_state = types.ModuleType("app.state")
    fake_app_state.init_session_state = MagicMock()
    fake_app_state.get_backend = MagicMock(return_value=types.SimpleNamespace())
    fake_app_state.get_backend_name = MagicMock(return_value="python_apy")
    fake_app_state.mark_economics_changed = MagicMock()
    fake_app_state.mark_economics_completed = MagicMock()
    fake_app_state.record_message = MagicMock()
    fake_app_state.sync_backend_status = MagicMock()

    monkeypatch.setitem(sys.modules, "streamlit", fake_streamlit)
    monkeypatch.setitem(sys.modules, "app.state", fake_app_state)

    path = Path(__file__).resolve().parents[1] / "pages" / "4_Economics.py"
    spec = importlib.util.spec_from_file_location(
        "economics_page_import_for_tests",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None

    try:
        spec.loader.exec_module(module)
    except RuntimeError as exc:
        assert str(exc) == "st.stop"
    else:
        assert True

    fake_streamlit.stop.assert_called_once()
