from __future__ import annotations

import builtins
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock


def _cache_resource(*args, **kwargs):
    def decorator(func):
        return func

    return decorator


def _run_model_page_path() -> Path:
    return Path(__file__).resolve().parents[1] / "pages" / "2_Run_Model.py"


def _load_run_model_page(monkeypatch, session_state: dict[str, object]):
    fake_streamlit = SimpleNamespace(
        cache_resource=_cache_resource,
        session_state=session_state,
        title=MagicMock(),
        caption=MagicMock(),
        info=MagicMock(),
        warning=MagicMock(),
        success=MagicMock(),
        error=MagicMock(),
        subheader=MagicMock(),
        json=MagicMock(),
        button=MagicMock(return_value=False),
        stop=MagicMock(side_effect=RuntimeError("st.stop")),
    )

    monkeypatch.setitem(sys.modules, "streamlit", fake_streamlit)
    sys.modules.pop("app.state", None)

    spec = importlib.util.spec_from_file_location(
        "run_model_page_import_for_tests",
        _run_model_page_path(),
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module, fake_streamlit


def test_run_model_page_import_reaches_empty_config_stop_without_matlab(monkeypatch) -> None:
    sys.modules.pop("matlab.engine", None)
    original_import = builtins.__import__

    def fail_on_matlab_engine_import(name, *args, **kwargs):
        if name == "matlab.engine":
            raise AssertionError("Run Model page attempted to import matlab.engine")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fail_on_matlab_engine_import)

    try:
        _load_run_model_page(monkeypatch, {})
    except RuntimeError as exc:
        assert str(exc) == "st.stop"
    else:
        raise AssertionError("Run Model page should stop when no scenario is loaded.")

    assert "matlab.engine" not in sys.modules


def test_store_python_apy_export_helpers_is_stable_and_json_serialisable(
    monkeypatch,
) -> None:
    module, _fake_streamlit = _load_run_model_page(
        monkeypatch,
        {"config": {"scenarioLabel": "Example"}},
    )

    class FakeBackend:
        def status(self):
            return {"name": "python_apy", "started": True, "error": ""}

        def json_export_payload(self, bundle):
            return {"metadata": dict(bundle["metadata"])}

        def summary_csv_export(self, bundle):
            return {"available": True, "columns": ["Metric"], "rows": [{"Metric": "nScreened"}]}

        def headline_display_tables(self, bundle):
            return {"summaryRows": {"available": True, "columns": ["Metric"], "rows": []}}

        def chart_numeric_series(self, bundle):
            return {"available": True, "series": [{"name": "nScreened", "points": []}]}

    bundle = {"metadata": {"scenarioLabel": "Example"}}
    state_a: dict[str, object] = {}
    state_b: dict[str, object] = {}

    exports_a = module.store_python_apy_export_helpers(
        state_a,
        FakeBackend(),
        bundle,
        "python_apy",
    )
    exports_b = module.store_python_apy_export_helpers(
        state_b,
        FakeBackend(),
        bundle,
        "python_apy",
    )

    assert exports_a == exports_b
    assert state_a == state_b
    assert set(exports_a) == {
        "apy_json_export_payload",
        "apy_summary_csv_export",
        "apy_headline_display_tables",
        "apy_chart_numeric_series",
    }
    json.dumps(exports_a, allow_nan=False, sort_keys=True)
