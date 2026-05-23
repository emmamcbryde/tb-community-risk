from __future__ import annotations

import importlib
import json
import sys
import builtins
from types import SimpleNamespace

from adapters.python_apy_backend import PythonApyBackend


def _import_app_state_with_fake_streamlit(monkeypatch):
    session_state: dict[str, object] = {}

    def cache_resource(*args, **kwargs):
        def decorator(func):
            return func

        return decorator

    fake_streamlit = SimpleNamespace(
        cache_resource=cache_resource,
        session_state=session_state,
    )
    monkeypatch.setitem(sys.modules, "streamlit", fake_streamlit)
    sys.modules.pop("app.state", None)
    return importlib.import_module("app.state"), session_state


def test_init_session_state_defaults_to_python_apy(monkeypatch) -> None:
    app_state, session_state = _import_app_state_with_fake_streamlit(monkeypatch)

    app_state.init_session_state()

    assert session_state["apy_backend_name"] == "python_apy"
    assert app_state.get_backend_name() == "python_apy"
    assert session_state["backend_status"]["name"] == "python_apy"


def test_get_backend_name_defaults_to_python_apy_without_init(monkeypatch) -> None:
    app_state, session_state = _import_app_state_with_fake_streamlit(monkeypatch)

    assert "apy_backend_name" not in session_state
    assert app_state.get_backend_name() == "python_apy"


def test_default_get_backend_returns_python_backend_without_matlab_engine(
    monkeypatch,
) -> None:
    app_state, _session_state = _import_app_state_with_fake_streamlit(monkeypatch)
    sys.modules.pop("matlab", None)
    sys.modules.pop("matlab.engine", None)
    original_import = builtins.__import__

    def fail_on_matlab_engine_import(name, *args, **kwargs):
        if name == "matlab.engine":
            raise AssertionError("get_backend() attempted to import matlab.engine")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fail_on_matlab_engine_import)

    app_state.init_session_state()
    backend = app_state.get_backend()

    assert isinstance(backend, PythonApyBackend)
    assert app_state.get_backend_name() == "python_apy"
    assert "matlab.engine" not in sys.modules


def test_get_apy_status_rows_returns_table_rows_without_matlab_import(
    monkeypatch,
) -> None:
    sys.modules.pop("matlab", None)
    sys.modules.pop("matlab.engine", None)
    original_import = builtins.__import__

    def fail_on_matlab_import(name, *args, **kwargs):
        if name in {"matlab", "matlab.engine"}:
            raise AssertionError(
                "get_apy_status_rows() attempted to import MATLAB modules"
            )
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fail_on_matlab_import)
    app_state, _session_state = _import_app_state_with_fake_streamlit(monkeypatch)

    rows = app_state.get_apy_status_rows()

    assert rows
    assert all(isinstance(row, dict) for row in rows)
    assert {
        "Component",
        "Migration status",
        "Coverage status",
        "MATLAB reference",
        "Reference fixtures",
        "Parity tests",
        "Distributional tests",
        "Python tests",
    }.issubset(rows[0])
    json.dumps(rows, allow_nan=False, sort_keys=True)
    assert "matlab" not in sys.modules
    assert "matlab.engine" not in sys.modules
