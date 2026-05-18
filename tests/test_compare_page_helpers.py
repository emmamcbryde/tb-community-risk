from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import types
from unittest.mock import MagicMock


def load_compare_page_module(monkeypatch):
    root = Path(__file__).resolve().parents[1]
    monkeypatch.syspath_prepend(str(root))

    fake_streamlit = types.ModuleType("streamlit")
    fake_streamlit.session_state = {}
    fake_streamlit.title = MagicMock()
    fake_streamlit.caption = MagicMock()
    fake_streamlit.columns = MagicMock(
        return_value=[
            types.SimpleNamespace(button=MagicMock(return_value=False)),
            types.SimpleNamespace(button=MagicMock(return_value=False)),
            types.SimpleNamespace(button=MagicMock(return_value=False)),
            types.SimpleNamespace(button=MagicMock(return_value=False)),
        ]
    )
    fake_streamlit.info = MagicMock()
    fake_streamlit.stop = MagicMock(side_effect=RuntimeError("st.stop"))

    fake_app_state = types.ModuleType("app.state")
    fake_app_state.init_session_state = MagicMock()
    fake_app_state.get_backend = MagicMock(return_value=types.SimpleNamespace())
    fake_app_state.record_message = MagicMock()
    fake_app_state.sync_backend_status = MagicMock()

    monkeypatch.setitem(sys.modules, "streamlit", fake_streamlit)
    monkeypatch.setitem(sys.modules, "app.state", fake_app_state)

    path = root / "pages" / "6_Compare.py"
    spec = importlib.util.spec_from_file_location("compare_page_for_tests", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None

    try:
        spec.loader.exec_module(module)
    except RuntimeError as exc:
        assert str(exc) == "st.stop"

    return module


def test_compare_economics_rows_accepts_python_lowercase_summary_rows(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)
    baseline_rows = [
        {
            "metric": "testingCost",
            "value": 100.0,
            "category": "directCost",
            "status": "implemented",
            "includedInTotal": True,
        },
        {
            "metric": "totalImplementedCost",
            "value": 300.0,
            "category": "directCostTotal",
            "status": "implemented",
            "includedInTotal": False,
        },
    ]
    comparator_rows = [
        {
            "metric": "testingCost",
            "value": 130.0,
            "category": "directCost",
            "status": "implemented",
            "includedInTotal": True,
        },
        {
            "metric": "totalImplementedCost",
            "value": 360.0,
            "category": "directCostTotal",
            "status": "implemented",
            "includedInTotal": False,
        },
    ]

    rows, warnings = module.compare_economics_rows(baseline_rows, comparator_rows)

    assert warnings == []
    by_metric = {row["metric"]: row for row in rows}
    assert by_metric["testingCost"]["baseline"] == 100.0
    assert by_metric["testingCost"]["comparator"] == 130.0
    assert by_metric["testingCost"]["absoluteDifference"] == 30.0
    assert by_metric["testingCost"]["relativeDifference"] == 0.3
    assert by_metric["testingCost"]["category"] == "directCost"
    assert by_metric["testingCost"]["status"] == "implemented"
    assert by_metric["testingCost"]["includedInTotal"] is True
    assert by_metric["totalImplementedCost"]["absoluteDifference"] == 60.0
    assert by_metric["totalImplementedCost"]["relativeDifference"] == 0.2
    assert by_metric["totalImplementedCost"]["includedInTotal"] is False


def test_compare_economics_rows_can_key_by_component(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)
    baseline_rows = [
        {
            "component": "programSetupCost",
            "value": 50.0,
            "category": "directCost",
            "status": "implemented",
            "includedInTotal": True,
        }
    ]
    comparator_rows = [
        {
            "component": "programSetupCost",
            "value": 75.0,
            "category": "directCost",
            "status": "implemented",
            "includedInTotal": True,
        }
    ]

    rows, warnings = module.compare_economics_rows(baseline_rows, comparator_rows)

    assert warnings == []
    assert rows == [
        {
            "metric": "programSetupCost",
            "baseline": 50.0,
            "comparator": 75.0,
            "absoluteDifference": 25.0,
            "relativeDifference": 0.5,
            "category": "directCost",
            "status": "implemented",
            "includedInTotal": True,
        }
    ]


def test_compare_economics_rows_retains_keyed_rows_without_numeric_values(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)
    baseline_rows = [
        {
            "component": "tbDiseaseCostsAverted",
            "status": "unsupported",
            "category": "diseaseCost",
            "includedInTotal": False,
        },
        {
            "metric": "programSetupCost",
            "value": None,
            "status": "missing",
            "category": "directCost",
            "includedInTotal": True,
        },
    ]
    comparator_rows = [
        {
            "component": "tbDiseaseCostsAverted",
            "status": "unsupported",
            "category": "diseaseCost",
            "includedInTotal": False,
        },
        {
            "metric": "programSetupCost",
            "value": 75.0,
            "status": "implemented",
            "category": "directCost",
            "includedInTotal": True,
        },
    ]

    rows, warnings = module.compare_economics_rows(baseline_rows, comparator_rows)

    assert warnings == []
    by_metric = {row["metric"]: row for row in rows}
    assert by_metric["tbDiseaseCostsAverted"] == {
        "metric": "tbDiseaseCostsAverted",
        "baseline": None,
        "comparator": None,
        "absoluteDifference": None,
        "relativeDifference": None,
        "category": "diseaseCost",
        "status": "unsupported",
        "includedInTotal": False,
    }
    assert by_metric["programSetupCost"] == {
        "metric": "programSetupCost",
        "baseline": None,
        "comparator": 75.0,
        "absoluteDifference": None,
        "relativeDifference": None,
        "category": "directCost",
        "status": "missing",
        "includedInTotal": True,
    }


def test_compare_outcome_rows_still_skips_rows_without_values(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)
    baseline_rows = [
        {"Metric": "nScreened"},
        {"Metric": "nTotalCoursesStarted", "Median": 40.0},
    ]
    comparator_rows = [
        {"Metric": "nScreened", "Median": 125.0},
        {"Metric": "nTotalCoursesStarted", "Median": 50.0},
    ]

    rows, warnings = module.compare_outcome_rows(baseline_rows, comparator_rows)

    assert [row["metric"] for row in rows] == ["nTotalCoursesStarted", "nScreened"]
    assert rows[0]["absoluteDifference"] == 10.0
    assert warnings == [
        "Skipped rows missing Metric + Median: baseline=1, comparator=0."
    ]
