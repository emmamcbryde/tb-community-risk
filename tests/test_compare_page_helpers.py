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


def test_attributable_risk_payload_for_bundle_calls_backend_and_json_safes(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)
    bundle = {"technical": {"dynamicComparison": {"activeTbCases": 4}}}

    class FakeBackend:
        def __init__(self) -> None:
            self.seen_bundle = None

        def run_attributable_risk(self, result_bundle):
            self.seen_bundle = result_bundle
            return {
                "status": "ok",
                "calculatedRows": [
                    {
                        "metric": "ExpectedAttributableCases20y_Per1500",
                        "value": float("inf"),
                    }
                ],
                "missingInputs": (),
                "unsupportedMetrics": [],
                "messages": ["calculated"],
            }

    backend = FakeBackend()

    payload = module.attributable_risk_payload_for_bundle(backend, bundle)

    assert backend.seen_bundle is bundle
    assert payload == {
        "status": "ok",
        "calculatedRows": [
            {
                "metric": "ExpectedAttributableCases20y_Per1500",
                "value": None,
            }
        ],
        "missingInputs": [],
        "unsupportedMetrics": [],
        "messages": ["calculated"],
    }


def test_attributable_risk_payload_for_bundle_reports_unsupported_backend(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)

    payload = module.attributable_risk_payload_for_bundle(types.SimpleNamespace(), {"id": "bundle"})

    assert payload["status"] == "unsupported"
    assert payload["source"] == "compare_page_backend_capability_check"
    assert payload["calculatedRows"] == []
    assert payload["missingInputs"] == []
    assert payload["unsupportedMetrics"] == [
        {
            "metric": "attributableRisk",
            "reason": "The selected backend does not expose run_attributable_risk.",
            "backendCapability": "run_attributable_risk",
        }
    ]
    assert payload["messages"] == [
        (
            "Attributable-risk comparison is unsupported for the selected "
            "backend because it does not expose run_attributable_risk."
        )
    ]


def test_store_attributable_risk_compare_payloads_writes_session_state(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)
    baseline_bundle = {"scenarioLabel": "baseline"}
    comparator_bundle = {"scenarioLabel": "comparator"}

    class FakeBackend:
        def run_attributable_risk(self, result_bundle):
            return {
                "status": "ok",
                "calculatedRows": [
                    {"metric": result_bundle["scenarioLabel"], "value": 1.0}
                ],
                "missingInputs": [],
                "unsupportedMetrics": [],
                "messages": [result_bundle["scenarioLabel"]],
            }

    baseline_payload, comparator_payload = module.store_attributable_risk_compare_payloads(
        FakeBackend(),
        baseline_bundle,
        comparator_bundle,
    )

    assert baseline_payload["messages"] == ["baseline"]
    assert comparator_payload["messages"] == ["comparator"]
    assert (
        module.st.session_state["compare_baseline_attributable_risk_payload"]
        is baseline_payload
    )
    assert (
        module.st.session_state["compare_comparator_attributable_risk_payload"]
        is comparator_payload
    )


def test_reset_compare_outputs_clears_attributable_risk_payloads(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)
    module.st.session_state["compare_baseline_attributable_risk_payload"] = {
        "status": "ok"
    }
    module.st.session_state["compare_comparator_attributable_risk_payload"] = {
        "status": "unsupported"
    }

    assert module.has_compare_outputs() is True

    module.reset_compare_outputs(outputs_cleared=True)

    assert module.st.session_state["compare_baseline_attributable_risk_payload"] is None
    assert module.st.session_state["compare_comparator_attributable_risk_payload"] is None
    assert module.st.session_state["compare_outputs_cleared"] is True


def test_attributable_risk_compare_rows_preserves_missing_input_status(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)
    payload = {
        "status": "missing-input",
        "calculatedRows": [],
        "missingInputs": [
            {
                "field": "technical.dynamicComparison",
                "message": "technical.dynamicComparison is required.",
            }
        ],
        "unsupportedMetrics": [],
        "messages": ["Missing technical.dynamicComparison."],
    }

    rows = module.attributable_risk_compare_rows(payload)

    by_field = {}
    for row in rows:
        by_field.setdefault(row["payloadField"], []).append(row)
    assert by_field["status"] == [
        {
            "payloadField": "status",
            "rowIndex": None,
            "status": "missing-input",
            "empty": False,
            "itemJson": None,
        }
    ]
    assert by_field["calculatedRows"][0]["empty"] is True
    assert by_field["missingInputs"][0]["status"] == "missing-input"
    assert by_field["missingInputs"][0]["field"] == "technical.dynamicComparison"
    assert by_field["missingInputs"][0]["message"] == "technical.dynamicComparison is required."
    assert by_field["unsupportedMetrics"][0]["empty"] is True
    assert by_field["messages"][0]["message"] == "Missing technical.dynamicComparison."


def test_attributable_risk_compare_rows_preserves_unsupported_metrics(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)
    payload = {
        "status": "unsupported",
        "calculatedRows": [],
        "missingInputs": [],
        "unsupportedMetrics": [
            {
                "metric": "PopulationAttributableFraction20y_ReactivationOnly",
                "reason": "Not implemented in the Python APY port yet.",
                "matlabSource": "run_tb_screening_reactivation_attributable_v9",
            }
        ],
        "messages": ["Reactivation attributable-risk metrics are unsupported."],
    }

    rows = module.attributable_risk_compare_rows(payload)

    unsupported_rows = [
        row for row in rows if row["payloadField"] == "unsupportedMetrics"
    ]
    assert unsupported_rows == [
        {
            "payloadField": "unsupportedMetrics",
            "rowIndex": 0,
            "status": "unsupported",
            "empty": False,
            "itemJson": (
                '{"matlabSource": "run_tb_screening_reactivation_attributable_v9", '
                '"metric": "PopulationAttributableFraction20y_ReactivationOnly", '
                '"reason": "Not implemented in the Python APY port yet."}'
            ),
            "metric": "PopulationAttributableFraction20y_ReactivationOnly",
            "reason": "Not implemented in the Python APY port yet.",
            "matlabSource": "run_tb_screening_reactivation_attributable_v9",
        }
    ]


def test_attributable_risk_compare_rows_preserves_future_calculated_rows(monkeypatch) -> None:
    module = load_compare_page_module(monkeypatch)
    payload = {
        "status": "ok",
        "calculatedRows": [
            {
                "metric": "ExpectedAttributableCases20y_Per1500",
                "value": 12.5,
                "units": "cases per 1500",
            }
        ],
        "missingInputs": [],
        "unsupportedMetrics": [],
        "messages": [],
    }

    rows = module.attributable_risk_compare_rows(payload)

    calculated_rows = [row for row in rows if row["payloadField"] == "calculatedRows"]
    assert calculated_rows == [
        {
            "payloadField": "calculatedRows",
            "rowIndex": 0,
            "status": "ok",
            "empty": False,
            "itemJson": (
                '{"metric": "ExpectedAttributableCases20y_Per1500", '
                '"units": "cases per 1500", "value": 12.5}'
            ),
            "metric": "ExpectedAttributableCases20y_Per1500",
            "value": 12.5,
            "units": "cases per 1500",
        }
    ]
    empty_fields = {
        row["payloadField"]
        for row in rows
        if row["payloadField"] in {"missingInputs", "unsupportedMetrics", "messages"}
        and row["empty"]
    }
    assert empty_fields == {"missingInputs", "unsupportedMetrics", "messages"}
