from __future__ import annotations

import builtins
import json
import math
import sys
from pathlib import Path

from adapters.paths import repo_root
from adapters.python_apy_backend import PythonApyBackend
from engine.apy.economics import calculate_economics


FIXTURE_PATH = Path(__file__).parent / "fixtures" / "apy_economics_hand_check_fixture.json"
EXPECTED_ALIGNMENT = "MATLAB-formula-aligned, not MATLAB-exported"
UNSUPPORTED_ECONOMICS_FIELDS = ("DALYs", "QALYs", "discounting", "ICERs")


def _load_fixture() -> dict:
    with FIXTURE_PATH.open(encoding="utf-8") as fixture_file:
        return json.load(fixture_file)


def _assert_supported_subset_matches_fixture(payload: dict, fixture: dict) -> None:
    for section, expected_values in fixture["expectedSupportedSubset"].items():
        for metric, expected_value in expected_values.items():
            assert metric in payload[section]
            assert math.isclose(
                payload[section][metric],
                expected_value,
                rel_tol=0.0,
                abs_tol=1e-9,
            )


def test_hand_check_fixture_provenance_and_supported_outputs() -> None:
    fixture = _load_fixture()

    assert fixture["metadata"]["alignment"] == EXPECTED_ALIGNMENT

    payload = calculate_economics(
        fixture["resultBundle"],
        fixture["economicsConfig"],
    )

    assert payload["status"] == "partial"
    assert payload["coverage"]["status"] == "partial"
    _assert_supported_subset_matches_fixture(payload, fixture)


def test_hand_check_fixture_leaves_unsupported_economics_absent_or_unsupported() -> None:
    fixture = _load_fixture()
    payload = calculate_economics(
        fixture["resultBundle"],
        fixture["economicsConfig"],
    )

    unsupported_components = {
        item["component"]
        for item in payload["coverage"]["unsupportedComponents"]
    }
    unsupported_summary_metrics = {
        row["metric"]
        for row in payload["summaryRows"]
        if row["status"] != "implemented"
    }

    for field in UNSUPPORTED_ECONOMICS_FIELDS:
        assert field not in payload["costs"]
        assert field not in payload["costEffectiveness"]
        assert field not in unsupported_summary_metrics
        assert (
            field in unsupported_components
            or field not in {row["metric"] for row in payload["summaryRows"]}
        )


def test_hand_check_fixture_payload_is_json_serialisable_with_stable_summary_metrics() -> None:
    fixture = _load_fixture()
    payload = calculate_economics(
        fixture["resultBundle"],
        fixture["economicsConfig"],
    )

    json.dumps(payload, allow_nan=False, sort_keys=True)

    assert [row["metric"] for row in payload["summaryRows"]] == [
        "testingCost",
        "treatmentCost",
        "programSetupCost",
        "programRunningCost",
        "falsePositiveIncrementalCost",
        "totalProgramCost",
        "baselineTBDiseaseCost",
        "interventionTBDiseaseCost",
        "tbDiseaseCostsAverted",
        "netCostVsBaseline",
        "costPerInfectionCured",
        "costPerTBCasePrevented",
        "totalImplementedCost",
    ]


def test_python_backend_run_economics_uses_fixture_without_matlab_engine() -> None:
    fixture = _load_fixture()
    sys.modules.pop("matlab", None)
    sys.modules.pop("matlab.engine", None)
    original_import = builtins.__import__

    def fail_on_matlab_import(name, *args, **kwargs):
        if name == "matlab" or name.startswith("matlab."):
            raise AssertionError(f"PythonApyBackend.run_economics imported {name}")
        return original_import(name, *args, **kwargs)

    try:
        builtins.__import__ = fail_on_matlab_import
        payload = PythonApyBackend(repo_root()).run_economics(
            fixture["resultBundle"],
            fixture["economicsConfig"],
        )
    finally:
        builtins.__import__ = original_import

    assert "matlab" not in sys.modules
    assert "matlab.engine" not in sys.modules
    _assert_supported_subset_matches_fixture(payload, fixture)
    json.dumps(payload, allow_nan=False, sort_keys=True)
