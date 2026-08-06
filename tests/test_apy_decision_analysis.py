from __future__ import annotations

import json
import unittest

from adapters.serialization import to_json_like
from engine.apy.decision_analysis import (
    DECISION_ANALYSIS_CONTRACT_VERSION,
    build_scenario_config,
    run_scenario_comparison,
)
from engine.apy.sensitivity import (
    SENSITIVITY_CONTRACT_VERSION,
    THRESHOLD_CONTRACT_VERSION,
    crossings_from_evaluated_grid,
    load_sensitivity_specs,
    monotonicity_diagnostic,
    run_one_way_sensitivity,
    run_threshold_analysis,
)
from tests.test_apy_event_ledger_economics import _reviewed_epi, _synthetic_econ


class ApyDecisionAnalysisScenarioComparisonTests(unittest.TestCase):
    def test_expected_value_scenario_comparison_contract_and_repeatability(self) -> None:
        base = _reviewed_epi({"N": 80, "screeningStrategy": "prevent", "screenCoverage": 0.25})
        scenarios = [
            {"scenarioId": "igra_prevent", "label": "IGRA prevent", "changes": {"test": "IGRA"}},
            {"scenarioId": "tst_prevent", "label": "TST prevent", "changes": {"test": "TST"}},
        ]

        first = run_scenario_comparison(base, _synthetic_econ(), scenarios, model_type="expected_value")
        second = run_scenario_comparison(base, _synthetic_econ(), scenarios, model_type="expected_value")

        self.assertEqual(first["contractVersion"], DECISION_ANALYSIS_CONTRACT_VERSION)
        self.assertTrue(first["validation"]["isValid"])
        self.assertEqual(first["scenarios"][0]["valueType"], "expected")
        self.assertEqual(first["scenarioSummaries"][0]["configurationHash"], second["scenarioSummaries"][0]["configurationHash"])
        self.assertEqual(
            first["scenarioSummaries"][0]["active_tb_cases_prevented"],
            second["scenarioSummaries"][0]["active_tb_cases_prevented"],
        )
        json.dumps(to_json_like(first))

    def test_scenario_changes_are_validated_adapters_not_arbitrary_paths(self) -> None:
        base = _reviewed_epi({"N": 20})

        with self.assertRaisesRegex(ValueError, "Unsupported or unvalidated"):
            build_scenario_config(
                base,
                {"scenarioId": "bad", "changes": {"ltbiStateAssumptions.status": "configured_reviewed"}},
            )

        cfg, changed, _ = build_scenario_config(
            base,
            {"scenarioId": "good", "changes": {"testType": "TST", "tstSpecificityBCG": 0.5}},
        )
        self.assertEqual(cfg["testType"], "TST")
        self.assertEqual(cfg["testCharacteristics"]["TST"]["specificityBCG"], 0.5)
        self.assertEqual(changed["tstSpecificityBCG"], 0.5)

    def test_incomplete_economics_does_not_create_strategy_economic_comparison(self) -> None:
        base = _reviewed_epi({"N": 60})
        econ = _synthetic_econ()
        for item in econ["costItems"]:
            if item["costItemId"] == "test_igra":
                item["originalCost"] = None
                item["conversionStatus"] = "unresolved_source_price_year"
                item["convertedTargetYearCost"] = None
        result = run_scenario_comparison(
            base,
            econ,
            [{"scenarioId": "igra", "changes": {"test": "IGRA"}}],
            model_type="expected_value",
        )

        summary = result["scenarioSummaries"][0]
        self.assertIsNone(summary["incrementalCost"])
        self.assertEqual(summary["icerClassification"], "incomplete / not calculated")
        self.assertTrue(result["validation"]["warnings"])

    def test_stochastic_strategy_comparison_preserves_baseline_cohort_fingerprint(self) -> None:
        base = _reviewed_epi({"N": 120, "screenCoverage": 0.25, "screeningStrategy": "prevent"})
        scenarios = [
            {"scenarioId": "igra", "changes": {"test": "IGRA"}},
            {"scenarioId": "tst", "changes": {"test": "TST"}},
        ]

        result = run_scenario_comparison(
            base,
            _synthetic_econ(),
            scenarios,
            model_type="agent_based",
            n_reps=3,
            master_seed=123,
        )

        self.assertTrue(result["validation"]["isValid"])
        self.assertTrue(result["pairedComparisons"][0]["pairedCohortFingerprintMatch"])
        self.assertEqual(result["scenarios"][0]["valueType"], "simulated_count")
        self.assertIsNotNone(result["scenarios"][0]["seed"])


class ApyDecisionAnalysisSensitivityTests(unittest.TestCase):
    def test_default_sensitivity_specs_are_unresolved_without_invented_ranges(self) -> None:
        specs = load_sensitivity_specs()

        self.assertTrue(specs)
        self.assertTrue(all(spec["lowValue"] is None and spec["highValue"] is None for spec in specs))

        result = run_one_way_sensitivity(
            _reviewed_epi({"N": 40}),
            _synthetic_econ(),
            specs[:1],
            ["active_tb_cases_prevented"],
        )

        self.assertEqual(result["contractVersion"], SENSITIVITY_CONTRACT_VERSION)
        self.assertFalse(result["results"][0]["complete"])
        self.assertIn("Reviewed low/high range required", result["results"][0]["unresolvedReason"])

    def test_one_way_sensitivity_runs_only_explicit_supplied_ranges(self) -> None:
        spec = {
            "parameterId": "epi.tpt_initiation.synthetic",
            "label": "Synthetic TPT initiation",
            "adapter": "pStartTPT",
            "baseValue": 0.85,
            "lowValue": 0.0,
            "highValue": 1.0,
            "unit": "probability",
            "scale": "probability",
            "source": "Synthetic unit-test range",
            "reviewStatus": "configured_reviewed",
            "provisional": False,
        }

        result = run_one_way_sensitivity(
            _reviewed_epi({"N": 60, "screeningStrategy": "random"}),
            _synthetic_econ(),
            [spec],
            ["tpt_started_total", "active_tb_cases_prevented", "incrementalCost"],
        )

        rows = result["results"]
        self.assertTrue(all(row["complete"] for row in rows))
        starts = next(row for row in rows if row["outcome"] == "tpt_started_total")
        self.assertEqual(starts["lowOutcome"], 0)
        self.assertGreater(starts["highOutcome"], starts["baseOutcome"])
        self.assertTrue(result["tornadoRows"])

    def test_sensitivity_cost_adapter_updates_authoritative_cost_item(self) -> None:
        spec = {
            "parameterId": "cost.program_running.synthetic",
            "label": "Synthetic running cost",
            "adapter": "programRunningCost",
            "baseValue": 0,
            "lowValue": 0,
            "highValue": 100,
            "unit": "AUD",
            "scale": "monetary",
            "source": "Synthetic unit-test range",
            "reviewStatus": "configured_reviewed",
            "provisional": False,
        }

        result = run_one_way_sensitivity(
            _reviewed_epi({"N": 50, "screeningWindowYears": 2}),
            _synthetic_econ({"program_running": 0}),
            [spec],
            ["incrementalCost"],
        )

        row = result["results"][0]
        self.assertTrue(row["complete"])
        self.assertGreater(row["highOutcome"], row["lowOutcome"])

    def test_threshold_analysis_preserves_grid_and_requires_valid_threshold_for_nmb(self) -> None:
        spec = {
            "parameterId": "epi.tpt_initiation.synthetic",
            "adapter": "pStartTPT",
            "label": "Synthetic TPT initiation",
        }

        unavailable = run_threshold_analysis(
            _reviewed_epi({"N": 40}),
            _synthetic_econ({"threshold": None}),
            spec,
            "nmb",
            {"low": 0, "high": 1, "gridPoints": 5},
        )
        self.assertFalse(unavailable["validation"]["isValid"])

        result = run_threshold_analysis(
            _reviewed_epi({"N": 20}),
            _synthetic_econ({"threshold": 50000}),
            spec,
            "active_tb_cases_prevented",
            {"low": 0, "high": 1, "gridPoints": 3, "target": -999},
        )

        self.assertEqual(result["contractVersion"], THRESHOLD_CONTRACT_VERSION)
        self.assertEqual(len(result["grid"]), 3)
        self.assertEqual(result["crossingCount"], 0)
        self.assertIn("isMonotonic", result["monotonicity"])

    def test_nonmonotonic_and_multiple_crossing_grid_is_detected(self) -> None:
        grid = [
            {"parameterValue": 0, "difference": -1, "complete": True},
            {"parameterValue": 1, "difference": 1, "complete": True},
            {"parameterValue": 2, "difference": -1, "complete": True},
            {"parameterValue": 3, "difference": 1, "complete": True},
        ]

        self.assertFalse(monotonicity_diagnostic(grid)["isMonotonic"])
        crossings = crossings_from_evaluated_grid(grid)
        self.assertEqual(len(crossings), 3)
        self.assertTrue(all(crossing["bracketed"] for crossing in crossings))


if __name__ == "__main__":
    unittest.main()
