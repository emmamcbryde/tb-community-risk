from __future__ import annotations

import json
import unittest

from adapters.serialization import to_json_like
from engine.apy.decision_analysis import (
    DECISION_ANALYSIS_CONTRACT_VERSION,
    build_scenario_config,
    run_scenario_comparison,
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


if __name__ == "__main__":
    unittest.main()
