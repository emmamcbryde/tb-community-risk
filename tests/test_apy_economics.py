from __future__ import annotations

import copy
import json
import unittest

from engine.apy.economics import calculate_economics


EXPECTED_TOP_LEVEL_KEYS = [
    "available",
    "source",
    "contractVersion",
    "message",
    "metadata",
    "inputs",
    "strategy",
    "quantities",
    "unitCosts",
    "costs",
    "costEffectiveness",
    "summaryRows",
    "coverage",
    "status",
]
EXPECTED_SUMMARY_ROW_KEYS = [
    "metric",
    "label",
    "value",
    "category",
    "status",
    "includedInTotal",
]


def minimal_result_bundle() -> dict:
    return {
        "metadata": {
            "backend": "python",
            "contractVersion": "apy_results_bundle_v9_python_port",
        },
        "results": {
            "interfaceConfig": {
                "testType": "IGRA",
                "regimen": "3HP",
            },
            "summary": [
                {"Metric": "nScreened", "Median": 125.0},
                {"Metric": "nTotalCoursesStarted", "Median": 40.0},
            ],
        },
    }


def minimal_economics_config() -> dict:
    return {
        "metadata": {
            "currencyCode": "AUD",
            "priceYear": 2019,
            "locationLabel": "Australia",
            "sourceNotes": "unit-test fixture",
            "programCostBasis": "unit",
        },
        "costs": {
            "test": {
                "IGRA": 113.48,
                "TST": None,
            },
            "regimen": {
                "x3HP": 165.5072,
                "x4R": None,
                "x3HR": None,
                "x6H": None,
                "x9H": None,
            },
            "falsePositiveIncrementalPerPerson": None,
            "activeTBDiseasePerCase": None,
            "programSetupTotal": None,
            "programRunningTotal": None,
        },
    }


class ApyEconomicsTests(unittest.TestCase):
    def test_minimal_bundle_computes_supported_screening_and_treatment_costs(self) -> None:
        payload = calculate_economics(
            minimal_result_bundle(),
            minimal_economics_config(),
        )

        self.assertEqual(payload["strategy"], {"testType": "IGRA", "regimen": "3HP"})
        self.assertEqual(payload["quantities"]["nScreened"], 125.0)
        self.assertEqual(payload["quantities"]["nTotalCoursesStarted"], 40.0)
        self.assertEqual(payload["unitCosts"]["testPerPerson"], 113.48)
        self.assertEqual(payload["unitCosts"]["treatmentPerStarted"], 165.5072)
        self.assertAlmostEqual(payload["costs"]["testingCost"], 125.0 * 113.48)
        self.assertAlmostEqual(payload["costs"]["treatmentCost"], 40.0 * 165.5072)

    def test_minimal_payload_has_stable_json_contract_and_not_ported_message(self) -> None:
        payload = calculate_economics(
            minimal_result_bundle(),
            minimal_economics_config(),
        )

        self.assertEqual(list(payload.keys()), EXPECTED_TOP_LEVEL_KEYS)
        self.assertIn("full", payload["message"].lower())
        self.assertIn("not ported", payload["message"].lower())
        self.assertEqual(payload["metadata"]["currencyCode"], "AUD")
        self.assertEqual(payload["inputs"], minimal_economics_config())
        for row in payload["summaryRows"]:
            self.assertEqual(list(row.keys()), EXPECTED_SUMMARY_ROW_KEYS)
        json.dumps(payload)

    def test_summary_rows_are_implemented_direct_cost_table(self) -> None:
        payload = calculate_economics(
            minimal_result_bundle(),
            minimal_economics_config(),
        )

        rows = payload["summaryRows"]
        self.assertEqual(
            [row["metric"] for row in rows],
            ["testingCost", "treatmentCost", "totalImplementedCost"],
        )
        for row in rows:
            self.assertEqual(row["status"], "implemented")
            self.assertNotIn(
                row["metric"],
                [
                    item["component"]
                    for item in payload["coverage"]["unsupportedComponents"]
                ],
            )

        component_rows = [row for row in rows if row["includedInTotal"]]
        self.assertEqual(
            [row["metric"] for row in component_rows],
            ["testingCost", "treatmentCost"],
        )
        self.assertTrue(
            all(row["category"] == "directCost" for row in component_rows)
        )

        total_row = rows[-1]
        self.assertEqual(total_row["metric"], "totalImplementedCost")
        self.assertEqual(total_row["category"], "directCostTotal")
        self.assertFalse(total_row["includedInTotal"])
        self.assertAlmostEqual(
            total_row["value"],
            sum(row["value"] for row in component_rows),
        )

    def test_minimal_payload_reports_python_economics_coverage(self) -> None:
        payload = calculate_economics(
            minimal_result_bundle(),
            minimal_economics_config(),
        )

        self.assertEqual(payload["status"], "partial")

        coverage = payload["coverage"]
        self.assertEqual(coverage["status"], "partial")
        self.assertEqual(
            coverage["calculatedComponents"],
            ["testingCost", "treatmentCost"],
        )
        self.assertIn("unsupportedComponents", coverage)
        self.assertIn("notCalculated", coverage)
        self.assertIn("missingInputs", coverage)
        self.assertIn("messages", coverage)
        self.assertEqual(
            coverage["missingInputs"],
            [
                "costs.programSetupTotal",
                "costs.programRunningTotal",
                "summary.nFalsePositiveTreated",
                "costs.falsePositiveIncrementalPerPerson",
            ],
        )
        self.assertIn("falsePositiveIncrementalCost", coverage["notCalculated"])
        self.assertIn("programSetupCost", coverage["notCalculated"])
        self.assertIn("programRunningCost", coverage["notCalculated"])
        self.assertIn("totalProgramCost", coverage["notCalculated"])
        self.assertIn("costPerTBCasePrevented", coverage["notCalculated"])
        self.assertTrue(
            all(
                "not calculated" in item["message"].lower()
                for item in coverage["unsupportedComponents"]
            )
        )
        self.assertTrue(
            any("not ported" in msg.lower() for msg in coverage["messages"])
        )
        json.dumps(coverage)

    def test_direct_program_costs_are_calculated_when_inputs_are_present(self) -> None:
        bundle = minimal_result_bundle()
        bundle["results"]["summary"].append(
            {"Metric": "nFalsePositiveTreated", "Median": 3.0}
        )
        config = minimal_economics_config()
        config["costs"]["falsePositiveIncrementalPerPerson"] = 10.0
        config["costs"]["programSetupTotal"] = 100.0
        config["costs"]["programRunningTotal"] = 5.0

        payload = calculate_economics(bundle, config)

        testing_cost = 125.0 * 113.48
        treatment_cost = 40.0 * 165.5072
        false_positive_cost = 3.0 * 10.0
        expected_total = (
            testing_cost + treatment_cost + false_positive_cost + 100.0 + 5.0
        )
        self.assertEqual(payload["quantities"]["nFalsePositiveTreated"], 3.0)
        self.assertEqual(
            payload["unitCosts"]["falsePositiveIncrementalPerPerson"],
            10.0,
        )
        self.assertEqual(payload["costs"]["falsePositiveIncrementalCost"], 30.0)
        self.assertEqual(payload["costs"]["programSetupCost"], 100.0)
        self.assertEqual(payload["costs"]["programRunningCost"], 5.0)
        self.assertAlmostEqual(payload["costs"]["totalProgramCost"], expected_total)
        self.assertEqual(payload["status"], "partial")
        self.assertEqual(payload["coverage"]["missingInputs"], [])
        for component in [
            "falsePositiveIncrementalCost",
            "programSetupCost",
            "programRunningCost",
            "totalProgramCost",
        ]:
            self.assertIn(component, payload["coverage"]["calculatedComponents"])
            self.assertNotIn(component, payload["coverage"]["notCalculated"])
            self.assertNotIn(
                component,
                [
                    item["component"]
                    for item in payload["coverage"]["unsupportedComponents"]
                ],
            )
        self.assertEqual(
            [row["metric"] for row in payload["summaryRows"]],
            [
                "testingCost",
                "treatmentCost",
                "programSetupCost",
                "programRunningCost",
                "falsePositiveIncrementalCost",
                "totalProgramCost",
                "totalImplementedCost",
            ],
        )
        summary_by_metric = {
            row["metric"]: row
            for row in payload["summaryRows"]
        }
        for component in [
            "testingCost",
            "treatmentCost",
            "programSetupCost",
            "programRunningCost",
            "falsePositiveIncrementalCost",
        ]:
            self.assertEqual(summary_by_metric[component]["category"], "directCost")
            self.assertEqual(summary_by_metric[component]["status"], "implemented")
            self.assertTrue(summary_by_metric[component]["includedInTotal"])

        self.assertFalse(summary_by_metric["totalProgramCost"]["includedInTotal"])
        self.assertEqual(
            summary_by_metric["totalProgramCost"]["category"],
            "directCostTotal",
        )
        self.assertFalse(
            summary_by_metric["totalImplementedCost"]["includedInTotal"]
        )
        self.assertAlmostEqual(
            summary_by_metric["totalImplementedCost"]["value"],
            expected_total,
        )

    def test_false_positive_cost_is_optional_and_missing_metric_does_not_raise(self) -> None:
        config = minimal_economics_config()
        config["costs"]["falsePositiveIncrementalPerPerson"] = 10.0
        config["costs"]["programSetupTotal"] = 100.0
        config["costs"]["programRunningTotal"] = 5.0

        payload = calculate_economics(minimal_result_bundle(), config)

        self.assertNotIn("nFalsePositiveTreated", payload["quantities"])
        self.assertNotIn("falsePositiveIncrementalCost", payload["costs"])
        self.assertNotIn("totalProgramCost", payload["costs"])
        self.assertIn(
            "summary.nFalsePositiveTreated",
            payload["coverage"]["missingInputs"],
        )
        self.assertIn(
            "falsePositiveIncrementalCost",
            payload["coverage"]["notCalculated"],
        )
        self.assertIn("totalProgramCost", payload["coverage"]["notCalculated"])
        self.assertIn("programSetupCost", payload["coverage"]["calculatedComponents"])
        self.assertIn("programRunningCost", payload["coverage"]["calculatedComponents"])
        self.assertNotIn("programSetupCost", payload["coverage"]["notCalculated"])
        self.assertNotIn("programRunningCost", payload["coverage"]["notCalculated"])

    def test_missing_required_summary_metrics_raise_explicit_error(self) -> None:
        bundle = minimal_result_bundle()
        bundle["results"]["summary"] = [
            {"Metric": "nScreened", "Median": 125.0},
        ]

        with self.assertRaisesRegex(ValueError, "nTotalCoursesStarted"):
            calculate_economics(bundle, minimal_economics_config())

    def test_missing_selected_unit_costs_raise_explicit_errors(self) -> None:
        config_without_test_cost = minimal_economics_config()
        config_without_test_cost["costs"]["test"]["IGRA"] = None

        with self.subTest("selected screening test cost"):
            with self.assertRaisesRegex(ValueError, "costs\\.test\\.IGRA"):
                calculate_economics(minimal_result_bundle(), config_without_test_cost)

        config_without_regimen_cost = minimal_economics_config()
        config_without_regimen_cost["costs"]["regimen"]["x3HP"] = None

        with self.subTest("selected preventive-treatment regimen cost"):
            with self.assertRaisesRegex(ValueError, "costs\\.regimen\\.x3HP"):
                calculate_economics(minimal_result_bundle(), config_without_regimen_cost)

    def test_input_types_and_required_sections_are_checked_before_calculation(self) -> None:
        with self.assertRaisesRegex(TypeError, "result_bundle must be a mapping"):
            calculate_economics([], minimal_economics_config())

        with self.assertRaisesRegex(TypeError, "economics_config must be a mapping"):
            calculate_economics(minimal_result_bundle(), [])

        bundle_without_results = copy.deepcopy(minimal_result_bundle())
        bundle_without_results.pop("results")
        with self.assertRaisesRegex(ValueError, "results"):
            calculate_economics(bundle_without_results, minimal_economics_config())

        bundle_with_bad_results = minimal_result_bundle()
        bundle_with_bad_results["results"] = []
        with self.assertRaisesRegex(ValueError, "result_bundle.results must be a mapping"):
            calculate_economics(bundle_with_bad_results, minimal_economics_config())


if __name__ == "__main__":
    unittest.main()
