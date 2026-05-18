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
    "status",
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
        json.dumps(payload)

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
