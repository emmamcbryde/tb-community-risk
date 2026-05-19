from __future__ import annotations

import copy
import json
import unittest
from pathlib import Path

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
EXPECTED_SOURCE = "apy_python_partial_economics"
EXPECTED_CONTRACT_VERSION = "apy_economics_partial_v1"
EXPECTED_LEGACY_SOURCE = "apy_python_minimal_economics"
EXPECTED_LEGACY_CONTRACT_VERSION = "apy_economics_minimal_v1"
FIXTURE_DIR = Path(__file__).parent / "fixtures"


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


def _load_hand_check_fixture() -> dict:
    fixture_path = FIXTURE_DIR / "apy_economics_hand_check_fixture.json"
    with fixture_path.open(encoding="utf-8") as fixture_file:
        return json.load(fixture_file)


def _supported_economics_subset(payload: dict) -> dict:
    supported_costs = [
        "testingCost",
        "treatmentCost",
        "falsePositiveIncrementalCost",
        "programSetupCost",
        "programRunningCost",
        "totalProgramCost",
        "baselineTBDiseaseCost",
        "interventionTBDiseaseCost",
        "tbDiseaseCostsAverted",
        "netCostVsBaseline",
    ]
    supported_ratios = [
        "costPerInfectionCured",
        "costPerTBCasePrevented",
    ]
    return {
        "costs": {metric: payload["costs"][metric] for metric in supported_costs},
        "costEffectiveness": {
            metric: payload["costEffectiveness"][metric]
            for metric in supported_ratios
        },
    }


def test_hand_checkable_python_economics_fixture_matches_supported_subset() -> None:
    fixture = _load_hand_check_fixture()
    payload = calculate_economics(
        fixture["resultBundle"],
        fixture["economicsConfig"],
    )

    assert fixture["source"] == "Python unit-test fixture, not MATLAB-exported"
    assert _supported_economics_subset(payload) == fixture["expectedSupportedSubset"]
    assert "results.dynamicComparison" not in payload["coverage"]["missingInputs"]


class ApyEconomicsTests(unittest.TestCase):
    def test_partial_bundle_computes_supported_screening_and_treatment_costs(self) -> None:
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

    def test_partial_payload_has_stable_json_contract_and_not_ported_message(self) -> None:
        payload = calculate_economics(
            minimal_result_bundle(),
            minimal_economics_config(),
        )

        self.assertEqual(list(payload.keys()), EXPECTED_TOP_LEVEL_KEYS)
        self.assertEqual(payload["source"], EXPECTED_SOURCE)
        self.assertEqual(payload["contractVersion"], EXPECTED_CONTRACT_VERSION)
        self.assertNotIn("legacySource", payload)
        self.assertNotIn("legacyContractVersion", payload)
        self.assertIn("partial", payload["message"].lower())
        self.assertNotIn("minimal", payload["message"].lower())
        self.assertIn("full", payload["message"].lower())
        self.assertIn("not ported", payload["message"].lower())
        self.assertEqual(payload["metadata"]["currencyCode"], "AUD")
        self.assertEqual(
            payload["metadata"]["legacyIdentifiers"],
            {
                "source": EXPECTED_LEGACY_SOURCE,
                "contractVersion": EXPECTED_LEGACY_CONTRACT_VERSION,
            },
        )
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

    def test_partial_payload_reports_python_economics_coverage(self) -> None:
        payload = calculate_economics(
            minimal_result_bundle(),
            minimal_economics_config(),
        )

        self.assertEqual(payload["status"], "partial")
        self.assertEqual(payload["source"], EXPECTED_SOURCE)
        self.assertEqual(payload["contractVersion"], EXPECTED_CONTRACT_VERSION)

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
                "results.dynamicComparison",
                "costs.activeTBDiseasePerCase",
                "summary.nPreventedActiveTB",
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
        self.assertTrue(
            all("minimal" not in msg.lower() for msg in coverage["messages"])
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
        self.assertNotIn(
            "summary.nFalsePositiveTreated",
            payload["coverage"]["missingInputs"],
        )
        self.assertNotIn(
            "costs.falsePositiveIncrementalPerPerson",
            payload["coverage"]["missingInputs"],
        )
        self.assertNotIn(
            "costs.programSetupTotal",
            payload["coverage"]["missingInputs"],
        )
        self.assertNotIn(
            "costs.programRunningTotal",
            payload["coverage"]["missingInputs"],
        )
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

    def test_tb_disease_benefit_costs_are_calculated_when_explicit_inputs_are_present(self) -> None:
        bundle = minimal_result_bundle()
        bundle["results"]["summary"].extend(
            [
                {"Metric": "nFalsePositiveTreated", "Median": 3.0},
                {"Metric": "nPreventedActiveTB", "Median": 2.0},
                {"Metric": "nCuredInfection", "Median": 8.0},
            ]
        )
        bundle["results"]["dynamicComparison"] = {
            "cumulative_baseline_active_tb_cases": 10.0,
            "cumulative_intervention_active_tb_cases": 8.0,
        }
        config = minimal_economics_config()
        config["costs"]["falsePositiveIncrementalPerPerson"] = 10.0
        config["costs"]["programSetupTotal"] = 100.0
        config["costs"]["programRunningTotal"] = 5.0
        config["costs"]["activeTBDiseasePerCase"] = 1000.0

        payload = calculate_economics(bundle, config)

        total_program_cost = (
            125.0 * 113.48
            + 40.0 * 165.5072
            + 3.0 * 10.0
            + 100.0
            + 5.0
        )
        expected_net_cost = total_program_cost - 2.0 * 1000.0
        self.assertEqual(payload["status"], "partial")
        self.assertEqual(payload["unitCosts"]["activeTBDiseasePerCase"], 1000.0)
        self.assertEqual(payload["quantities"]["nPreventedActiveTB"], 2.0)
        self.assertEqual(payload["quantities"]["nCuredInfection"], 8.0)
        self.assertEqual(payload["quantities"]["cumulativeBaselineActiveTBCases"], 10.0)
        self.assertEqual(
            payload["quantities"]["cumulativeInterventionActiveTBCases"],
            8.0,
        )
        self.assertEqual(payload["costs"]["baselineTBDiseaseCost"], 10000.0)
        self.assertEqual(payload["costs"]["interventionTBDiseaseCost"], 8000.0)
        self.assertEqual(payload["costs"]["tbDiseaseCostsAverted"], 2000.0)
        self.assertAlmostEqual(payload["costs"]["netCostVsBaseline"], expected_net_cost)
        self.assertAlmostEqual(
            payload["costEffectiveness"]["costPerInfectionCured"],
            expected_net_cost / 8.0,
        )
        self.assertAlmostEqual(
            payload["costEffectiveness"]["costPerTBCasePrevented"],
            expected_net_cost / 2.0,
        )

        for component in [
            "baselineTBDiseaseCost",
            "interventionTBDiseaseCost",
            "tbDiseaseCostsAverted",
            "netCostVsBaseline",
            "costPerInfectionCured",
            "costPerTBCasePrevented",
        ]:
            self.assertIn(component, payload["coverage"]["calculatedComponents"])
            self.assertNotIn(component, payload["coverage"]["notCalculated"])
        self.assertNotIn("results.dynamicComparison", payload["coverage"]["missingInputs"])
        self.assertEqual(payload["coverage"]["unsupportedComponents"], [])

        summary_by_metric = {
            row["metric"]: row
            for row in payload["summaryRows"]
        }
        for component in [
            "baselineTBDiseaseCost",
            "interventionTBDiseaseCost",
            "tbDiseaseCostsAverted",
        ]:
            self.assertEqual(summary_by_metric[component]["category"], "benefitCost")
            self.assertFalse(summary_by_metric[component]["includedInTotal"])
        self.assertEqual(summary_by_metric["netCostVsBaseline"]["category"], "netCost")
        self.assertEqual(
            summary_by_metric["costPerTBCasePrevented"]["category"],
            "costEffectiveness",
        )
        json.dumps(payload)

    def test_missing_dynamic_comparison_leaves_tb_case_costs_not_calculated(self) -> None:
        bundle = minimal_result_bundle()
        bundle["results"]["summary"].append(
            {"Metric": "nPreventedActiveTB", "Median": 2.0}
        )
        config = minimal_economics_config()
        config["costs"]["activeTBDiseasePerCase"] = 1000.0

        payload = calculate_economics(bundle, config)

        self.assertNotIn("baselineTBDiseaseCost", payload["costs"])
        self.assertNotIn("interventionTBDiseaseCost", payload["costs"])
        self.assertNotIn(
            "baselineTBDiseaseCost",
            [row["metric"] for row in payload["summaryRows"]],
        )
        self.assertNotIn(
            "interventionTBDiseaseCost",
            [row["metric"] for row in payload["summaryRows"]],
        )
        self.assertEqual(payload["costs"]["tbDiseaseCostsAverted"], 2000.0)
        self.assertEqual(payload["costEffectiveness"], {})
        self.assertIn("results.dynamicComparison", payload["coverage"]["missingInputs"])
        self.assertIn("baselineTBDiseaseCost", payload["coverage"]["notCalculated"])
        self.assertIn("interventionTBDiseaseCost", payload["coverage"]["notCalculated"])
        self.assertIn("netCostVsBaseline", payload["coverage"]["notCalculated"])
        self.assertIn("costPerInfectionCured", payload["coverage"]["notCalculated"])
        self.assertIn("costPerTBCasePrevented", payload["coverage"]["notCalculated"])
        self.assertTrue(
            any(
                "dynamiccomparison is missing" in message.lower()
                for message in payload["coverage"]["messages"]
            )
        )
        json.dumps(payload)

    def test_absent_ratio_denominator_is_explicitly_missing_and_not_calculated(self) -> None:
        bundle = minimal_result_bundle()
        bundle["results"]["summary"].extend(
            [
                {"Metric": "nFalsePositiveTreated", "Median": 3.0},
                {"Metric": "nPreventedActiveTB", "Median": 2.0},
            ]
        )
        bundle["results"]["dynamicComparison"] = {
            "cumulative_baseline_active_tb_cases": 10.0,
            "cumulative_intervention_active_tb_cases": 8.0,
        }
        config = minimal_economics_config()
        config["costs"]["falsePositiveIncrementalPerPerson"] = 10.0
        config["costs"]["programSetupTotal"] = 100.0
        config["costs"]["programRunningTotal"] = 5.0
        config["costs"]["activeTBDiseasePerCase"] = 1000.0

        payload = calculate_economics(bundle, config)

        self.assertIn("netCostVsBaseline", payload["costs"])
        self.assertIn("costPerTBCasePrevented", payload["costEffectiveness"])
        self.assertNotIn("costPerInfectionCured", payload["costEffectiveness"])
        self.assertIn("summary.nCuredInfection", payload["coverage"]["missingInputs"])
        self.assertIn("costPerInfectionCured", payload["coverage"]["notCalculated"])
        self.assertNotIn(
            "costPerInfectionCured",
            [row["metric"] for row in payload["summaryRows"]],
        )
        self.assertTrue(
            any(
                "summary.ncuredinfection is missing" in message.lower()
                for message in payload["coverage"]["messages"]
            )
        )
        json.dumps(payload)

    def test_simple_ratios_require_net_cost_and_positive_denominators(self) -> None:
        bundle = minimal_result_bundle()
        bundle["results"]["summary"].extend(
            [
                {"Metric": "nFalsePositiveTreated", "Median": 3.0},
                {"Metric": "nPreventedActiveTB", "Median": 0.0},
                {"Metric": "nCuredInfection", "Median": 0.0},
            ]
        )
        bundle["results"]["dynamicComparison"] = {
            "cumulative_baseline_active_tb_cases": 10.0,
            "cumulative_intervention_active_tb_cases": 10.0,
        }
        config = minimal_economics_config()
        config["costs"]["falsePositiveIncrementalPerPerson"] = 10.0
        config["costs"]["programSetupTotal"] = 100.0
        config["costs"]["programRunningTotal"] = 5.0
        config["costs"]["activeTBDiseasePerCase"] = 1000.0

        payload = calculate_economics(bundle, config)

        self.assertIn("netCostVsBaseline", payload["costs"])
        self.assertEqual(payload["costEffectiveness"], {})
        self.assertIn("costPerInfectionCured", payload["coverage"]["notCalculated"])
        self.assertIn("costPerTBCasePrevented", payload["coverage"]["notCalculated"])
        self.assertNotIn(
            "costPerInfectionCured",
            [row["metric"] for row in payload["summaryRows"]],
        )
        self.assertNotIn(
            "costPerTBCasePrevented",
            [row["metric"] for row in payload["summaryRows"]],
        )
        self.assertTrue(
            any(
                "summary.ncuredinfection is not positive" in message.lower()
                for message in payload["coverage"]["messages"]
            )
        )
        self.assertTrue(
            any(
                "summary.npreventedactivetb is not positive" in message.lower()
                for message in payload["coverage"]["messages"]
            )
        )
        self.assertNotIn("summary.nCuredInfection", payload["coverage"]["missingInputs"])
        self.assertNotIn("summary.nPreventedActiveTB", payload["coverage"]["missingInputs"])
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
