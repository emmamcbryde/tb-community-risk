from __future__ import annotations

import unittest

from engine.apy.economics_compare import compare_economics_rows, economics_rows


class ApyEconomicsCompareTests(unittest.TestCase):
    def test_economics_rows_extracts_summary_rows(self) -> None:
        summary_rows = [{"metric": "testingCost", "value": 100.0}]

        self.assertEqual(economics_rows({"summaryRows": summary_rows}), summary_rows)
        self.assertEqual(economics_rows(None), [])
        self.assertEqual(economics_rows({"summaryRows": None}), [])

    def test_matches_rows_by_metric_and_component_keys(self) -> None:
        baseline_rows = [
            {
                "metric": "programSetupCost",
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

        rows, warnings = compare_economics_rows(baseline_rows, comparator_rows)

        self.assertEqual(warnings, [])
        self.assertEqual(
            rows,
            [
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
            ],
        )

    def test_calculates_numeric_baseline_comparator_and_differences(self) -> None:
        baseline_rows = [
            {"Metric": "testingCost", "Value": "100"},
            {"metric": "treatmentCost", "value": 0},
        ]
        comparator_rows = [
            {"Metric": "testingCost", "Value": 130.0},
            {"metric": "treatmentCost", "value": 25.0},
        ]

        rows, warnings = compare_economics_rows(baseline_rows, comparator_rows)

        self.assertEqual(warnings, [])
        by_metric = {row["metric"]: row for row in rows}
        self.assertEqual(by_metric["testingCost"]["baseline"], "100")
        self.assertEqual(by_metric["testingCost"]["comparator"], 130.0)
        self.assertEqual(by_metric["testingCost"]["absoluteDifference"], 30.0)
        self.assertEqual(by_metric["testingCost"]["relativeDifference"], 0.3)
        self.assertEqual(by_metric["treatmentCost"]["absoluteDifference"], 25.0)
        self.assertIsNone(by_metric["treatmentCost"]["relativeDifference"])

    def test_retains_non_numeric_keyed_rows_with_missing_numeric_differences(self) -> None:
        baseline_rows = [
            {
                "metric": "costPerTBCasePrevented",
                "value": "not calculated",
                "status": "unsupported",
            },
            {
                "component": "totalProgramCost",
                "status": "missing",
            },
        ]
        comparator_rows = [
            {
                "metric": "costPerTBCasePrevented",
                "value": "not calculated",
                "status": "unsupported",
            },
            {
                "component": "totalProgramCost",
                "value": 400.0,
                "status": "implemented",
            },
        ]

        rows, warnings = compare_economics_rows(baseline_rows, comparator_rows)

        self.assertEqual(warnings, [])
        by_metric = {row["metric"]: row for row in rows}
        self.assertEqual(
            by_metric["costPerTBCasePrevented"],
            {
                "metric": "costPerTBCasePrevented",
                "baseline": "not calculated",
                "comparator": "not calculated",
                "absoluteDifference": None,
                "relativeDifference": None,
                "status": "unsupported",
            },
        )
        self.assertEqual(
            by_metric["totalProgramCost"],
            {
                "metric": "totalProgramCost",
                "baseline": None,
                "comparator": 400.0,
                "absoluteDifference": None,
                "relativeDifference": None,
                "status": "missing",
            },
        )

    def test_preserves_unsupported_and_missing_component_status(self) -> None:
        baseline_rows = [
            {
                "component": "tbDiseaseCostsAverted",
                "category": "diseaseCost",
                "status": "unsupported",
                "includedInTotal": False,
            },
            {
                "component": "falsePositiveIncrementalCost",
                "category": "directCost",
                "status": "missing",
                "includedInTotal": True,
            },
        ]
        comparator_rows = [
            {
                "component": "tbDiseaseCostsAverted",
                "category": "diseaseCost",
                "status": "unsupported",
                "includedInTotal": False,
            },
            {
                "component": "falsePositiveIncrementalCost",
                "category": "directCost",
                "status": "missing",
                "includedInTotal": True,
            },
        ]

        rows, warnings = compare_economics_rows(baseline_rows, comparator_rows)

        self.assertEqual(warnings, [])
        by_metric = {row["metric"]: row for row in rows}
        self.assertEqual(by_metric["tbDiseaseCostsAverted"]["status"], "unsupported")
        self.assertEqual(by_metric["tbDiseaseCostsAverted"]["category"], "diseaseCost")
        self.assertFalse(by_metric["tbDiseaseCostsAverted"]["includedInTotal"])
        self.assertEqual(by_metric["falsePositiveIncrementalCost"]["status"], "missing")
        self.assertEqual(by_metric["falsePositiveIncrementalCost"]["category"], "directCost")
        self.assertTrue(by_metric["falsePositiveIncrementalCost"]["includedInTotal"])


if __name__ == "__main__":
    unittest.main()
