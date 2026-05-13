from __future__ import annotations

import unittest

from engine.integration.compare_dynamic_abm_v9 import compare_dynamic_abm_v9


def dynamic_bundle() -> dict:
    return {
        "modelVersion": "dynamic_python_v1",
        "metadata": {"population": 1000},
        "headline": {
            "keyMetricsRows": [
                {"Metric": "projection_horizon", "Value": 20},
                {"Metric": "cumulative_baseline_active_tb_cases", "Value": 100},
                {"Metric": "cumulative_intervention_active_tb_cases", "Value": 75},
                {"Metric": "cumulative_cases_averted", "Value": 25},
                {"Metric": "relative_reduction_cumulative_active_tb_cases", "Value": 0.25},
            ]
        },
    }


def abm_bundle() -> dict:
    return {
        "metadata": {"modelVersion": "apy_v9"},
        "technical": {"interfaceConfig": {"N": 1000, "followHorizon": 20}},
        "headline": {
            "keyMetricsRows": [
                {"Metric": "cumulative_baseline_active_tb_cases", "Median": 90},
                {"Metric": "cumulative_intervention_active_tb_cases", "Median": 70},
                {"Metric": "nPreventedActiveTB", "Median": 20},
                {"Metric": "relative_reduction", "Median": 0.2222},
            ]
        },
    }


def abm_bundle_with_dynamic_comparison() -> dict:
    return {
        "metadata": {"modelVersion": "apy_v9"},
        "technical": {
            "interfaceConfig": {"N": 1000, "followHorizon": 20},
            "dynamicComparison": {
                "available": True,
                "population": 1000,
                "followHorizon": 20,
                "cumulative_baseline_active_tb_cases": 90,
                "cumulative_intervention_active_tb_cases": 70,
                "cumulative_cases_averted": 20,
                "relative_reduction_cumulative_active_tb_cases": 0.2222,
                "metricRows": [
                    {"Metric": "cumulative_baseline_active_tb_cases", "Median": 90},
                    {"Metric": "cumulative_intervention_active_tb_cases", "Median": 70},
                    {"Metric": "cumulative_cases_averted", "Median": 20},
                    {"Metric": "relative_reduction_cumulative_active_tb_cases", "Median": 0.2222},
                ],
            },
        },
        "headline": {"keyMetricsRows": [{"Metric": "nPreventedActiveTB", "Median": 20}]},
    }


def abm_bundle_cases_averted_only() -> dict:
    return {
        "metadata": {"modelVersion": "apy_v9"},
        "technical": {"interfaceConfig": {"N": 1000, "followHorizon": 20}},
        "headline": {
            "keyMetricsRows": [
                {"Metric": "nPreventedActiveTB", "Median": 20},
            ]
        },
    }


class CompareDynamicAbmV9Tests(unittest.TestCase):
    def test_comparison_aligns_available_metrics(self) -> None:
        comparison = compare_dynamic_abm_v9(dynamic_bundle(), abm_bundle())
        rows = comparison["comparisonRows"]
        by_metric = {row["metric"]: row for row in rows}

        self.assertTrue(by_metric["cumulative cases averted"]["comparable"])
        self.assertEqual(by_metric["cumulative cases averted"]["dynamic_value"], 25)
        self.assertEqual(by_metric["cumulative cases averted"]["abm_value"], 20)
        self.assertEqual(by_metric["cumulative cases averted"]["absolute_difference"], 5)

    def test_missing_metrics_are_handled_gracefully(self) -> None:
        comparison = compare_dynamic_abm_v9({"headline": {"keyMetricsRows": []}}, abm_bundle())

        self.assertTrue(comparison["comparisonRows"])
        self.assertTrue(comparison["warnings"])
        self.assertTrue(comparison["missing_dynamic_metrics"])
        self.assertTrue(any(not row["comparable"] for row in comparison["comparisonRows"]))

    def test_non_comparable_metrics_are_warned_not_silently_dropped(self) -> None:
        comparison = compare_dynamic_abm_v9(dynamic_bundle(), abm_bundle())
        warning_text = "\n".join(comparison["warnings"])

        self.assertIn("NNS/NNT", warning_text)
        self.assertIn("False positives treated", warning_text)
        self.assertIn("ABM economics outputs", warning_text)
        self.assertTrue(comparison["structurally_non_comparable_metrics"])

    def test_dynamic_comparison_block_produces_six_comparable_rows(self) -> None:
        comparison = compare_dynamic_abm_v9(dynamic_bundle(), abm_bundle_with_dynamic_comparison())

        self.assertEqual(sum(1 for row in comparison["comparisonRows"] if row["comparable"]), 6)
        self.assertEqual(comparison["missing_abm_metrics"], [])

    def test_bundle_without_baseline_intervention_counts_warns(self) -> None:
        comparison = compare_dynamic_abm_v9(dynamic_bundle(), abm_bundle_cases_averted_only())

        self.assertTrue(comparison["missing_abm_metrics"])
        self.assertTrue(any(not row["comparable"] for row in comparison["comparisonRows"]))

    def test_cases_averted_fallback_keeps_three_comparable_metrics(self) -> None:
        comparison = compare_dynamic_abm_v9(dynamic_bundle(), abm_bundle_cases_averted_only())
        comparable = [row["metric"] for row in comparison["comparisonRows"] if row["comparable"]]

        self.assertEqual(comparable, ["horizon", "population", "cumulative cases averted"])

    def test_horizon_mismatch_warns(self) -> None:
        abm = abm_bundle_with_dynamic_comparison()
        abm["technical"]["interfaceConfig"]["followHorizon"] = 10
        comparison = compare_dynamic_abm_v9(dynamic_bundle(), abm)

        self.assertTrue(any("horizon" in warning for warning in comparison["warnings"]))

    def test_population_mismatch_warns(self) -> None:
        abm = abm_bundle_with_dynamic_comparison()
        abm["technical"]["interfaceConfig"]["N"] = 2000
        comparison = compare_dynamic_abm_v9(dynamic_bundle(), abm)

        self.assertTrue(any("population" in warning for warning in comparison["warnings"]))
        self.assertTrue(any("different population sizes" in warning for warning in comparison["warnings"]))

    def test_relative_reduction_requires_both_models(self) -> None:
        dynamic = dynamic_bundle()
        dynamic["headline"]["keyMetricsRows"] = [
            row
            for row in dynamic["headline"]["keyMetricsRows"]
            if row["Metric"] != "relative_reduction_cumulative_active_tb_cases"
        ]
        comparison = compare_dynamic_abm_v9(dynamic, abm_bundle_with_dynamic_comparison())
        by_metric = {row["metric"]: row for row in comparison["comparisonRows"]}

        self.assertFalse(by_metric["relative reduction"]["comparable"])
        self.assertIn("relative reduction (relative_reduction_cumulative_active_tb_cases)", comparison["missing_dynamic_metrics"])


if __name__ == "__main__":
    unittest.main()
