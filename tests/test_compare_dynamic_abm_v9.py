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
        self.assertTrue(any(not row["comparable"] for row in comparison["comparisonRows"]))

    def test_non_comparable_metrics_are_warned_not_silently_dropped(self) -> None:
        comparison = compare_dynamic_abm_v9(dynamic_bundle(), abm_bundle())
        warning_text = "\n".join(comparison["warnings"])

        self.assertIn("NNS/NNT", warning_text)
        self.assertIn("False positives treated", warning_text)
        self.assertIn("ABM economics outputs", warning_text)


if __name__ == "__main__":
    unittest.main()
