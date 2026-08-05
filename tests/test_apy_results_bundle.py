from __future__ import annotations

import unittest

from engine.apy.results_bundle import KEY_METRICS, build_results_bundle
from engine.apy.natural_history import run_do_nothing_summary
from engine.apy.ltbi_state import enable_development_compatibility_mode
from engine.apy.runner import run_replicates


class ApyResultsBundleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = run_replicates(
            enable_development_compatibility_mode({}), n=100, n_reps=3, seed=20
        )
        cls.bundle = build_results_bundle(cls.results)

    def test_bundle_has_core_sections(self) -> None:
        self.assertIn("metadata", self.bundle)
        self.assertIn("headline", self.bundle)
        self.assertIn("technical", self.bundle)
        self.assertIn("downloads", self.bundle)

    def test_metadata_identifies_python_backend(self) -> None:
        self.assertEqual(self.bundle["metadata"]["modelVersion"], "python_apy_v9_port")
        self.assertEqual(self.bundle["metadata"]["backend"], "python")

    def test_key_metrics_rows_include_expected_metrics(self) -> None:
        metrics = {
            row["Metric"] for row in self.bundle["headline"]["keyMetricsRows"]
        }

        for metric in KEY_METRICS:
            self.assertIn(metric, metrics)

    def test_interface_config_and_calibration_are_present(self) -> None:
        self.assertIn("interfaceConfig", self.bundle["technical"])
        self.assertIn("calibration", self.bundle["technical"])
        self.assertIn("ageInfGamma", self.bundle["technical"]["calibration"])

    def test_dynamic_comparison_uses_event_ledger_when_available(self) -> None:
        dynamic = self.bundle["technical"]["dynamicComparison"]

        self.assertIs(dynamic["available"], True)
        self.assertEqual(dynamic["source"], "technical.eventLedger")
        self.assertIn("cumulative_cases_averted", dynamic)
        self.assertAlmostEqual(
            dynamic["cumulative_baseline_active_tb_cases"],
            dynamic["cumulative_intervention_active_tb_cases"]
            + dynamic["cumulative_cases_averted"],
        )

    def test_bundle_with_do_nothing_has_complete_dynamic_comparison(self) -> None:
        do_nothing = run_do_nothing_summary(self.results)
        bundle = build_results_bundle(self.results, do_nothing=do_nothing)
        dynamic = bundle["technical"]["dynamicComparison"]

        self.assertIs(dynamic["available"], True)
        self.assertEqual(dynamic["source"], "doNothing.derived")
        self.assertEqual(dynamic["missingFields"], [])

    def test_complete_dynamic_comparison_metric_rows_include_required_metrics(self) -> None:
        do_nothing = run_do_nothing_summary(self.results)
        bundle = build_results_bundle(self.results, do_nothing=do_nothing)
        rows = bundle["technical"]["dynamicComparison"]["metricRows"]
        metrics = {row["Metric"] for row in rows}

        self.assertEqual(
            metrics,
            {
                "cumulative_baseline_active_tb_cases",
                "cumulative_intervention_active_tb_cases",
                "cumulative_cases_averted",
                "relative_reduction_cumulative_active_tb_cases",
            },
        )


if __name__ == "__main__":
    unittest.main()
