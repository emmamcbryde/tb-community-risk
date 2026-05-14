from __future__ import annotations

import unittest

from engine.apy.results_bundle import KEY_METRICS, build_results_bundle
from engine.apy.runner import run_replicates


class ApyResultsBundleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = run_replicates(n=100, n_reps=3, seed=20)
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

    def test_dynamic_comparison_is_explicitly_partial(self) -> None:
        dynamic = self.bundle["technical"]["dynamicComparison"]

        self.assertEqual(dynamic["available"], "partial")
        self.assertIn("missingFields", dynamic)
        self.assertIn("do-nothing", dynamic["notes"])
        self.assertIn("cumulative_cases_averted", dynamic)


if __name__ == "__main__":
    unittest.main()
