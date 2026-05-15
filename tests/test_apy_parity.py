from __future__ import annotations

import unittest
from pathlib import Path

import pandas as pd

from engine.apy.parity import (
    compare_dynamic_comparison_to_reference,
    compare_python_summary_to_reference,
)
from engine.apy.reference_loader import load_reference_dir
from engine.apy.runner import run_replicates, run_scenario_with_do_nothing


FIXTURE_DIR = (
    Path("validation") / "matlab_reference" / "default_random_igra_3hp_N1500_seed1"
)


class ApyParityTests(unittest.TestCase):
    def test_compact_matlab_reference_fixture_loads(self) -> None:
        loaded = load_reference_dir(FIXTURE_DIR)

        self.assertIn("summary", loaded)
        self.assertFalse(loaded["summary"].empty)

    def test_python_summary_compares_to_reference_without_crashing(self) -> None:
        reference = load_reference_dir(FIXTURE_DIR)
        python_results = run_replicates(n=100, n_reps=2, seed=1)

        comparison = compare_python_summary_to_reference(
            python_results["summary"],
            reference["summary"],
            metrics=["nScreened", "nPreventedActiveTB", "missing_metric"],
        )

        self.assertEqual(
            list(comparison.columns),
            [
                "Metric",
                "PythonMedian",
                "MatlabMedian",
                "AbsoluteDifference",
                "RelativeDifference",
                "Comparable",
                "Notes",
            ],
        )
        self.assertEqual(len(comparison), 3)
        self.assertFalse(
            bool(comparison.loc[comparison["Metric"] == "missing_metric", "Comparable"].iloc[0])
        )

    def test_missing_metrics_are_handled_gracefully(self) -> None:
        python_summary = pd.DataFrame(
            [{"Metric": "a", "Median": 1, "Low95": 1, "High95": 1}]
        )
        reference_summary = pd.DataFrame(
            [{"Metric": "b", "Median": 2, "Low95": 2, "High95": 2}]
        )

        comparison = compare_python_summary_to_reference(
            python_summary, reference_summary, metrics=["a", "b"]
        )

        self.assertFalse(comparison["Comparable"].any())

    def test_dynamic_comparison_compares_to_reference_without_crashing(self) -> None:
        reference = load_reference_dir(FIXTURE_DIR)
        python_out = run_scenario_with_do_nothing({"N": 100, "nReps": 2, "seed": 13})

        comparison = compare_dynamic_comparison_to_reference(
            python_out["bundle"]["technical"]["dynamicComparison"],
            reference["dynamic_comparison"],
        )

        self.assertEqual(
            list(comparison.columns),
            [
                "Metric",
                "PythonMedian",
                "MatlabMedian",
                "AbsoluteDifference",
                "RelativeDifference",
                "Comparable",
                "Notes",
            ],
        )
        self.assertEqual(len(comparison), 4)
        self.assertTrue(comparison["Notes"].str.contains("strict equality").any())


if __name__ == "__main__":
    unittest.main()
