from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from engine.apy.distributional_validation import (
    build_distributional_validation_table,
    portable_config_from_reference,
    run_reference_distributional_validation,
)
from engine.apy.reference_loader import load_reference_dir


FIXTURE_DIR = (
    Path("validation") / "matlab_reference" / "default_random_igra_3hp_N1500_seed1"
)


class ApyDistributionalValidationTests(unittest.TestCase):
    def test_portable_config_preserves_reference_controls(self) -> None:
        reference = load_reference_dir(FIXTURE_DIR)
        config = portable_config_from_reference(reference["scenario_config"])

        self.assertEqual(config["N"], 1500)
        self.assertEqual(config["nReps"], 2000)
        self.assertEqual(config["seed"], 1)
        self.assertEqual(config["followHorizon"], 20)
        self.assertEqual(config["screenWindow"], 2)
        self.assertEqual(config["screeningStrategy"], "random")
        self.assertEqual(config["testType"], "IGRA")
        self.assertEqual(config["regimen"], "3HP")

    def test_validation_table_has_required_columns(self) -> None:
        python_summary = pd.DataFrame(
            [{"Metric": "nScreened", "Median": 10, "Low95": 10, "High95": 10}]
        )
        matlab_summary = pd.DataFrame(
            [{"Metric": "nScreened", "Median": 12, "Low95": 12, "High95": 12}]
        )
        python_dynamic = {"cumulative_cases_averted": 2}
        matlab_dynamic = {"cumulative_cases_averted": 3}

        table = build_distributional_validation_table(
            python_summary,
            matlab_summary,
            python_dynamic,
            matlab_dynamic,
            summary_metrics=["nScreened"],
        )

        self.assertEqual(
            list(table.columns),
            [
                "Metric",
                "PythonMedian",
                "MatlabMedian",
                "AbsoluteDifference",
                "RelativeDifference",
                "Tolerance",
                "Pass",
                "Notes",
            ],
        )
        self.assertEqual(len(table), 5)

    def test_reference_validation_runs_with_small_override(self) -> None:
        out = run_reference_distributional_validation(
            FIXTURE_DIR,
            config_overrides={"N": 100, "nReps": 2},
            summary_metrics=["nScreened", "nPreventedActiveTB"],
        )

        self.assertIn("validation", out)
        self.assertEqual(out["config"]["seed"], 1)
        self.assertEqual(out["config"]["screeningStrategy"], "random")
        self.assertEqual(len(out["validation"]), 6)


if __name__ == "__main__":
    unittest.main()
