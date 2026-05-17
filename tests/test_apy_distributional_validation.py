from __future__ import annotations

import json
from pathlib import Path
import shutil
import unittest

import pandas as pd

from engine.apy.distributional_validation import (
    build_distributional_validation_table,
    portable_config_from_reference,
    run_reference_distributional_validation,
    run_reference_suite_distributional_validation,
)
from engine.apy.reference_loader import load_reference_dir
from scripts.run_apy_distributional_validation import main as run_validation_cli


FIXTURE_DIR = (
    Path("validation") / "matlab_reference" / "default_random_igra_3hp_N1500_seed1"
)
REFERENCE_ROOT = Path("validation") / "matlab_reference"
SUITE_FILE = REFERENCE_ROOT / "scenario_suite_v1.json"


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

    def test_batch_validation_runs_on_committed_default_fixture(self) -> None:
        out = run_reference_suite_distributional_validation(
            REFERENCE_ROOT,
            suite_file=SUITE_FILE,
            scenario_ids=["default_random_igra_3hp_N1500_seed1"],
            config_overrides={"N": 100, "nReps": 2},
        )

        self.assertIn("scenarioRows", out)
        self.assertIn("metricRows", out)
        self.assertEqual(len(out["scenarioRows"]), 1)
        self.assertFalse(out["metricRows"].empty)

    def test_batch_validation_reports_missing_suite_fixtures(self) -> None:
        suite_path = Path("validation") / "output" / "test_missing_suite.json"
        suite_path.parent.mkdir(parents=True, exist_ok=True)
        suite_path.write_text(
            json.dumps(
                [
                    {
                        "scenario_id": "missing_fixture_for_test",
                        "description": "Missing fixture for unit test.",
                        "config_overrides": {},
                        "expected_focus": [],
                        "notes": "",
                    }
                ]
            ),
            encoding="utf-8",
        )
        try:
            out = run_reference_suite_distributional_validation(
                REFERENCE_ROOT,
                suite_file=suite_path,
                scenario_ids=["missing_fixture_for_test"],
                config_overrides={"N": 100, "nReps": 2},
            )

            row = out["scenarioRows"].iloc[0]
            self.assertEqual(row["n_metrics"], 0)
            self.assertIn("missing", row["notes"])
        finally:
            suite_path.unlink(missing_ok=True)

    def test_cli_runs_quick_default_fixture(self) -> None:
        out_dir = Path("validation") / "output" / "test_cli_distributional_validation"
        shutil.rmtree(out_dir, ignore_errors=True)
        try:
            exit_code = run_validation_cli(
                [
                    "--reference-root",
                    str(REFERENCE_ROOT),
                    "--suite-file",
                    str(SUITE_FILE),
                    "--scenario-id",
                    "default_random_igra_3hp_N1500_seed1",
                    "--quick",
                    "2",
                    "--output-dir",
                    str(out_dir),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue((out_dir / "scenario_summary.csv").is_file())
            self.assertTrue((out_dir / "metric_validation_rows.csv").is_file())
            self.assertTrue((out_dir / "validation_report.md").is_file())
        finally:
            shutil.rmtree(out_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
