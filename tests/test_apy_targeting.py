from __future__ import annotations

import json
import subprocess
import sys
import unittest
from types import SimpleNamespace

from engine.apy.targeting import (
    build_targeting_comparison,
    compare_summary_metrics,
    extract_targeting_metadata,
)


class ApyTargetingTests(unittest.TestCase):
    def test_metadata_extracts_common_fields_from_camel_and_snake_aliases(self) -> None:
        strategy = {
            "strategyLabel": "Target high risk",
            "screening_test": "IGRA",
            "treatmentRegimen": "3HP",
        }
        config = SimpleNamespace(
            screeningCoverage=0.75,
            sample_size=5000,
            follow_horizon=20,
            randomSeed=7,
        )

        metadata = extract_targeting_metadata(strategy=strategy, config=config)

        self.assertEqual(metadata["strategy"], "Target high risk")
        self.assertEqual(metadata["screeningTest"], "IGRA")
        self.assertEqual(metadata["regimen"], "3HP")
        self.assertEqual(metadata["coverage"], 0.75)
        self.assertEqual(metadata["population"], 5000)
        self.assertEqual(metadata["followHorizon"], 20)
        self.assertEqual(metadata["seed"], 7)

    def test_comparison_rows_are_strict_json_serialisable(self) -> None:
        comparison = build_targeting_comparison(
            {"cases": {"Median": 10}, "cost": 100.5, "bad": float("nan")},
            ["cases", "cost", "bad", "overflow"],
            baseline={
                "cases": {"Median": 8},
                "cost": 90,
                "bad": float("inf"),
                "overflow": float("-inf"),
            },
            strategy={"strategy": "random", "coverage": float("nan")},
            config={"N": 100, "followHorizon": float("inf")},
            include_relative=True,
        )

        encoded = json.dumps(comparison, allow_nan=False, sort_keys=True)

        self.assertIn('"metricRows"', encoded)

    def test_metrics_rejects_plain_string(self) -> None:
        with self.assertRaisesRegex(
            TypeError,
            "metrics must be a sequence of metric names, not a string",
        ):
            compare_summary_metrics({"cases": 10}, "cases")

    def test_numeric_metrics_get_absolute_differences(self) -> None:
        rows = compare_summary_metrics(
            {"nScreened": 125, "nPreventedActiveTB": {"Median": 4.5}},
            ["nScreened", "nPreventedActiveTB"],
            baseline={"nScreened": 100, "nPreventedActiveTB": {"Median": 1.5}},
        )

        rows_by_metric = {row["metric"]: row for row in rows}
        self.assertEqual(rows_by_metric["nScreened"]["status"], "ok")
        self.assertEqual(rows_by_metric["nScreened"]["absoluteDifference"], 25.0)
        self.assertEqual(
            rows_by_metric["nPreventedActiveTB"]["absoluteDifference"],
            3.0,
        )
        self.assertIsNone(rows_by_metric["nScreened"]["relativeDifference"])

    def test_relative_difference_is_explicit_and_requires_nonzero_baseline(self) -> None:
        without_relative = compare_summary_metrics(
            {"metric": 12},
            ["metric"],
            baseline={"metric": 8},
        )[0]
        with_relative = compare_summary_metrics(
            {"metric": 12},
            ["metric"],
            baseline={"metric": 8},
            include_relative=True,
        )[0]
        zero_baseline = compare_summary_metrics(
            {"metric": 12},
            ["metric"],
            baseline={"metric": 0},
            include_relative=True,
        )[0]

        self.assertIsNone(without_relative["relativeDifference"])
        self.assertEqual(with_relative["relativeDifference"], 0.5)
        self.assertIsNone(zero_baseline["relativeDifference"])

    def test_missing_and_non_numeric_metrics_are_reported_clearly(self) -> None:
        rows = compare_summary_metrics(
            {"unsupported": "not available", "bad": float("nan"), "ok": 3},
            ["missing", "unsupported", "bad", "ok"],
            baseline={"ok": "also unsupported"},
        )

        rows_by_metric = {row["metric"]: row for row in rows}
        self.assertEqual(rows_by_metric["missing"]["status"], "missing")
        self.assertIn("not present", rows_by_metric["missing"]["notes"])
        self.assertEqual(rows_by_metric["unsupported"]["status"], "non_numeric")
        self.assertIn("not a finite numeric", rows_by_metric["unsupported"]["notes"])
        self.assertEqual(rows_by_metric["bad"]["status"], "non_numeric")
        self.assertEqual(rows_by_metric["ok"]["baselineStatus"], "non_numeric")
        self.assertIn("Baseline metric", rows_by_metric["ok"]["notes"])

    def test_import_has_no_matlab_engine_dependency(self) -> None:
        command = (
            "import engine.apy.targeting, sys; "
            "print([name for name in sys.modules "
            "if name == 'matlab' or name.startswith('matlab.')])"
        )
        completed = subprocess.run(
            [sys.executable, "-c", command],
            check=True,
            capture_output=True,
            text=True,
        )

        self.assertEqual(completed.stdout.strip(), "[]")


if __name__ == "__main__":
    unittest.main()
