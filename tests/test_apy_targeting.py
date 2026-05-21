from __future__ import annotations

import json
import subprocess
import sys
import unittest
from types import SimpleNamespace

from engine.apy.targeting import (
    apply_strategy_config_overrides,
    build_targeting_comparison,
    compare_summary_metrics,
    compare_targeting_result_bundles,
    extract_targeting_metadata,
    make_strategy_spec,
)


class ApyTargetingTests(unittest.TestCase):
    def test_strategy_spec_is_explicit_and_json_serialisable(self) -> None:
        spec = make_strategy_spec(
            "IGRA high-risk",
            {
                "screeningStrategy": "ltbi",
                "testType": "IGRA",
                "screenCoverage": 0.8,
            },
        )

        encoded = json.dumps(spec, allow_nan=False, sort_keys=True)

        self.assertEqual(spec["name"], "IGRA high-risk")
        self.assertEqual(
            spec["config_overrides"],
            {
                "screeningStrategy": "ltbi",
                "testType": "IGRA",
                "screenCoverage": 0.8,
            },
        )
        self.assertIn('"config_overrides"', encoded)

    def test_strategy_spec_rejects_implicit_or_non_json_values(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-empty string"):
            make_strategy_spec("", {"screeningStrategy": "random"})

        with self.assertRaisesRegex(ValueError, "finite numeric"):
            make_strategy_spec("bad", {"screenCoverage": float("nan")})

        with self.assertRaisesRegex(TypeError, "JSON-compatible"):
            make_strategy_spec("bad", {"unsupported": object()})

    def test_apply_strategy_config_overrides_does_not_mutate_base_config(self) -> None:
        base_config = {
            "screeningStrategy": "random",
            "testType": "TST",
            "nested": {"kept": True},
        }
        spec = make_strategy_spec(
            "IGRA targeted",
            {
                "screeningStrategy": "ltbi",
                "testType": "IGRA",
                "nested": {"changed": True},
            },
        )

        config = apply_strategy_config_overrides(base_config, spec)
        config["nested"]["changed"] = False

        self.assertEqual(base_config["screeningStrategy"], "random")
        self.assertEqual(base_config["testType"], "TST")
        self.assertEqual(base_config["nested"], {"kept": True})
        self.assertEqual(config["screeningStrategy"], "ltbi")
        self.assertEqual(config["testType"], "IGRA")

    def test_applied_strategy_config_preserves_explicit_strategy_spec(self) -> None:
        base_config = {
            "screeningStrategy": "random",
            "strategySpec": {
                "name": "Base strategy",
                "config_overrides": {"nested": {"base": True}},
            },
            "nested": {"kept": True},
        }
        spec = make_strategy_spec(
            "Explicit high-risk",
            {
                "screeningStrategy": "ltbi",
                "testType": "IGRA",
                "nested": {"changed": True},
            },
        )

        config = apply_strategy_config_overrides(base_config, spec)
        metadata = extract_targeting_metadata(config=config)
        config["strategySpec"]["config_overrides"]["nested"]["changed"] = False
        config["nested"]["changed"] = False

        self.assertEqual(metadata["strategy"], "Explicit high-risk")
        self.assertEqual(config["strategySpec"]["name"], "Explicit high-risk")
        self.assertIsNot(config["strategySpec"], spec)
        self.assertEqual(
            base_config,
            {
                "screeningStrategy": "random",
                "strategySpec": {
                    "name": "Base strategy",
                    "config_overrides": {"nested": {"base": True}},
                },
                "nested": {"kept": True},
            },
        )
        self.assertEqual(
            spec,
            {
                "name": "Explicit high-risk",
                "config_overrides": {
                    "screeningStrategy": "ltbi",
                    "testType": "IGRA",
                    "nested": {"changed": True},
                },
            },
        )

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

    def test_metadata_extracts_strategy_spec_name_from_config_dict(self) -> None:
        config = {
            "strategySpec": make_strategy_spec(
                "Explicit strategy",
                {"screeningStrategy": "prevent"},
            ),
            "testType": "IGRA",
            "screenCoverage": 0.9,
            "N": 1200,
            "followHorizon": 15,
            "seed": 11,
        }

        metadata = extract_targeting_metadata(config=config)

        self.assertEqual(metadata["strategy"], "Explicit strategy")
        self.assertEqual(metadata["screeningTest"], "IGRA")
        self.assertEqual(metadata["coverage"], 0.9)
        self.assertEqual(metadata["population"], 1200)
        self.assertEqual(metadata["followHorizon"], 15)
        self.assertEqual(metadata["seed"], 11)

    def test_metadata_extracts_from_runner_result_bundle_shape(self) -> None:
        result_bundle = {
            "results": {
                "strategy": {
                    "screeningStrategy": "cure",
                    "testType": "TST",
                    "regimen": "6H",
                    "screenCoverage": 0.65,
                },
                "interfaceConfig": {"N": 900, "followHorizon": 20, "seed": 5},
            }
        }

        metadata = extract_targeting_metadata(result_bundle=result_bundle)

        self.assertEqual(metadata["strategy"], "cure")
        self.assertEqual(metadata["screeningTest"], "TST")
        self.assertEqual(metadata["regimen"], "6H")
        self.assertEqual(metadata["coverage"], 0.65)
        self.assertEqual(metadata["population"], 900)
        self.assertEqual(metadata["followHorizon"], 20)
        self.assertEqual(metadata["seed"], 5)

    def test_metadata_extracts_from_results_bundle_output_shape(self) -> None:
        result_bundle = {
            "headline": {
                "strategy": {
                    "screeningStrategy": "ltbi",
                    "testType": "IGRA",
                    "regimen": "3HP",
                    "screenCoverage": 0.7,
                },
            },
            "technical": {
                "interfaceConfig": {"N": 1000, "followHorizon": 10, "seed": 3},
            },
        }

        metadata = extract_targeting_metadata(result_bundle=result_bundle)

        self.assertEqual(metadata["strategy"], "ltbi")
        self.assertEqual(metadata["screeningTest"], "IGRA")
        self.assertEqual(metadata["regimen"], "3HP")
        self.assertEqual(metadata["coverage"], 0.7)
        self.assertEqual(metadata["population"], 1000)
        self.assertEqual(metadata["followHorizon"], 10)
        self.assertEqual(metadata["seed"], 3)

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

    def test_bundle_comparison_uses_existing_summary_rows_and_metadata(self) -> None:
        baseline = {
            "headline": {
                "strategy": {"screeningStrategy": "baseline"},
                "summaryRows": [
                    {"Metric": "nScreened", "Median": 100},
                    {"Metric": "nPreventedActiveTB", "Median": 4},
                ],
            },
        }
        comparator = {
            "headline": {
                "strategy": {"strategyLabel": "targeted"},
                "summaryRows": [
                    {"Metric": "nScreened", "Median": 125},
                    {"Metric": "nPreventedActiveTB", "Median": 6},
                ],
            },
        }

        rows = compare_targeting_result_bundles(
            baseline,
            [comparator],
            ["nScreened", "nPreventedActiveTB"],
        )

        rows_by_metric = {row["metric"]: row for row in rows}
        self.assertEqual(rows_by_metric["nScreened"]["baselineStrategy"], "baseline")
        self.assertEqual(rows_by_metric["nScreened"]["comparatorStrategy"], "targeted")
        self.assertEqual(rows_by_metric["nScreened"]["baselineValue"], 100.0)
        self.assertEqual(rows_by_metric["nScreened"]["comparatorValue"], 125.0)
        self.assertEqual(rows_by_metric["nScreened"]["absoluteDifference"], 25.0)
        self.assertEqual(rows_by_metric["nScreened"]["relativeDifference"], 0.25)
        self.assertEqual(
            rows_by_metric["nPreventedActiveTB"]["absoluteDifference"],
            2.0,
        )
        json.dumps(rows, allow_nan=False, sort_keys=True)

    def test_bundle_comparison_reads_results_summary_without_running_models(self) -> None:
        baseline = {
            "results": {
                "strategy": {"screeningStrategy": "baseline"},
                "summary": [{"Metric": "nScreened", "Median": 0}],
            },
        }
        comparator = {
            "results": {
                "strategy": {"screeningStrategy": "comparator"},
                "summary": [{"Metric": "nScreened", "Median": 12}],
            },
        }

        row = compare_targeting_result_bundles(
            baseline,
            [comparator],
            ["nScreened"],
        )[0]

        self.assertEqual(row["baselineValue"], 0.0)
        self.assertEqual(row["comparatorValue"], 12.0)
        self.assertEqual(row["absoluteDifference"], 12.0)
        self.assertIsNone(row["relativeDifference"])

    def test_bundle_comparison_reports_missing_and_unsupported_metrics(self) -> None:
        rows = compare_targeting_result_bundles(
            {"headline": {"strategy": {"screeningStrategy": "baseline"}}},
            [
                {
                    "headline": {
                        "strategy": {"screeningStrategy": "comparator"},
                        "summaryRows": [
                            {"Metric": "unsupported", "Value": "not available"},
                        ],
                    },
                }
            ],
            ["missing", "unsupported"],
        )

        rows_by_metric = {row["metric"]: row for row in rows}
        self.assertEqual(rows_by_metric["missing"]["status"], "missing")
        self.assertEqual(rows_by_metric["missing"]["baselineStatus"], "missing")
        self.assertEqual(rows_by_metric["unsupported"]["status"], "non_numeric")
        self.assertEqual(rows_by_metric["unsupported"]["baselineStatus"], "missing")
        self.assertIsNone(rows_by_metric["unsupported"]["comparatorValue"])

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
