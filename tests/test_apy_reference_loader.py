from __future__ import annotations

import json
import shutil
import unittest
import uuid
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import pandas as pd

from engine.apy.reference_loader import (
    load_reference_dir,
    load_reference_bundle,
    load_reference_dynamic_comparison,
    load_reference_manifest,
    load_reference_scenario_config,
    load_reference_summary,
)


class ApyReferenceLoaderTests(unittest.TestCase):
    @contextmanager
    def temporary_workspace(self) -> Iterator[Path]:
        path = Path.cwd() / f".tmp_apy_reference_loader_{uuid.uuid4().hex}"
        path.mkdir()
        try:
            yield path
        finally:
            shutil.rmtree(path, ignore_errors=True)

    def test_reference_bundle_json_loads(self) -> None:
        with self.temporary_workspace() as tmp:
            path = tmp / "bundle.json"
            path.write_text(json.dumps({"metadata": {"modelVersion": "v9"}}), encoding="utf-8")

            bundle = load_reference_bundle(path)

        self.assertEqual(bundle["metadata"]["modelVersion"], "v9")

    def test_reference_manifest_json_loads(self) -> None:
        with self.temporary_workspace() as tmp:
            path = tmp / "manifest.json"
            path.write_text(json.dumps({"scenario_id": "fixture"}), encoding="utf-8")

            manifest = load_reference_manifest(path)

        self.assertEqual(manifest["scenario_id"], "fixture")

    def test_reference_scenario_config_json_loads(self) -> None:
        with self.temporary_workspace() as tmp:
            path = tmp / "scenario_config.json"
            path.write_text(json.dumps({"N": 1500}), encoding="utf-8")

            config = load_reference_scenario_config(path)

        self.assertEqual(config["N"], 1500)

    def test_reference_summary_csv_loads(self) -> None:
        with self.temporary_workspace() as tmp:
            path = tmp / "summary.csv"
            pd.DataFrame(
                [{"Metric": "nPreventedActiveTB", "Median": 12.0}]
            ).to_csv(path, index=False)

            summary = load_reference_summary(path)

        self.assertEqual(summary.loc[0, "Metric"], "nPreventedActiveTB")
        self.assertEqual(float(summary.loc[0, "Median"]), 12.0)

    def test_reference_dynamic_comparison_json_loads(self) -> None:
        with self.temporary_workspace() as tmp:
            path = tmp / "dynamic_comparison.json"
            path.write_text(
                json.dumps({"available": True, "source": "doNothing.derived"}),
                encoding="utf-8",
            )

            comparison = load_reference_dynamic_comparison(path)

        self.assertTrue(comparison["available"])
        self.assertEqual(comparison["source"], "doNothing.derived")

    def test_missing_bundle_raises_clear_error(self) -> None:
        with self.assertRaisesRegex(FileNotFoundError, "MATLAB reference bundle not found"):
            load_reference_bundle("missing_bundle.json")

    def test_missing_summary_raises_clear_error(self) -> None:
        with self.assertRaisesRegex(FileNotFoundError, "MATLAB reference summary CSV not found"):
            load_reference_summary("missing_summary.csv")

    def test_missing_dynamic_comparison_raises_clear_error(self) -> None:
        with self.assertRaisesRegex(
            FileNotFoundError, "MATLAB dynamic-comparison reference not found"
        ):
            load_reference_dynamic_comparison("missing_dynamic_comparison.json")

    def test_committed_compact_reference_dir_loads(self) -> None:
        reference_dir = (
            Path("validation")
            / "matlab_reference"
            / "default_random_igra_3hp_N1500_seed1"
        )
        if not reference_dir.exists():
            self.skipTest("Compact MATLAB reference fixture has not been generated.")

        loaded = load_reference_dir(reference_dir)

        expected_files = set(loaded["manifest"]["expected_files"])
        excluded_files = set(loaded["manifest"]["excluded_large_files"])
        self.assertIn("scenario_config.json", expected_files)
        self.assertIn("matlab_dynamic_comparison.json", expected_files)
        self.assertIn("matlab_summary.csv", expected_files)
        self.assertIn("matlab_results_bundle.json", excluded_files)
        self.assertIn("matlab_raw_replicates.csv", excluded_files)
        self.assertEqual(loaded["scenario_config"]["N"], 1500)
        self.assertTrue(loaded["dynamic_comparison"]["available"])
        self.assertFalse(loaded["summary"].empty)
        self.assertFalse((reference_dir / "matlab_results_bundle.json").is_file())
        self.assertFalse((reference_dir / "matlab_raw_replicates.csv").is_file())


if __name__ == "__main__":
    unittest.main()
