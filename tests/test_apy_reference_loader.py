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
    load_reference_bundle,
    load_reference_dynamic_comparison,
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


if __name__ == "__main__":
    unittest.main()
