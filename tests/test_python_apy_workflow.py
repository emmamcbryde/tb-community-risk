from __future__ import annotations

import builtins
import copy
import importlib
import json
import sys
import unittest
from pathlib import Path

from engine.apy.economics_compare import compare_economics_rows, economics_rows


class PythonApyWorkflowTests(unittest.TestCase):
    def test_small_python_only_economics_compare_workflow(self) -> None:
        module_sentinel = object()
        original_matlab_module = sys.modules.pop("matlab", module_sentinel)
        original_matlab_engine_module = sys.modules.pop("matlab.engine", module_sentinel)
        original_import = builtins.__import__

        def fail_on_matlab_import(name, *args, **kwargs):
            if name == "matlab" or name.startswith("matlab."):
                raise AssertionError(f"Python APY workflow imported {name}")
            return original_import(name, *args, **kwargs)

        try:
            builtins.__import__ = fail_on_matlab_import
            backend_module = importlib.import_module("adapters.python_apy_backend")
            backend = backend_module.PythonApyBackend(Path(__file__).resolve().parents[1])

            baseline_config = backend.default_config()
            baseline_config.update({"N": 50, "nReps": 1, "seed": 11})

            comparator_config = copy.deepcopy(baseline_config)
            comparator_config["scenarioLabel"] = "Comparator"
            comparator_config["screenCoverage"] = baseline_config["screenCoverage"] * 0.5

            economics_config = backend.economics_preset_kwab150()
            baseline = backend.run_economics_for_config(
                baseline_config,
                economics_config,
            )
            comparator = backend.run_economics_for_config(
                comparator_config,
                economics_config,
            )
        finally:
            builtins.__import__ = original_import
            if original_matlab_module is module_sentinel:
                sys.modules.pop("matlab", None)
            else:
                sys.modules["matlab"] = original_matlab_module
            if original_matlab_engine_module is module_sentinel:
                sys.modules.pop("matlab.engine", None)
            else:
                sys.modules["matlab.engine"] = original_matlab_engine_module

        json.dumps(baseline)
        json.dumps(comparator)

        self.assertEqual(baseline["status"], "partial")
        self.assertEqual(comparator["status"], "partial")

        baseline_rows = economics_rows(baseline)
        comparator_rows = economics_rows(comparator)
        self.assertGreater(len(baseline_rows), 0)
        self.assertGreater(len(comparator_rows), 0)
        expected_row_fields = {
            "metric",
            "label",
            "value",
            "category",
            "status",
            "includedInTotal",
        }
        for rows in (baseline_rows, comparator_rows):
            for row in rows:
                self.assertEqual(set(row), expected_row_fields)
            self.assertTrue(any(row["metric"] for row in rows))
            self.assertTrue(any(row["label"] for row in rows))

        compare_rows, warnings = compare_economics_rows(
            baseline_rows,
            comparator_rows,
        )

        self.assertEqual(warnings, [])
        self.assertGreater(len(compare_rows), 0)
        self.assertIn("testingCost", {row["metric"] for row in compare_rows})
        json.dumps(compare_rows)


if __name__ == "__main__":
    unittest.main()
