from __future__ import annotations

import builtins
import json
import sys
import unittest

from adapters.paths import repo_root
from adapters.python_apy_backend import PythonApyBackend


class PythonApyBackendAttributableRiskTests(unittest.TestCase):
    def setUp(self) -> None:
        self.backend = PythonApyBackend(repo_root())

    def test_accepts_minimal_app_style_result_bundle_and_returns_json(self) -> None:
        bundle = {
            "metadata": {"backend": "python"},
            "headline": {"summaryRows": []},
            "technical": {"dynamicComparison": {"available": True}},
        }

        payload = self.backend.run_attributable_risk(bundle)

        self.assertEqual(payload["status"], "unsupported")
        self.assertEqual(payload["missingInputs"], [])
        self.assertEqual(payload["calculatedRows"], [])
        self.assertIsInstance(payload["unsupportedMetrics"], list)
        json.dumps(payload, allow_nan=False, sort_keys=True)

    def test_reports_missing_technical_dynamic_comparison_clearly(self) -> None:
        payload = self.backend.run_attributable_risk(
            {
                "metadata": {"backend": "python"},
                "headline": {"summaryRows": []},
                "technical": {},
            }
        )

        self.assertEqual(payload["status"], "missing-input")
        self.assertEqual(payload["calculatedRows"], [])
        self.assertEqual(len(payload["missingInputs"]), 1)
        self.assertEqual(payload["missingInputs"][0]["field"], "technical.dynamicComparison")
        self.assertIn("technical.dynamicComparison", payload["missingInputs"][0]["message"])
        self.assertTrue(
            any("technical.dynamicComparison" in message for message in payload["messages"])
        )
        json.dumps(payload, allow_nan=False, sort_keys=True)

    def test_reports_current_unsupported_status_when_dynamic_comparison_is_present(
        self,
    ) -> None:
        payload = self.backend.run_attributable_risk(
            {"technical": {"dynamicComparison": {"available": True}}},
            requested_metrics=["ExpectedAttributableCases20y_Per1500"],
        )

        self.assertEqual(payload["status"], "unsupported")
        self.assertEqual(payload["missingInputs"], [])
        self.assertEqual(payload["calculatedRows"], [])
        self.assertEqual(
            payload["unsupportedMetrics"],
            [
                {
                    "metric": "ExpectedAttributableCases20y_Per1500",
                    "reason": (
                        "Reactivation attributable-risk calculation is not implemented "
                        "in the Python APY port yet."
                    ),
                    "matlabSource": "run_tb_screening_reactivation_attributable_v9",
                }
            ],
        )
        json.dumps(payload, allow_nan=False, sort_keys=True)

    def test_passes_through_precomputed_calculated_rows_from_artificial_bundle(self) -> None:
        payload = self.backend.run_attributable_risk(
            {
                "metadata": {"backend": "python"},
                "headline": {"summaryRows": []},
                "technical": {"dynamicComparison": {"available": True}},
                "attributableRisk": {
                    "calculatedRows": [
                        {
                            "metric": "artificialSimpleRow",
                            "value": 2.25,
                            "status": "precomputed",
                        }
                    ]
                },
            },
            requested_metrics=["PopulationAttributableFraction20y_ReactivationOnly"],
        )

        self.assertEqual(payload["status"], "partial")
        self.assertEqual(
            payload["calculatedRows"],
            [
                {
                    "metric": "artificialSimpleRow",
                    "value": 2.25,
                    "status": "precomputed",
                }
            ],
        )
        self.assertEqual(payload["missingInputs"], [])
        self.assertEqual(
            [row["metric"] for row in payload["unsupportedMetrics"]],
            ["PopulationAttributableFraction20y_ReactivationOnly"],
        )
        self.assertTrue(any("passed through" in message for message in payload["messages"]))
        self.assertTrue(any("unsupported" in message for message in payload["messages"]))
        json.dumps(payload, allow_nan=False, sort_keys=True)

    def test_passes_through_existing_attributable_rows_as_calculated_rows(self) -> None:
        payload = self.backend.run_attributable_risk(
            {
                "technical": {"dynamicComparison": {"available": True}},
                "attributableRisk": {
                    "attributableRows": [
                        {
                            "Metric": "Artificial existing row",
                            "Value": 9,
                        }
                    ]
                },
            },
            requested_metrics=[],
        )

        self.assertEqual(payload["status"], "partial")
        self.assertEqual(
            payload["calculatedRows"],
            [{"Metric": "Artificial existing row", "Value": 9}],
        )
        self.assertEqual(payload["unsupportedMetrics"], [])
        self.assertTrue(
            any("attributableRisk.attributableRows" in message for message in payload["messages"])
        )
        json.dumps(payload, allow_nan=False, sort_keys=True)

    def test_does_not_require_or_import_matlab_engine(self) -> None:
        sys.modules.pop("matlab.engine", None)
        original_import = builtins.__import__

        def fail_on_matlab_engine_import(name, *args, **kwargs):
            if name == "matlab.engine":
                raise AssertionError("PythonApyBackend.run_attributable_risk imported matlab.engine")
            return original_import(name, *args, **kwargs)

        try:
            builtins.__import__ = fail_on_matlab_engine_import
            payload = self.backend.run_attributable_risk(
                {"technical": {"dynamicComparison": {"available": True}}}
            )
        finally:
            builtins.__import__ = original_import

        self.assertNotIn("matlab.engine", sys.modules)
        self.assertEqual(payload["status"], "unsupported")


if __name__ == "__main__":
    unittest.main()
