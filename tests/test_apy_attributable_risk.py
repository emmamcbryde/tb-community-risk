from __future__ import annotations

import builtins
import json
import sys
import unittest

from engine.apy.attributable_risk import (
    KNOWN_MATLAB_ATTRIBUTABLE_METRICS,
    run_attributable_risk,
)


class ApyAttributableRiskTests(unittest.TestCase):
    def test_contract_is_json_serialisable_with_exact_top_level_keys(self) -> None:
        payload = run_attributable_risk({})

        self.assertEqual(
            set(payload),
            {
                "status",
                "source",
                "calculatedRows",
                "missingInputs",
                "unsupportedMetrics",
                "messages",
            },
        )
        json.dumps(payload, allow_nan=False, sort_keys=True)

    def test_missing_dynamic_comparison_is_reported(self) -> None:
        payload = run_attributable_risk({"technical": {}})

        self.assertEqual(payload["status"], "missing-input")
        self.assertEqual(payload["missingInputs"][0]["field"], "technical.dynamicComparison")
        self.assertIn("technical.dynamicComparison", payload["missingInputs"][0]["message"])
        self.assertEqual(payload["calculatedRows"], [])

    def test_known_matlab_attributable_metrics_are_explicitly_unsupported(self) -> None:
        payload = run_attributable_risk(
            {"technical": {"dynamicComparison": {"available": True}}}
        )

        unsupported = {row["metric"] for row in payload["unsupportedMetrics"]}
        self.assertEqual(unsupported, set(KNOWN_MATLAB_ATTRIBUTABLE_METRICS))
        self.assertTrue(
            all("not implemented" in row["reason"] for row in payload["unsupportedMetrics"])
        )

    def test_requested_metrics_limit_the_unsupported_list(self) -> None:
        payload = run_attributable_risk(
            {"technical": {"dynamicComparison": {"available": True}}},
            requested_metrics=["ExpectedAttributableCases20y_Per1500"],
        )

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

    def test_does_not_import_matlab_engine(self) -> None:
        sys.modules.pop("matlab.engine", None)
        original_import = builtins.__import__

        def fail_on_matlab_engine_import(name, *args, **kwargs):
            if name == "matlab.engine":
                raise AssertionError("attributable-risk scaffold imported matlab.engine")
            return original_import(name, *args, **kwargs)

        try:
            builtins.__import__ = fail_on_matlab_engine_import
            payload = run_attributable_risk({})
        finally:
            builtins.__import__ = original_import

        self.assertNotIn("matlab.engine", sys.modules)
        self.assertEqual(payload["status"], "missing-input")

    def test_dynamic_comparison_active_tb_fields_do_not_create_rows(self) -> None:
        payload = run_attributable_risk(
            {
                "technical": {
                    "dynamicComparison": {
                        "available": True,
                        "cumulative_baseline_active_tb_cases": 12.0,
                        "cumulative_intervention_active_tb_cases": 9.0,
                        "cumulative_cases_averted": 3.0,
                        "relative_reduction_cumulative_active_tb_cases": 0.25,
                    }
                }
            }
        )

        unsupported = {row["metric"] for row in payload["unsupportedMetrics"]}
        self.assertEqual(payload["status"], "unsupported")
        self.assertEqual(payload["missingInputs"], [])
        self.assertEqual(payload["calculatedRows"], [])
        self.assertNotIn("cumulative_baseline_active_tb_cases", unsupported)
        self.assertNotIn("cumulative_cases_averted", unsupported)

    def test_precomputed_calculated_rows_are_passed_through_without_formulas(self) -> None:
        rows = [
            {
                "metric": "simpleArtificialMetric",
                "value": 7.5,
                "status": "precomputed",
                "notes": ("from", "bundle"),
            }
        ]

        payload = run_attributable_risk(
            {
                "technical": {"dynamicComparison": {"available": True}},
                "attributableRisk": {"calculatedRows": rows},
            },
            requested_metrics=["ExpectedAttributableCases20y_Per1500"],
        )

        self.assertEqual(payload["status"], "partial")
        self.assertEqual(
            payload["calculatedRows"],
            [
                {
                    "metric": "simpleArtificialMetric",
                    "value": 7.5,
                    "status": "precomputed",
                    "notes": ["from", "bundle"],
                }
            ],
        )
        self.assertEqual(payload["missingInputs"], [])
        self.assertEqual(
            [row["metric"] for row in payload["unsupportedMetrics"]],
            ["ExpectedAttributableCases20y_Per1500"],
        )
        self.assertTrue(any("passed through" in message for message in payload["messages"]))
        self.assertTrue(any("unsupported" in message for message in payload["messages"]))
        json.dumps(payload, allow_nan=False, sort_keys=True)

    def test_precomputed_attributable_rows_are_fallback_rows(self) -> None:
        payload = run_attributable_risk(
            {
                "technical": {"dynamicComparison": {"available": True}},
                "attributableRisk": {
                    "attributableRows": [
                        {
                            "Metric": "Artificial attributable row",
                            "Value": 3,
                        }
                    ]
                },
            },
            requested_metrics=[],
        )

        self.assertEqual(payload["status"], "partial")
        self.assertEqual(
            payload["calculatedRows"],
            [{"Metric": "Artificial attributable row", "Value": 3}],
        )
        self.assertEqual(payload["unsupportedMetrics"], [])
        self.assertTrue(
            any("attributableRisk.attributableRows" in message for message in payload["messages"])
        )
        json.dumps(payload, allow_nan=False, sort_keys=True)

    def test_non_empty_attributable_rows_win_over_empty_calculated_rows(self) -> None:
        payload = run_attributable_risk(
            {
                "technical": {"dynamicComparison": {"available": True}},
                "attributableRisk": {
                    "calculatedRows": [],
                    "attributableRows": [{"metric": "fallbackRow", "value": 4}],
                },
            },
            requested_metrics=[],
        )

        self.assertEqual(payload["status"], "partial")
        self.assertEqual(payload["calculatedRows"], [{"metric": "fallbackRow", "value": 4}])


if __name__ == "__main__":
    unittest.main()
