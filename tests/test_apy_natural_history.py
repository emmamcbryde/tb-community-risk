from __future__ import annotations

import math
import json
import unittest

import numpy as np
import pandas as pd

from engine.apy.natural_history import (
    build_natural_history_addon_report,
    run_do_nothing_summary,
    safe_fraction_vector,
)
from engine.apy.results_bundle import build_results_bundle
from engine.apy.runner import run_replicates


class ApyNaturalHistoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = run_replicates(n=100, n_reps=3, seed=30)
        cls.do_nothing = run_do_nothing_summary(cls.results)

    def test_run_do_nothing_summary_returns_expected_sections(self) -> None:
        self.assertIn("summary", self.do_nothing)
        self.assertIn("derived", self.do_nothing)
        self.assertIn("resultsUsed", self.do_nothing)
        self.assertIs(self.do_nothing["resultsUsed"], self.results)

    def test_summary_dataframe_has_expected_columns(self) -> None:
        summary = self.do_nothing["summary"]

        self.assertEqual(
            list(summary.columns),
            ["Metric", "Median", "Low95", "High95", "DisplayScale"],
        )

    def test_derived_dataframe_has_matlab_compatible_fields(self) -> None:
        derived = self.do_nothing["derived"]
        required = {
            "nInfected",
            "nActiveBy2y_DoNothing",
            "nActiveBy20y_DoNothing",
            "nActiveBy20y_AfterStrategy",
            "nActiveBy20y_Prevented",
            "relReduction20y",
            "ltbiPrev_DoNothing",
            "activePrev2y_DoNothing",
            "activePrev20y_DoNothing",
            "activePrev20y_AfterStrategy",
        }

        self.assertTrue(required.issubset(set(derived.columns)))

    def test_after_strategy_cases_are_baseline_minus_prevented(self) -> None:
        derived = self.do_nothing["derived"]
        expected = (
            derived["nActiveBy20y_DoNothing"]
            - derived["nActiveBy20y_Prevented"]
        )

        pd.testing.assert_series_equal(
            derived["nActiveBy20y_AfterStrategy"],
            expected,
            check_names=False,
        )

    def test_relative_reduction_uses_safe_fraction(self) -> None:
        derived = self.do_nothing["derived"]
        expected = safe_fraction_vector(
            derived["nActiveBy20y_Prevented"],
            derived["nActiveBy20y_DoNothing"],
        )

        np.testing.assert_allclose(
            derived["relReduction20y"].to_numpy(),
            expected,
            equal_nan=True,
        )
        self.assertTrue(math.isnan(safe_fraction_vector([1], [0])[0]))

    def test_prevalence_fields_use_config_population(self) -> None:
        derived = self.do_nothing["derived"]
        n = self.results["interfaceConfig"]["N"]

        np.testing.assert_allclose(
            derived["ltbiPrev_DoNothing"],
            derived["nInfected"] / n,
        )
        np.testing.assert_allclose(
            derived["activePrev20y_AfterStrategy"],
            derived["nActiveBy20y_AfterStrategy"] / n,
        )

    def test_nns_rows_are_present_when_source_fields_exist(self) -> None:
        metrics = set(self.do_nothing["summary"]["Metric"])

        self.assertIn(
            "NNS to prevent one active TB case (current strategy)",
            metrics,
        )
        self.assertIn(
            "NNS to cure one infection (current strategy)",
            metrics,
        )

    def test_addon_report_contract_is_json_serialisable(self) -> None:
        bundle = build_results_bundle(self.results, do_nothing=self.do_nothing)
        payload = build_natural_history_addon_report(
            bundle,
            do_nothing=self.do_nothing,
            requested_attributable_metrics=[],
        )

        self.assertEqual(
            list(payload),
            [
                "status",
                "source",
                "summaryRows",
                "doNothing",
                "attributableRisk",
                "missingInputs",
                "unsupportedMetrics",
                "messages",
            ],
        )
        self.assertEqual(payload["status"], "available")
        json.dumps(payload, allow_nan=False, sort_keys=True)

    def test_addon_report_converts_do_nothing_dataframes_to_rows(self) -> None:
        bundle = build_results_bundle(self.results, do_nothing=self.do_nothing)
        payload = build_natural_history_addon_report(
            bundle,
            do_nothing=self.do_nothing,
            requested_attributable_metrics=[],
        )

        expected_summary = self.do_nothing["summary"].to_dict(orient="records")
        expected_derived_columns = list(self.do_nothing["derived"].columns)

        self.assertEqual(payload["summaryRows"], expected_summary)
        self.assertEqual(payload["doNothing"]["summaryRows"], expected_summary)
        self.assertEqual(
            len(payload["doNothing"]["derivedRows"]),
            len(self.do_nothing["derived"]),
        )
        self.assertEqual(
            list(payload["doNothing"]["derivedRows"][0]),
            expected_derived_columns,
        )
        self.assertTrue(
            any(
                row["NNS_preventActiveTB"] is None
                for row in payload["doNothing"]["derivedRows"]
            )
        )

    def test_addon_report_reports_missing_do_nothing_input(self) -> None:
        bundle = build_results_bundle(self.results, do_nothing=self.do_nothing)
        payload = build_natural_history_addon_report(
            bundle,
            requested_attributable_metrics=[],
        )

        self.assertEqual(payload["status"], "missing-input")
        self.assertEqual(payload["summaryRows"], [])
        self.assertEqual(payload["doNothing"]["status"], "missing-input")
        self.assertIn(
            {
                "field": "doNothing",
                "message": (
                    "Do-nothing natural-history summary was not provided; "
                    "run_do_nothing_summary output is required for these rows."
                ),
            },
            payload["missingInputs"],
        )
        json.dumps(payload, allow_nan=False, sort_keys=True)

    def test_addon_report_preserves_attributable_risk_missing_and_unsupported_status(
        self,
    ) -> None:
        payload = build_natural_history_addon_report(
            {"technical": {}},
            do_nothing=self.do_nothing,
            requested_attributable_metrics=["ExpectedAttributableCases20y_Per1500"],
        )

        self.assertEqual(payload["status"], "missing-input")
        self.assertEqual(
            payload["attributableRisk"]["status"],
            "missing-input",
        )
        self.assertEqual(
            [row["field"] for row in payload["missingInputs"]],
            ["technical.dynamicComparison"],
        )
        self.assertEqual(
            [row["metric"] for row in payload["unsupportedMetrics"]],
            ["ExpectedAttributableCases20y_Per1500"],
        )


if __name__ == "__main__":
    unittest.main()
