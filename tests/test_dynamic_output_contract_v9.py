from __future__ import annotations

import json
import unittest

import numpy as np
import pandas as pd

from engine.integration.dynamic_output_contract_v9 import build_dynamic_results_bundle_v9


class DynamicOutputContractV9Tests(unittest.TestCase):
    def test_fake_projection_builds_bundle(self) -> None:
        df_future = pd.DataFrame(
            {
                "Year": [1, 2],
                "Baseline_inc_per100k": [20.0, 18.0],
                "Intervention_inc_per100k": [19.0, 15.0],
                "Baseline_annual_count": [10.0, 9.0],
                "Intervention_annual_count": [9.5, 7.5],
                "Baseline_cum_count": [10.0, 19.0],
                "Intervention_cum_count": [9.5, 17.0],
                "Cases_averted_cum_count": [0.5, 2.0],
            }
        )
        bundle = build_dynamic_results_bundle_v9(
            df_future,
            params_base={"beta": np.float64(1.2)},
            metadata={"population": 50000},
        )

        self.assertEqual(bundle["contractVersion"], "dynamic_output_contract_v9")
        self.assertTrue(bundle["headline"]["summaryRows"])
        self.assertEqual(len(bundle["projection"]["annualRows"]), 2)

    def test_cumulative_cases_averted_is_calculated_when_missing(self) -> None:
        df_future = pd.DataFrame(
            {
                "Year": [1, 2],
                "Baseline_cum_count": [10.0, 25.0],
                "Intervention_cum_count": [8.0, 19.0],
            }
        )
        bundle = build_dynamic_results_bundle_v9(df_future)
        metrics = {row["Metric"]: row["Value"] for row in bundle["headline"]["summaryRows"]}

        self.assertEqual(metrics["cumulative_cases_averted"], 6)
        self.assertAlmostEqual(metrics["relative_reduction_cumulative_active_tb_cases"], 0.24)

    def test_empty_dataframe_is_handled(self) -> None:
        bundle = build_dynamic_results_bundle_v9(pd.DataFrame())

        self.assertEqual(bundle["projection"]["annualRows"], [])
        self.assertEqual(bundle["headline"]["keyMetricsRows"], [])

    def test_two_epoch_calibration_diagnostics_are_preserved_json_safe(self) -> None:
        df_future = pd.DataFrame(
            {
                "Year": [1],
                "Baseline_cum_count": [10.0],
                "Intervention_cum_count": [8.0],
            }
        )
        calibration = {
            "calibration_mode": "Two-epoch beta: historical + recent 10 years",
            "beta_historical": np.float64(2.0),
            "beta_recent": np.float64(1.0),
            "beta_forward": np.float64(1.0),
            "beta_series": [2.0, 2.0, 1.0, 1.0],
            "calibration_years": [2021, 2022, 2023, 2024],
            "beta_historical_years": [2021, 2022],
            "beta_recent_years": [2023, 2024],
            "fitted_incidence": [30.0, 28.0, 20.0, 18.0],
            "target_incidence": [31.0, 27.0, 21.0, 17.0],
            "residuals": [-1.0, 1.0, -1.0, 1.0],
        }

        bundle = build_dynamic_results_bundle_v9(df_future, calibration=calibration)
        technical_calibration = bundle["technical"]["calibration"]

        for key in [
            "beta_series",
            "calibration_years",
            "beta_historical_years",
            "beta_recent_years",
            "fitted_incidence",
            "target_incidence",
            "residuals",
        ]:
            self.assertIn(key, technical_calibration)

        self.assertEqual(technical_calibration["beta_series"], [2.0, 2.0, 1.0, 1.0])
        self.assertEqual(
            technical_calibration["calibration_years"], [2021, 2022, 2023, 2024]
        )
        json.dumps(bundle, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
