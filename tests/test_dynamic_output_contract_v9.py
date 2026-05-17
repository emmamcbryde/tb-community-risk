from __future__ import annotations

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


if __name__ == "__main__":
    unittest.main()
