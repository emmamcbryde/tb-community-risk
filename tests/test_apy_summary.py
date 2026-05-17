from __future__ import annotations

import unittest

import pandas as pd

from engine.apy.summary import empirical_quantile, summarise_numeric_rows


class ApySummaryTests(unittest.TestCase):
    def test_empirical_quantile_matches_matlab_style_interpolation(self) -> None:
        values = [1, 2, 3, 4]

        self.assertAlmostEqual(empirical_quantile(values, 0.25), 1.75)
        self.assertAlmostEqual(empirical_quantile(values, 0.5), 2.5)
        self.assertAlmostEqual(empirical_quantile(values, 0.975), 3.925)

    def test_numeric_summary_rows_include_expected_columns(self) -> None:
        raw = pd.DataFrame(
            {
                "numeric_metric": [1, 2, 3, 4],
                "logical_metric": [True, False, True, True],
                "label": ["a", "b", "c", "d"],
            }
        )

        summary = summarise_numeric_rows(raw)

        self.assertEqual(list(summary.columns), ["Metric", "Median", "Low95", "High95"])
        self.assertEqual(summary["Metric"].tolist(), ["numeric_metric", "logical_metric"])
        numeric_row = summary.loc[summary["Metric"] == "numeric_metric"].iloc[0]
        self.assertEqual(float(numeric_row["Median"]), 2.5)
        self.assertAlmostEqual(float(numeric_row["Low95"]), 1.075)
        self.assertAlmostEqual(float(numeric_row["High95"]), 3.925)


if __name__ == "__main__":
    unittest.main()
