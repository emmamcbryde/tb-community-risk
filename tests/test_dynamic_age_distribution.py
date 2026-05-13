from __future__ import annotations

import unittest

import pandas as pd

from ui.dynamic_ui import (
    default_five_year_age_distribution,
    expand_five_year_age_distribution,
    normalise_age_distribution,
)


class DynamicAgeDistributionTests(unittest.TestCase):
    def test_default_five_year_distribution_sums_to_one(self) -> None:
        df = default_five_year_age_distribution()

        self.assertEqual(len(df), 21)
        self.assertAlmostEqual(float(df["Proportion"].sum()), 1.0)
        self.assertEqual(df.iloc[-1]["AgeGroup"], "100+")

    def test_manual_five_year_distribution_expands_to_single_year_ages(self) -> None:
        df = pd.DataFrame(
            [
                {"AgeGroup": "0-4", "AgeStart": 0, "AgeEnd": 4, "Proportion": 2.0},
                {"AgeGroup": "5-9", "AgeStart": 5, "AgeEnd": 9, "Proportion": 1.0},
            ]
        )
        expanded = expand_five_year_age_distribution(df, population=300)

        self.assertEqual(expanded["age"].tolist(), list(range(10)))
        self.assertAlmostEqual(float(expanded["population"].sum()), 300.0)
        self.assertAlmostEqual(float(expanded.loc[expanded["age"] == 0, "population"].iloc[0]), 40.0)
        self.assertAlmostEqual(float(expanded.loc[expanded["age"] == 5, "population"].iloc[0]), 20.0)

    def test_normalise_clips_negative_values(self) -> None:
        df = pd.DataFrame(
            [
                {"AgeGroup": "0-4", "AgeStart": 0, "AgeEnd": 4, "Proportion": -1.0},
                {"AgeGroup": "5-9", "AgeStart": 5, "AgeEnd": 9, "Proportion": 2.0},
            ]
        )
        normalised = normalise_age_distribution(df)

        self.assertEqual(float(normalised.iloc[0]["Proportion"]), 0.0)
        self.assertEqual(float(normalised.iloc[1]["Proportion"]), 1.0)


if __name__ == "__main__":
    unittest.main()
