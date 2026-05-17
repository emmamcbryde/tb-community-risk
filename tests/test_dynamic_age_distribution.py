from __future__ import annotations

import unittest

import pandas as pd

from ui.dynamic_age_distribution import (
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

    def test_non_numeric_and_missing_values_become_zero(self) -> None:
        df = pd.DataFrame(
            [
                {"AgeGroup": "0-4", "AgeStart": 0, "AgeEnd": 4, "Proportion": "bad"},
                {"AgeGroup": "5-9", "AgeStart": 5, "AgeEnd": 9, "Proportion": None},
                {"AgeGroup": "10-14", "AgeStart": 10, "AgeEnd": 14, "Proportion": 5},
            ]
        )
        normalised = normalise_age_distribution(df)

        self.assertEqual(float(normalised.iloc[0]["Proportion"]), 0.0)
        self.assertEqual(float(normalised.iloc[1]["Proportion"]), 0.0)
        self.assertEqual(float(normalised.iloc[2]["Proportion"]), 1.0)

    def test_values_that_sum_to_100_are_normalised(self) -> None:
        df = pd.DataFrame(
            [
                {"AgeGroup": "0-4", "AgeStart": 0, "AgeEnd": 4, "Proportion": 60},
                {"AgeGroup": "5-9", "AgeStart": 5, "AgeEnd": 9, "Proportion": 40},
            ]
        )
        normalised = normalise_age_distribution(df)

        self.assertAlmostEqual(float(normalised["Proportion"].sum()), 1.0)
        self.assertAlmostEqual(float(normalised.iloc[0]["Proportion"]), 0.6)
        self.assertAlmostEqual(float(normalised.iloc[1]["Proportion"]), 0.4)

    def test_all_zero_five_year_values_fall_back_safely(self) -> None:
        df = default_five_year_age_distribution()
        df["Proportion"] = 0.0

        normalised = normalise_age_distribution(df)

        self.assertEqual(len(normalised), 21)
        self.assertAlmostEqual(float(normalised["Proportion"].sum()), 1.0)
        self.assertEqual(normalised.iloc[-1]["AgeGroup"], "100+")

    def test_default_five_year_distribution_expands_to_age_100(self) -> None:
        expanded = expand_five_year_age_distribution(
            default_five_year_age_distribution(), population=10100
        )

        self.assertEqual(int(expanded["age"].min()), 0)
        self.assertEqual(int(expanded["age"].max()), 100)
        self.assertEqual(len(expanded), 101)
        self.assertAlmostEqual(float(expanded["population"].sum()), 10100.0)


if __name__ == "__main__":
    unittest.main()
