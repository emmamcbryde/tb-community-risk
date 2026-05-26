from __future__ import annotations

import json
import unittest

import numpy as np

from engine.dynamic.dynamic_model import (
    build_two_epoch_beta_diagnostics,
    build_two_epoch_beta_series,
    simulate_dynamic,
)


class BuildTwoEpochBetaSeriesTests(unittest.TestCase):
    def test_returns_historical_values_followed_by_recent_values(self) -> None:
        beta_series = build_two_epoch_beta_series(
            beta_historical=0.25,
            beta_recent=0.75,
            calib_years=6,
            recent_years=2,
        )

        np.testing.assert_array_equal(
            beta_series,
            np.array([0.25, 0.25, 0.25, 0.25, 0.75, 0.75], dtype=float),
        )
        self.assertEqual(beta_series.dtype, float)

    def test_length_and_change_point_use_recent_years(self) -> None:
        beta_series = build_two_epoch_beta_series(1, 2, calib_years=12, recent_years=5)

        self.assertEqual(len(beta_series), 12)
        self.assertTrue(np.all(beta_series[:7] == 1.0))
        self.assertTrue(np.all(beta_series[7:] == 2.0))

    def test_recent_years_can_cover_all_calibration_years(self) -> None:
        beta_series = build_two_epoch_beta_series(1, 2, calib_years=3, recent_years=3)

        np.testing.assert_array_equal(beta_series, np.array([2.0, 2.0, 2.0]))

    def test_invalid_inputs_raise_clear_value_errors(self) -> None:
        invalid_cases = [
            {"calib_years": 0, "recent_years": 1, "message": "calib_years"},
            {"calib_years": -1, "recent_years": 1, "message": "calib_years"},
            {"calib_years": 10.0, "recent_years": 1, "message": "calib_years"},
            {"calib_years": 10, "recent_years": 0, "message": "recent_years"},
            {"calib_years": 10, "recent_years": -1, "message": "recent_years"},
            {"calib_years": 10, "recent_years": 1.0, "message": "recent_years"},
            {"calib_years": 5, "recent_years": 6, "message": "less than or equal"},
        ]

        for case in invalid_cases:
            with self.subTest(case=case):
                with self.assertRaisesRegex(ValueError, case["message"]):
                    build_two_epoch_beta_series(
                        0.25,
                        0.75,
                        calib_years=case["calib_years"],
                        recent_years=case["recent_years"],
                    )


class BuildTwoEpochBetaDiagnosticsTests(unittest.TestCase):
    def test_projection_start_year_sets_recent_start_year(self) -> None:
        diagnostics = build_two_epoch_beta_diagnostics(
            beta_historical=0.2,
            beta_recent=0.5,
            calibration_years=list(range(2010, 2025)),
            recent_years=10,
            projection_start_year=2025,
        )

        self.assertEqual(diagnostics["beta_recent_start_year"], 2015)
        self.assertEqual(diagnostics["beta_change_index"], 5)
        self.assertEqual(diagnostics["beta_recent_years"], list(range(2015, 2025)))

    def test_beta_series_matches_projection_start_split(self) -> None:
        diagnostics = build_two_epoch_beta_diagnostics(
            beta_historical=1,
            beta_recent=2,
            calibration_years=[2020, 2021, 2022, 2023, 2024],
            recent_years=2,
            projection_start_year=2024,
        )

        self.assertEqual(diagnostics["beta_historical_years"], [2020, 2021])
        self.assertEqual(diagnostics["beta_recent_years"], [2022, 2023, 2024])
        self.assertEqual(diagnostics["beta_series"], [1.0, 1.0, 2.0, 2.0, 2.0])
        self.assertEqual(diagnostics["beta_ratio_recent_to_historical"], 2.0)

    def test_without_projection_start_year_uses_last_recent_years(self) -> None:
        diagnostics = build_two_epoch_beta_diagnostics(
            beta_historical=0.25,
            beta_recent=0.75,
            calibration_years=[2018, 2019, 2020, 2021, 2022],
            recent_years=2,
        )

        self.assertEqual(diagnostics["beta_change_index"], 3)
        self.assertEqual(diagnostics["beta_historical_years"], [2018, 2019, 2020])
        self.assertEqual(diagnostics["beta_recent_years"], [2021, 2022])
        self.assertEqual(
            diagnostics["beta_series"], [0.25, 0.25, 0.25, 0.75, 0.75]
        )

    def test_diagnostics_are_json_serialisable(self) -> None:
        diagnostics = build_two_epoch_beta_diagnostics(
            beta_historical=np.float64(0.4),
            beta_recent=np.float64(0.8),
            calibration_years=np.arange(2020, 2023),
            recent_years=1,
        )

        json.dumps(diagnostics)


class SimulateDynamicBetaSeriesTests(unittest.TestCase):
    def test_two_epoch_beta_series_runs_and_returns_requested_incidence_length(
        self,
    ) -> None:
        years = 2
        beta_series = build_two_epoch_beta_series(
            beta_historical=0.1,
            beta_recent=0.2,
            calib_years=years,
            recent_years=1,
        )
        params = {
            "age_counts": {"adult": 1000.0},
            "ltbi_ever": {"adult": 0.10},
            "ltbi_recent": {"adult": 0.02},
            "smoker_pct": 0.0,
            "alcohol_pct": 0.0,
            "diabetes_pct": 0.0,
            "renal_pct": 0.0,
            "HIV_treated_pct": 0.0,
            "HIV_untreated_pct": 0.0,
            "initial_incidence_per_100k": 0.0,
            "beta_series": beta_series,
        }

        result = simulate_dynamic(params, years=years, intervention=False)

        self.assertIn("annual_incidence", result)
        self.assertEqual(len(result["annual_incidence"]), years)


if __name__ == "__main__":
    unittest.main()
