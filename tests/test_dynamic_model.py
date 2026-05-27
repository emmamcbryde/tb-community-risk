from __future__ import annotations

import json
import unittest
from unittest.mock import patch

import numpy as np

from engine.dynamic.dynamic_model import (
    build_secular_multiplier_series,
    build_two_epoch_beta_diagnostics,
    build_two_epoch_beta_series,
    compute_annual_secular_decline_multiplier,
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
        self.assertEqual(diagnostics["beta_recent_end_year"], 2024)
        self.assertEqual(diagnostics["beta_change_year"], 2015)
        self.assertEqual(diagnostics["beta_historical_index_start"], 0)
        self.assertEqual(diagnostics["beta_historical_index_end"], 4)
        self.assertEqual(diagnostics["beta_recent_index_start"], 5)
        self.assertEqual(diagnostics["beta_recent_index_end"], 14)
        self.assertEqual(diagnostics["projection_start_year"], 2025)

    def test_projection_start_year_2025_recent_10_uses_2015_to_2024_recent_beta(
        self,
    ) -> None:
        diagnostics = build_two_epoch_beta_diagnostics(
            beta_historical=0.2,
            beta_recent=0.5,
            calibration_years=list(range(2015, 2025)),
            recent_years=10,
            projection_start_year=2025,
        )

        self.assertEqual(diagnostics["beta_recent_start_year"], 2015)
        self.assertEqual(diagnostics["beta_recent_end_year"], 2024)
        self.assertEqual(diagnostics["beta_change_index"], 0)
        self.assertEqual(diagnostics["beta_change_year"], 2015)
        self.assertEqual(diagnostics["beta_historical_years"], [])
        self.assertEqual(diagnostics["beta_recent_years"], list(range(2015, 2025)))
        self.assertEqual(diagnostics["beta_historical_index_start"], None)
        self.assertEqual(diagnostics["beta_historical_index_end"], None)
        self.assertEqual(diagnostics["beta_recent_index_start"], 0)
        self.assertEqual(diagnostics["beta_recent_index_end"], 9)
        self.assertEqual(diagnostics["beta_series"], [0.5] * 10)

    def test_beta_series_matches_projection_start_split(self) -> None:
        diagnostics = build_two_epoch_beta_diagnostics(
            beta_historical=1,
            beta_recent=2,
            calibration_years=[2020, 2021, 2022, 2023, 2024],
            recent_years=2,
            projection_start_year=2024,
            beta_bounds=(0.01, 50.0),
        )

        self.assertEqual(diagnostics["beta_historical_years"], [2020, 2021])
        self.assertEqual(diagnostics["beta_recent_years"], [2022, 2023, 2024])
        self.assertEqual(diagnostics["beta_series"], [1.0, 1.0, 2.0, 2.0, 2.0])
        self.assertEqual(diagnostics["beta1"], diagnostics["beta_historical"])
        self.assertEqual(diagnostics["beta2"], diagnostics["beta_recent"])
        self.assertEqual(diagnostics["beta_ratio_recent_to_historical"], 2.0)
        self.assertEqual(diagnostics["beta_lower_bound"], 0.01)
        self.assertEqual(diagnostics["beta_upper_bound"], 50.0)

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
            beta_bounds=(np.float64(0.01), np.float64(50.0)),
        )

        json.dumps(diagnostics, allow_nan=False)

    def test_bound_hit_flags_use_beta_bounds(self) -> None:
        diagnostics = build_two_epoch_beta_diagnostics(
            beta_historical=np.float64(0.01),
            beta_recent=np.float64(50.0),
            calibration_years=[2021, 2022],
            recent_years=1,
            beta_bounds=(0.01, 50.0),
        )

        self.assertIs(diagnostics["beta_historical_lower_bound_hit"], True)
        self.assertIs(diagnostics["beta_recent_lower_bound_hit"], False)
        self.assertIs(diagnostics["beta_historical_upper_bound_hit"], False)
        self.assertIs(diagnostics["beta_recent_upper_bound_hit"], True)


class SecularDeclineMultiplierTests(unittest.TestCase):
    def test_multiplier_is_one_before_start_and_declines_after_start(self) -> None:
        self.assertEqual(
            compute_annual_secular_decline_multiplier(
                2024,
                secular_decline_rate_annual=0.02,
                secular_decline_start_year=2025,
            ),
            1.0,
        )

        self.assertEqual(
            compute_annual_secular_decline_multiplier(
                2025,
                secular_decline_rate_annual=0.02,
                secular_decline_start_year=2025,
            ),
            1.0,
        )
        self.assertAlmostEqual(
            compute_annual_secular_decline_multiplier(
                2027,
                secular_decline_rate_annual=0.02,
                secular_decline_start_year=2025,
            ),
            np.exp(-0.02 * 2),
        )

    def test_hold_policy_holds_fixed_after_projection_start_year(self) -> None:
        self.assertAlmostEqual(
            compute_annual_secular_decline_multiplier(
                2028,
                secular_decline_rate_annual=0.02,
                secular_decline_start_year=2020,
                future_secular_trend_policy="hold",
                projection_start_year=2025,
            ),
            np.exp(-0.02 * 5),
        )

    def test_continue_capped_uses_cap_after_projection_start_year(self) -> None:
        self.assertAlmostEqual(
            compute_annual_secular_decline_multiplier(
                2027,
                secular_decline_rate_annual=0.05,
                secular_decline_start_year=2020,
                future_secular_trend_policy="continue_capped",
                future_secular_decline_cap_annual=0.02,
                projection_start_year=2025,
            ),
            np.exp(-0.05 * 5) * np.exp(-0.02 * 2),
        )

    def test_default_policy_holds_at_projection_start_year(self) -> None:
        self.assertAlmostEqual(
            compute_annual_secular_decline_multiplier(
                2028,
                secular_decline_rate_annual=0.02,
                secular_decline_start_year=2020,
                projection_start_year=2025,
            ),
            np.exp(-0.02 * 5),
        )

    def test_default_hold_uses_relative_boundary_when_projection_year_absent(
        self,
    ) -> None:
        boundary = 1
        expected = np.exp(-0.02 * (boundary - -9))
        self.assertAlmostEqual(
            compute_annual_secular_decline_multiplier(
                2,
                secular_decline_rate_annual=0.02,
                secular_decline_start_year=-9,
                future_secular_trend_policy="hold",
                projection_start_year=boundary,
            ),
            expected,
        )

    def test_continue_capped_defaults_future_cap_to_one_percent(self) -> None:
        self.assertAlmostEqual(
            compute_annual_secular_decline_multiplier(
                2027,
                secular_decline_rate_annual=0.05,
                secular_decline_start_year=2020,
                future_secular_trend_policy="continue_capped",
                projection_start_year=2025,
            ),
            np.exp(-0.05 * 5) * np.exp(-0.01 * 2),
        )


class BuildAnnualSecularMultiplierSeriesTests(unittest.TestCase):
    def test_simulation_start_year_maps_10_year_calibration_window(self) -> None:
        params = {
            "secular_decline_rate_annual": 0.02,
            "secular_decline_start_year": 2015,
            "future_secular_trend_policy": "hold",
            "projection_start_year": 2025,
            "simulation_start_year": 2015,
        }

        series = build_secular_multiplier_series(params, years=10)

        expected = [1.0] + [
            float(np.exp(-0.02 * elapsed)) for elapsed in range(1, 10)
        ]
        self.assertEqual(len(series), 10)
        np.testing.assert_allclose(series, expected)

    def test_relative_projection_boundary_holds_after_year_zero(self) -> None:
        params = {
            "secular_decline_rate_annual": 0.02,
            "secular_decline_start_year": -9,
            "future_secular_trend_policy": "hold",
            "projection_start_year": None,
            "secular_projection_boundary_year": 1,
            "simulation_start_year": 1,
        }

        series = build_secular_multiplier_series(params, years=3)

        expected = [float(np.exp(-0.02 * 10))] * 3
        np.testing.assert_allclose(series, expected)

    def test_series_is_json_serialisable(self) -> None:
        series = build_secular_multiplier_series(
            {
                "secular_decline_rate_annual": np.float64(0.02),
                "secular_decline_start_year": np.int64(2015),
                "future_secular_trend_policy": "continue_capped",
                "future_secular_decline_cap_annual": np.float64(0.01),
                "projection_start_year": np.int64(2025),
                "calendar_start_year": np.int64(2015),
            },
            years=10,
        )

        json.dumps(series, allow_nan=False)
        self.assertTrue(all(type(value) is float for value in series))


class SimulateDynamicBetaSeriesTests(unittest.TestCase):
    def _base_params(self) -> dict[str, object]:
        return {
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
        }

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
        params = self._base_params()
        params["beta_series"] = beta_series

        result = simulate_dynamic(params, years=years, intervention=False)

        self.assertIn("annual_incidence", result)
        self.assertEqual(len(result["annual_incidence"]), years)

    def test_zero_rate_reproduces_no_secular_decline_dynamic_output(self) -> None:
        years = 4
        base_params = self._base_params()
        base_params["beta"] = 0.2
        secular_params = dict(base_params)
        secular_params.update(
            {
                "secular_decline_rate_annual": 0.0,
                "secular_decline_start_year": 2020,
                "future_secular_trend_policy": "continue_capped",
                "future_secular_decline_cap_annual": 0.01,
                "projection_start_year": 2025,
            }
        )

        baseline = simulate_dynamic(base_params, years=years, intervention=False)
        with_zero_rate = simulate_dynamic(
            secular_params, years=years, intervention=False
        )

        trajectory_keys = [
            "annual_incidence",
            "annual_prevalence_I",
            "S",
            "L_fast",
            "L_slow",
            "I",
            "R",
        ]
        for key in trajectory_keys:
            np.testing.assert_allclose(with_zero_rate[key], baseline[key])
        self.assertEqual(with_zero_rate["final_state"], baseline["final_state"])

    def test_simulation_start_year_maps_backcast_secular_decline_to_historical_years(
        self,
    ) -> None:
        years = 10
        projection_start_year = 2025
        params = self._base_params()
        params.update(
            {
                "beta": 0.0,
                "secular_decline_rate_annual": 0.02,
                "secular_decline_start_year": 2015,
                "future_secular_trend_policy": "hold",
                "projection_start_year": projection_start_year,
                "simulation_start_year": projection_start_year - years,
            }
        )
        secular_years = []
        projection_boundaries = []

        def capture_multiplier(year, **kwargs):
            secular_years.append(float(year))
            projection_boundaries.append(kwargs["projection_start_year"])
            return 1.0

        with patch(
            "engine.dynamic.dynamic_model.compute_annual_secular_decline_multiplier",
            side_effect=capture_multiplier,
        ):
            simulate_dynamic(params, years=years, intervention=False)

        observed_calendar_years = sorted(
            {int(np.floor(year)) for year in secular_years}
        )
        self.assertEqual(
            observed_calendar_years,
            list(range(projection_start_year - years, projection_start_year)),
        )
        self.assertEqual(set(projection_boundaries), {projection_start_year})

    def test_absent_simulation_start_year_preserves_projection_start_year_mapping(
        self,
    ) -> None:
        years = 3
        projection_start_year = 2025
        params = self._base_params()
        params.update(
            {
                "beta": 0.0,
                "secular_decline_rate_annual": 0.02,
                "secular_decline_start_year": 2015,
                "projection_start_year": projection_start_year,
            }
        )
        secular_years = []

        def capture_multiplier(year, **kwargs):
            secular_years.append(float(year))
            return 1.0

        with patch(
            "engine.dynamic.dynamic_model.compute_annual_secular_decline_multiplier",
            side_effect=capture_multiplier,
        ):
            simulate_dynamic(params, years=years, intervention=False)

        observed_calendar_years = sorted(
            {int(np.floor(year)) for year in secular_years}
        )
        self.assertEqual(
            observed_calendar_years,
            list(range(projection_start_year, projection_start_year + years)),
        )

    def test_relative_projection_boundary_passed_to_secular_policy(self) -> None:
        years = 3
        params = self._base_params()
        params.update(
            {
                "beta": 0.0,
                "secular_decline_rate_annual": 0.02,
                "secular_decline_start_year": -9,
                "future_secular_trend_policy": "hold",
                "projection_start_year": None,
                "secular_projection_boundary_year": 1,
                "simulation_start_year": 1,
            }
        )
        secular_years = []
        projection_boundaries = []

        def capture_multiplier(year, **kwargs):
            secular_years.append(float(year))
            projection_boundaries.append(kwargs["projection_start_year"])
            return 1.0

        with patch(
            "engine.dynamic.dynamic_model.compute_annual_secular_decline_multiplier",
            side_effect=capture_multiplier,
        ):
            simulate_dynamic(params, years=years, intervention=False)

        observed_calendar_years = sorted(
            {int(np.floor(year)) for year in secular_years}
        )
        self.assertEqual(observed_calendar_years, [1, 2, 3])
        self.assertEqual(set(projection_boundaries), {1})


if __name__ == "__main__":
    unittest.main()
