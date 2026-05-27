from __future__ import annotations

import sys
import types
import unittest
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

fake_streamlit = types.ModuleType("streamlit")
fake_streamlit.session_state = {}


def _cache_resource(*args, **kwargs):
    def decorator(func):
        return func

    return decorator


fake_streamlit.cache_resource = _cache_resource
sys.modules.setdefault("streamlit", fake_streamlit)
sys.modules.setdefault("altair", types.ModuleType("altair"))

from ui.dynamic_ui import (
    CAL_MODE_RANDOM_WALK,
    CAL_MODE_TWO_EPOCH,
    _scalar_metadata,
    _secular_projection_run_years,
    calibrate_beta_two_epoch,
    calibrate_beta_two_epoch_with_secular_decline,
    calibration_success_text,
    two_epoch_diagnostic_warnings,
    two_epoch_diagnostics_display_tables,
    two_epoch_lower_bound_warning,
)


class CalibrateBetaTwoEpochTests(unittest.TestCase):
    def test_scalar_metadata_keeps_two_epoch_diagnostics_json_safe(self) -> None:
        metadata = {
            "beta_historical": np.float64(2.5),
            "beta_recent": np.float64(1.25),
            "beta_forward": np.float64(1.25),
            "beta_series": np.array([2.5, 2.5, 1.25, 1.25], dtype=float),
            "calibration_years": np.array([2021, 2022, 2023, 2024], dtype=int),
            "beta_historical_years": [2021, np.int64(2022)],
            "beta_recent_years": [np.int64(2023), 2024],
            "fitted_incidence": np.array([20.0, 18.0, 16.0, 14.0], dtype=float),
            "target_incidence": [20.0, np.float64(18.0), 15.0, 14.0],
            "residuals": np.array([0.0, 0.0, 1.0, 0.0], dtype=float),
            "solver_trace": np.array([99.0], dtype=float),
            "recent_years": np.int64(2),
            "beta_historical_lower_bound_hit": np.bool_(False),
        }

        scalar = _scalar_metadata(metadata)

        self.assertEqual(scalar["beta_series"], [2.5, 2.5, 1.25, 1.25])
        self.assertEqual(scalar["calibration_years"], [2021, 2022, 2023, 2024])
        self.assertEqual(scalar["beta_historical_years"], [2021, 2022])
        self.assertEqual(scalar["beta_recent_years"], [2023, 2024])
        self.assertEqual(scalar["fitted_incidence"], [20.0, 18.0, 16.0, 14.0])
        self.assertEqual(scalar["target_incidence"], [20.0, 18.0, 15.0, 14.0])
        self.assertEqual(scalar["residuals"], [0.0, 0.0, 1.0, 0.0])
        self.assertNotIn("solver_trace", scalar)
        self.assertEqual(scalar["beta_forward"], scalar["beta_recent"])
        self.assertEqual(scalar["beta_forward"], 1.25)
        self.assertIsInstance(scalar["beta_forward"], float)
        self.assertEqual(scalar["recent_years"], 2)
        self.assertIsInstance(scalar["recent_years"], int)
        self.assertIs(scalar["beta_historical_lower_bound_hit"], False)
        json.dumps(scalar, allow_nan=False)

    def test_two_epoch_success_text_labels_betas_and_projection_beta(self) -> None:
        text = calibration_success_text(
            CAL_MODE_TWO_EPOCH,
            beta_series_hat=np.array([2.5, 2.5, 1.25, 1.25], dtype=float),
            beta_forward=999.0,
            ari_adj_hat=1.0,
            rmse_rw=3.0,
            calibration_metadata={
                "calibration_mode": CAL_MODE_TWO_EPOCH,
                "beta_historical": 2.5,
                "beta_recent": 1.25,
            },
        )

        self.assertIn(f"calibration_mode={CAL_MODE_TWO_EPOCH}", text)
        self.assertIn("β1 / historical beta=2.50", text)
        self.assertIn("β2 / recent beta=1.25", text)
        self.assertIn("Using β2 / recent beta=1.25 for projections", text)
        self.assertNotIn("Using β=999.00 for projections", text)

    def test_non_two_epoch_success_text_keeps_generic_beta_message(self) -> None:
        text = calibration_success_text(
            CAL_MODE_RANDOM_WALK,
            beta_series_hat=np.array([1.0, 1.5, 2.0], dtype=float),
            beta_forward=2.0,
            ari_adj_hat=1.0,
            rmse_rw=3.0,
            calibration_metadata=None,
        )

        self.assertIn(f"calibration_mode={CAL_MODE_RANDOM_WALK}", text)
        self.assertIn("β(t) range 1.00–2.00", text)
        self.assertIn("Using β=2.00 for projections", text)
        self.assertNotIn("β1 / historical beta", text)
        self.assertNotIn("β2 / recent beta", text)

    def test_two_epoch_display_tables_include_requested_diagnostics(self) -> None:
        metadata = {
            "calibration_mode": "Two-epoch beta: historical + recent 10 years",
            "beta_historical": 2.0,
            "beta_recent": 1.0,
            "beta_forward": 1.0,
            "beta_ratio_recent_to_historical": 0.5,
            "beta_recent_start_year": 2023,
            "recent_years": 2,
            "beta_lower_bound": 0.01,
            "beta_upper_bound": 50.0,
            "beta_historical_lower_bound_hit": False,
            "beta_recent_lower_bound_hit": False,
            "beta_historical_upper_bound_hit": False,
            "beta_recent_upper_bound_hit": False,
            "beta_historical_near_lower_bound": False,
            "beta_recent_near_lower_bound": False,
            "beta_historical_near_upper_bound": False,
            "beta_recent_near_upper_bound": False,
            "secular_decline_rate_annual": 0.012,
            "secular_decline_rate_percent_annual": 1.2,
            "secular_decline_start_year": 2023,
            "future_secular_trend_policy": "hold",
            "future_secular_decline_cap_annual": 0.01,
            "secular_decline_rate_lower_bound": 0.0,
            "secular_decline_rate_upper_bound": 0.10,
            "secular_decline_rate_lower_bound_hit": False,
            "secular_decline_rate_upper_bound_hit": False,
            "rmse_overall": 3.0,
            "rmse_historical": 4.0,
            "rmse_recent": 1.5,
            "beta_change_index": 3,
            "beta_change_year": 2023,
            "beta_historical_index_start": 0,
            "beta_historical_index_end": 2,
            "beta_recent_index_start": 3,
            "beta_recent_index_end": 4,
            "beta_historical_years": [2020, 2021, 2022],
            "beta_recent_years": [2023, 2024],
        }

        summary_df, preview_df = two_epoch_diagnostics_display_tables(metadata)

        summary = dict(zip(summary_df["Metric"], summary_df["Value"]))
        metric_text = "\n".join(summary_df["Metric"].tolist())
        for field_name in (
            "beta_historical",
            "beta_recent",
            "beta_recent_start_year",
            "recent_years",
            "beta_lower_bound",
            "beta_upper_bound",
            "beta_forward",
            "projection beta",
            "secular_decline_rate_percent_annual",
            "background calibration trend",
            "secular_decline_start_year",
            "future_secular_trend_policy",
            "future_secular_decline_cap_annual",
            "secular_decline_rate_lower_bound",
            "secular_decline_rate_upper_bound",
            "secular_decline_rate_lower_bound_hit",
            "secular_decline_rate_upper_bound_hit",
            "beta_historical_near_lower_bound",
            "beta_recent_near_lower_bound",
            "beta_historical_near_upper_bound",
            "beta_recent_near_upper_bound",
        ):
            self.assertIn(field_name, metric_text)

        self.assertEqual(
            summary["beta_historical (Beta 1 / historical beta)"], 2.0
        )
        self.assertEqual(summary["beta_recent (Beta 2 / recent beta)"], 1.0)
        self.assertEqual(summary["beta_forward / projection beta"], 1.0)
        self.assertEqual(summary["Beta 2 / Beta 1 ratio"], 0.5)
        self.assertEqual(
            summary["beta_recent_start_year (recent beta start year)"], 2023
        )
        self.assertEqual(summary["recent_years (recent window length)"], 2)
        self.assertEqual(summary["beta_lower_bound"], 0.01)
        self.assertEqual(summary["beta_upper_bound"], 50.0)
        self.assertIs(summary["beta_historical_lower_bound_hit"], False)
        self.assertIs(summary["beta_recent_lower_bound_hit"], False)
        self.assertIs(summary["beta_historical_upper_bound_hit"], False)
        self.assertIs(summary["beta_recent_upper_bound_hit"], False)
        self.assertIs(summary["beta_historical_near_lower_bound"], False)
        self.assertIs(summary["beta_recent_near_lower_bound"], False)
        self.assertIs(summary["beta_historical_near_upper_bound"], False)
        self.assertIs(summary["beta_recent_near_upper_bound"], False)
        self.assertEqual(
            summary[
                "secular_decline_rate_percent_annual (background calibration trend, %/year)"
            ],
            1.2,
        )
        self.assertEqual(
            summary["secular_decline_start_year (background trend start year)"],
            2023,
        )
        self.assertEqual(
            summary["future_secular_trend_policy (background trend policy)"],
            "hold",
        )
        self.assertEqual(
            summary["future_secular_decline_cap_annual (background trend cap)"],
            0.01,
        )
        self.assertEqual(summary["secular_decline_rate_lower_bound"], 0.0)
        self.assertEqual(summary["secular_decline_rate_upper_bound"], 0.10)
        self.assertIs(summary["secular_decline_rate_lower_bound_hit"], False)
        self.assertIs(summary["secular_decline_rate_upper_bound_hit"], False)
        self.assertEqual(summary["RMSE overall"], 3.0)
        self.assertEqual(summary["RMSE historical"], 4.0)
        self.assertEqual(summary["RMSE recent"], 1.5)
        self.assertEqual(summary["Beta change index"], 3)
        self.assertEqual(summary["Beta change year"], 2023)
        self.assertEqual(
            preview_df["Epoch"].tolist(), ["Beta 1 / historical", "Beta 2 / recent"]
        )
        self.assertEqual(preview_df["Years"].tolist(), ["2020-2022", "2023-2024"])
        self.assertEqual(preview_df["Start year"].tolist(), [2020, 2023])
        self.assertEqual(preview_df["Start index"].tolist(), [0, 3])
        self.assertEqual(preview_df["End index"].tolist(), [2, 4])
        self.assertEqual(preview_df["Window length"].tolist(), [3, 2])
        self.assertEqual(preview_df["Beta"].tolist(), [2.0, 1.0])
        self.assertEqual(preview_df["RMSE"].tolist(), [4.0, 1.5])
        self.assertEqual(preview_df["Lower bound hit"].tolist(), [False, False])
        self.assertEqual(preview_df["Upper bound hit"].tolist(), [False, False])
        self.assertEqual(preview_df["Near lower bound"].tolist(), [False, False])
        self.assertEqual(preview_df["Near upper bound"].tolist(), [False, False])

    def test_two_epoch_display_tables_accept_beta_fields_without_mode(self) -> None:
        metadata = {
            "beta_historical": 2.0,
            "beta_recent": 1.0,
            "beta_forward": 1.0,
            "beta_historical_years": [2020, 2021],
            "beta_recent_years": [2022, 2023],
            "beta_recent_start_year": 2022,
        }

        summary_df, preview_df = two_epoch_diagnostics_display_tables(metadata)

        summary = dict(zip(summary_df["Metric"], summary_df["Value"]))
        self.assertEqual(
            summary["beta_historical (Beta 1 / historical beta)"], 2.0
        )
        self.assertEqual(summary["beta_recent (Beta 2 / recent beta)"], 1.0)
        self.assertEqual(summary["beta_forward / projection beta"], 1.0)
        self.assertEqual(preview_df["Epoch"].tolist(), ["Beta 1 / historical", "Beta 2 / recent"])
        self.assertEqual(preview_df["Years"].tolist(), ["2020-2021", "2022-2023"])

    def test_two_epoch_lower_bound_warning_requires_both_lower_hits(self) -> None:
        metadata = {
            "calibration_mode": "Two-epoch beta: historical + recent 10 years",
            "beta_historical_lower_bound_hit": np.bool_(True),
            "beta_recent_lower_bound_hit": True,
        }

        self.assertEqual(
            two_epoch_lower_bound_warning(metadata),
            "Both beta estimates are at the lower bound. Calibration may be "
            "constrained by beta bounds rather than data fit.",
        )

        metadata["beta_recent_lower_bound_hit"] = False
        self.assertIsNone(two_epoch_lower_bound_warning(metadata))

    def test_two_epoch_lower_bound_warning_falls_back_to_beta_values(self) -> None:
        metadata = {
            "calibration_mode": CAL_MODE_TWO_EPOCH,
            "beta_historical": np.float64(0.01),
            "beta_recent": 0.01,
            "beta_lower_bound": 0.01,
        }

        self.assertEqual(
            two_epoch_lower_bound_warning(metadata),
            "Both beta estimates are at the lower bound. Calibration may be "
            "constrained by beta bounds rather than data fit.",
        )

    def test_two_epoch_lower_bound_warning_uses_default_floor_without_bound(
        self,
    ) -> None:
        metadata = {
            "calibration_mode": CAL_MODE_TWO_EPOCH,
            "beta_historical": np.float64(0.01),
            "beta_recent": 0.01,
        }

        self.assertEqual(
            two_epoch_lower_bound_warning(metadata),
            "Both beta estimates are at the lower bound. Calibration may be "
            "constrained by beta bounds rather than data fit.",
        )

    def test_two_epoch_lower_bound_warning_values_override_false_flags(self) -> None:
        metadata = {
            "calibration_mode": CAL_MODE_TWO_EPOCH,
            "beta_historical": 0.01,
            "beta_recent": np.float64(0.01),
            "beta_lower_bound": np.float64(0.01),
            "beta_historical_lower_bound_hit": False,
            "beta_recent_lower_bound_hit": np.bool_(False),
        }

        self.assertIsNotNone(two_epoch_lower_bound_warning(metadata))

        metadata["beta_recent"] = 0.02
        self.assertIsNone(two_epoch_lower_bound_warning(metadata))

    def test_two_epoch_lower_bound_warning_accepts_beta_fields_without_mode(
        self,
    ) -> None:
        metadata = {
            "beta_historical": 0.01,
            "beta_recent": np.float64(0.01),
            "beta_lower_bound": 0.01,
        }

        self.assertEqual(
            two_epoch_lower_bound_warning(metadata),
            "Both beta estimates are at the lower bound. Calibration may be "
            "constrained by beta bounds rather than data fit.",
        )

    def test_two_epoch_diagnostic_warnings_include_bound_lists_and_secular_context(
        self,
    ) -> None:
        metadata = {
            "calibration_mode": CAL_MODE_TWO_EPOCH,
            "beta_historical": 1.0,
            "beta_recent": 2.0,
            "beta_lower_bound": 1.0,
            "beta_historical_near_lower_bound": True,
            "secular_decline_rate_annual": 0.006,
            "parameter_bound_warnings": [
                "One or more beta estimates are at or near the lower bound."
            ],
            "bound_warnings": [
                "One or more beta estimates are at or near the lower bound.",
                "Secular decline rate is at or near the upper bound.",
            ],
        }

        warnings = two_epoch_diagnostic_warnings(metadata)

        self.assertEqual(
            warnings.count(
                "One or more beta estimates are at or near the lower bound."
            ),
            1,
        )
        self.assertIn("Secular decline rate is at or near the upper bound.", warnings)
        secular_warnings = [
            warning for warning in warnings if "Background calibration trend" in warning
        ]
        self.assertEqual(len(secular_warnings), 1)
        self.assertIn("not an intervention effect", secular_warnings[0])

    def test_two_epoch_diagnostic_warnings_ignore_small_secular_decline(self) -> None:
        metadata = {
            "calibration_mode": CAL_MODE_TWO_EPOCH,
            "beta_historical_near_lower_bound": True,
            "secular_decline_rate_annual": 0.005,
        }

        warnings = two_epoch_diagnostic_warnings(metadata)

        self.assertFalse(
            any("Background calibration trend" in warning for warning in warnings)
        )

    def test_two_epoch_display_tables_ignore_other_modes(self) -> None:
        self.assertIsNone(
            two_epoch_diagnostics_display_tables(
                {"calibration_mode": "Random-walk beta"}
            )
        )

    def test_falling_target_calibrates_distinct_two_epoch_beta_series(self) -> None:
        age_counts = {"all": 100000.0}
        ages = list(age_counts)
        inc_hist = {
            -4: 40.0,
            -3: 40.0,
            -2: 38.0,
            -1: 24.0,
            0: 22.0,
        }
        risk_inputs = {
            "smoker_pct": 0.0,
            "alcohol_pct": 0.0,
            "diabetes_pct": 0.0,
            "renal_pct": 0.0,
            "HIV_treated_pct": 0.0,
            "HIV_untreated_pct": 0.0,
        }

        def fake_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
            return (
                {age: 0.10 for age in ages},
                {age: 0.02 for age in ages},
                {"ari_adjustment": ari_adjustment, "shift_years": shift_years},
            )

        def fake_run_dynamic_model(params, years, intervention=True):
            beta_series = np.asarray(params["beta_series"], dtype=float)
            self.assertEqual(len(beta_series), years)
            return {"annual_incidence": beta_series}

        with patch(
            "ui.dynamic_ui.compute_ltbi_from_inc_hist",
            side_effect=fake_ltbi_from_inc_hist,
        ), patch(
            "ui.dynamic_ui.run_dynamic_model",
            side_effect=fake_run_dynamic_model,
        ):
            beta_forward, adj, rmse, obs, metadata = calibrate_beta_two_epoch(
                age_counts=age_counts,
                ages=ages,
                inc_hist=inc_hist,
                calib_years=4,
                risk_inputs=risk_inputs,
                pre_det_months=12.0,
                delta_pre=0.0,
                beta_bounds=(1.0, 60.0),
                adj_bounds=(1.0, 1.0),
                adj_grid_points=1,
                recent_years=2,
            )

        self.assertTrue(np.isfinite(rmse))
        self.assertEqual(adj, 1.0)
        np.testing.assert_array_equal(obs, np.array([40.0, 38.0, 24.0, 22.0]))

        beta_historical = metadata["beta_historical"]
        beta_recent = metadata["beta_recent"]
        beta_series = metadata["beta_series"]

        self.assertNotEqual(beta_historical, beta_recent)
        self.assertLess(beta_recent, beta_historical)
        self.assertEqual(beta_forward, beta_recent)
        self.assertEqual(metadata["beta_forward"], beta_recent)
        np.testing.assert_allclose(
            beta_series,
            np.array(
                [beta_historical, beta_historical, beta_recent, beta_recent],
                dtype=float,
            ),
        )

    def test_projection_start_year_controls_two_epoch_split(self) -> None:
        age_counts = {"all": 100000.0}
        ages = list(age_counts)
        inc_hist = {
            -5: 30.0,
            -4: 30.0,
            -3: 30.0,
            -2: 18.0,
            -1: 18.0,
            0: 18.0,
        }
        risk_inputs = {
            "smoker_pct": 0.0,
            "alcohol_pct": 0.0,
            "diabetes_pct": 0.0,
            "renal_pct": 0.0,
            "HIV_treated_pct": 0.0,
            "HIV_untreated_pct": 0.0,
        }

        def fake_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
            return (
                {age: 0.10 for age in ages},
                {age: 0.02 for age in ages},
                {"ari_adjustment": ari_adjustment, "shift_years": shift_years},
            )

        def fake_run_dynamic_model(params, years, intervention=True):
            beta_series = np.asarray(params["beta_series"], dtype=float)
            self.assertEqual(len(beta_series), years)
            return {"annual_incidence": beta_series}

        with patch(
            "ui.dynamic_ui.compute_ltbi_from_inc_hist",
            side_effect=fake_ltbi_from_inc_hist,
        ), patch(
            "ui.dynamic_ui.run_dynamic_model",
            side_effect=fake_run_dynamic_model,
        ):
            beta_forward, _, _, obs, metadata = calibrate_beta_two_epoch(
                age_counts=age_counts,
                ages=ages,
                inc_hist=inc_hist,
                calib_years=5,
                risk_inputs=risk_inputs,
                pre_det_months=12.0,
                delta_pre=0.0,
                beta_bounds=(1.0, 60.0),
                adj_bounds=(1.0, 1.0),
                adj_grid_points=1,
                recent_years=2,
                projection_start_year=2025,
            )

        self.assertEqual(metadata["projection_start_year"], 2025)
        self.assertEqual(metadata["calibration_years"], [2020, 2021, 2022, 2023, 2024])
        self.assertEqual(metadata["beta_recent_start_year"], 2023)
        self.assertEqual(metadata["beta_change_index"], 3)
        self.assertEqual(metadata["beta_historical_years"], [2020, 2021, 2022])
        self.assertEqual(metadata["beta_recent_years"], [2023, 2024])
        self.assertEqual(beta_forward, metadata["beta_recent"])
        self.assertEqual(metadata["beta_forward"], metadata["beta_recent"])
        np.testing.assert_array_equal(obs, np.array([30.0, 30.0, 18.0, 18.0, 18.0]))

    def test_projection_start_year_2025_recent_10_uses_recent_beta_for_projection(
        self,
    ) -> None:
        age_counts = {"all": 100000.0}
        ages = list(age_counts)
        inc_hist = {offset: 12.0 for offset in range(-10, 1)}
        risk_inputs = {
            "smoker_pct": 0.0,
            "alcohol_pct": 0.0,
            "diabetes_pct": 0.0,
            "renal_pct": 0.0,
            "HIV_treated_pct": 0.0,
            "HIV_untreated_pct": 0.0,
        }

        def fake_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
            return (
                {age: 0.10 for age in ages},
                {age: 0.02 for age in ages},
                {"ari_adjustment": ari_adjustment, "shift_years": shift_years},
            )

        def fake_run_dynamic_model(params, years, intervention=True):
            beta_series = np.asarray(params["beta_series"], dtype=float)
            self.assertEqual(len(beta_series), years)
            np.testing.assert_allclose(beta_series, np.full(years, beta_series[-1]))
            self.assertEqual(float(params["beta"]), float(beta_series[-1]))
            return {"annual_incidence": beta_series}

        with patch(
            "ui.dynamic_ui.compute_ltbi_from_inc_hist",
            side_effect=fake_ltbi_from_inc_hist,
        ), patch(
            "ui.dynamic_ui.run_dynamic_model",
            side_effect=fake_run_dynamic_model,
        ):
            beta_forward, _, _, obs, metadata = calibrate_beta_two_epoch(
                age_counts=age_counts,
                ages=ages,
                inc_hist=inc_hist,
                calib_years=10,
                risk_inputs=risk_inputs,
                pre_det_months=12.0,
                delta_pre=0.0,
                beta_bounds=(1.0, 60.0),
                adj_bounds=(1.0, 1.0),
                adj_grid_points=1,
                recent_years=10,
                projection_start_year=2025,
            )

        self.assertEqual(metadata["projection_start_year"], 2025)
        self.assertEqual(metadata["calibration_years"], list(range(2015, 2025)))
        self.assertEqual(
            metadata["beta_recent_start_year"],
            metadata["projection_start_year"] - 10,
        )
        self.assertEqual(metadata["beta_recent_start_year"], 2015)
        self.assertEqual(metadata["beta_recent_end_year"], 2024)
        self.assertEqual(metadata["beta_change_index"], 0)
        self.assertEqual(metadata["beta_change_year"], 2015)
        self.assertEqual(metadata["beta_historical_years"], [])
        self.assertEqual(metadata["beta_recent_years"], list(range(2015, 2025)))
        self.assertIsNone(metadata["beta_historical_index_start"])
        self.assertIsNone(metadata["beta_historical_index_end"])
        self.assertEqual(metadata["beta_recent_index_start"], 0)
        self.assertEqual(metadata["beta_recent_index_end"], 9)
        self.assertEqual(beta_forward, metadata["beta_recent"])
        self.assertEqual(metadata["beta_forward"], metadata["beta_recent"])
        np.testing.assert_allclose(
            metadata["beta_series"], np.full(10, metadata["beta_recent"])
        )
        np.testing.assert_array_equal(obs, np.full(10, 12.0))

    def test_two_epoch_metadata_has_diagnostics_and_is_json_serialisable(self) -> None:
        age_counts = {"all": 100000.0}
        ages = list(age_counts)
        inc_hist = {-3: 12.0, -2: 12.0, -1: 6.0, 0: 6.0}
        risk_inputs = {
            "smoker_pct": 0.0,
            "alcohol_pct": 0.0,
            "diabetes_pct": 0.0,
            "renal_pct": 0.0,
            "HIV_treated_pct": 0.0,
            "HIV_untreated_pct": 0.0,
        }

        def fake_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
            return (
                {age: 0.10 for age in ages},
                {age: 0.02 for age in ages},
                {"ari_adjustment": ari_adjustment, "shift_years": shift_years},
            )

        def fake_run_dynamic_model(params, years, intervention=True):
            return {"annual_incidence": np.asarray(params["beta_series"], dtype=float)}

        with patch(
            "ui.dynamic_ui.compute_ltbi_from_inc_hist",
            side_effect=fake_ltbi_from_inc_hist,
        ), patch(
            "ui.dynamic_ui.run_dynamic_model",
            side_effect=fake_run_dynamic_model,
        ):
            _, _, _, _, metadata = calibrate_beta_two_epoch(
                age_counts=age_counts,
                ages=ages,
                inc_hist=inc_hist,
                calib_years=4,
                risk_inputs=risk_inputs,
                pre_det_months=12.0,
                delta_pre=0.0,
                beta_bounds=(1.0, 20.0),
                adj_bounds=(1.0, 1.0),
                adj_grid_points=1,
                recent_years=2,
            )

        for key in [
            "calibration_mode",
            "beta1",
            "beta2",
            "beta_historical",
            "beta_recent",
            "beta_lower_bound",
            "beta_upper_bound",
            "beta_historical_lower_bound_hit",
            "beta_recent_lower_bound_hit",
            "beta_historical_upper_bound_hit",
            "beta_recent_upper_bound_hit",
            "beta_ratio_recent_to_historical",
            "beta_ratio",
            "recent_years",
            "projection_start_year",
            "beta_recent_start_year",
            "beta_recent_end_year",
            "beta_change_year",
            "beta_historical_start_year",
            "beta_historical_end_year",
            "beta_historical_index_start",
            "beta_historical_index_end",
            "beta_recent_index_start",
            "beta_recent_index_end",
            "beta_historical_years",
            "beta_recent_years",
            "beta_change_index",
            "beta_series",
            "calibration_years",
            "fitted_incidence",
            "target_incidence",
            "residuals",
            "rmse_overall",
            "rmse_historical",
            "rmse_recent",
            "beta_forward",
        ]:
            self.assertIn(key, metadata)

        self.assertEqual(
            metadata["calibration_mode"],
            "Two-epoch beta: historical + recent 10 years",
        )
        self.assertEqual(metadata["beta1"], metadata["beta_historical"])
        self.assertEqual(metadata["beta2"], metadata["beta_recent"])
        self.assertEqual(metadata["beta_lower_bound"], 1.0)
        self.assertEqual(metadata["beta_upper_bound"], 20.0)
        self.assertFalse(metadata["beta_historical_lower_bound_hit"])
        self.assertFalse(metadata["beta_recent_lower_bound_hit"])
        self.assertFalse(metadata["beta_historical_upper_bound_hit"])
        self.assertFalse(metadata["beta_recent_upper_bound_hit"])
        self.assertEqual(
            metadata["beta_ratio"],
            metadata["beta_ratio_recent_to_historical"],
        )
        self.assertEqual(metadata["beta_forward"], metadata["beta_recent"])
        self.assertEqual(len(metadata["fitted_incidence"]), 4)
        self.assertEqual(len(metadata["target_incidence"]), 4)
        self.assertEqual(len(metadata["residuals"]), 4)
        json.dumps(metadata, allow_nan=False)

    def test_two_epoch_metadata_flags_beta_bounds(self) -> None:
        age_counts = {"all": 100000.0}
        ages = list(age_counts)
        inc_hist = {-3: 1.0, -2: 1.0, -1: 20.0, 0: 20.0}
        risk_inputs = {
            "smoker_pct": 0.0,
            "alcohol_pct": 0.0,
            "diabetes_pct": 0.0,
            "renal_pct": 0.0,
            "HIV_treated_pct": 0.0,
            "HIV_untreated_pct": 0.0,
        }

        def fake_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
            return (
                {age: 0.10 for age in ages},
                {age: 0.02 for age in ages},
                {"ari_adjustment": ari_adjustment, "shift_years": shift_years},
            )

        def fake_run_dynamic_model(params, years, intervention=True):
            return {"annual_incidence": np.asarray(params["beta_series"], dtype=float)}

        with patch(
            "ui.dynamic_ui.compute_ltbi_from_inc_hist",
            side_effect=fake_ltbi_from_inc_hist,
        ), patch(
            "ui.dynamic_ui.run_dynamic_model",
            side_effect=fake_run_dynamic_model,
        ):
            _, _, _, _, metadata = calibrate_beta_two_epoch(
                age_counts=age_counts,
                ages=ages,
                inc_hist=inc_hist,
                calib_years=4,
                risk_inputs=risk_inputs,
                pre_det_months=12.0,
                delta_pre=0.0,
                beta_bounds=(1.0, 20.0),
                adj_bounds=(1.0, 1.0),
                adj_grid_points=1,
                recent_years=2,
            )

        self.assertTrue(
            np.isclose(metadata["beta_historical"], metadata["beta_lower_bound"])
        )
        self.assertTrue(
            np.isclose(metadata["beta_recent"], metadata["beta_upper_bound"])
        )
        self.assertIs(metadata["beta_historical_lower_bound_hit"], True)
        self.assertIs(metadata["beta_recent_lower_bound_hit"], False)
        self.assertIs(metadata["beta_historical_upper_bound_hit"], False)
        self.assertIs(metadata["beta_recent_upper_bound_hit"], True)
        json.dumps(metadata, allow_nan=False)

    def test_two_epoch_secular_metadata_is_json_serialisable(self) -> None:
        age_counts = {"all": 100000.0}
        ages = list(age_counts)
        inc_hist = {-3: 12.0, -2: 12.0, -1: 6.0, 0: 6.0}
        risk_inputs = {
            "smoker_pct": 0.0,
            "alcohol_pct": 0.0,
            "diabetes_pct": 0.0,
            "renal_pct": 0.0,
            "HIV_treated_pct": 0.0,
            "HIV_untreated_pct": 0.0,
        }

        def fake_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
            return (
                {age: 0.10 for age in ages},
                {age: 0.02 for age in ages},
                {"ari_adjustment": ari_adjustment, "shift_years": shift_years},
            )

        def fake_run_dynamic_model(params, years, intervention=True):
            beta_series = np.asarray(params["beta_series"], dtype=float)
            rate = float(params.get("secular_decline_rate_annual", 0.0))
            start_year = int(params["secular_decline_start_year"])
            sim_start = int(params["simulation_start_year"])
            multipliers = np.array(
                [
                    np.exp(-rate * max(0, sim_start + idx - start_year))
                    for idx in range(years)
                ],
                dtype=float,
            )
            return {"annual_incidence": beta_series * multipliers}

        with patch(
            "ui.dynamic_ui.compute_ltbi_from_inc_hist",
            side_effect=fake_ltbi_from_inc_hist,
        ), patch(
            "ui.dynamic_ui.run_dynamic_model",
            side_effect=fake_run_dynamic_model,
        ):
            _, _, _, _, metadata = calibrate_beta_two_epoch_with_secular_decline(
                age_counts=age_counts,
                ages=ages,
                inc_hist=inc_hist,
                calib_years=4,
                risk_inputs=risk_inputs,
                pre_det_months=12.0,
                delta_pre=0.0,
                beta_bounds=(1.0, 20.0),
                adj_bounds=(1.0, 1.0),
                adj_grid_points=1,
                recent_years=2,
                projection_start_year=2025,
            )

        for key in [
            "calibration_mode",
            "beta_historical",
            "beta_recent",
            "beta_ratio",
            "beta_lower_bound",
            "beta_upper_bound",
            "beta_historical_lower_bound_hit",
            "beta_recent_lower_bound_hit",
            "beta_ratio_recent_to_historical",
            "secular_decline_rate_annual",
            "secular_decline_rate_percent_annual",
            "secular_decline_percent_annual",
            "secular_decline_start_year",
            "recent_years",
            "future_secular_trend_policy",
            "future_secular_decline_cap_annual",
            "secular_multiplier_series",
            "beta_series",
            "calibration_years",
            "fitted_incidence",
            "target_incidence",
            "residuals",
            "rmse_overall",
            "rmse_historical",
            "rmse_recent",
            "bound_warnings",
            "parameter_bound_warnings",
            "messages",
        ]:
            self.assertIn(key, metadata)

        self.assertEqual(
            metadata["calibration_mode"],
            "Two-epoch beta with secular decline: historical + recent 2 years",
        )
        self.assertEqual(metadata["projection_start_year"], 2025)
        self.assertEqual(metadata["simulation_start_year"], 2021)
        self.assertEqual(metadata["secular_projection_boundary_year"], 2025)
        self.assertEqual(metadata["secular_decline_start_year"], 2023)
        self.assertEqual(metadata["future_secular_trend_policy"], "hold")
        self.assertEqual(metadata["future_secular_decline_cap_annual"], 0.01)
        self.assertEqual(
            metadata["secular_decline_percent_annual"],
            metadata["secular_decline_rate_percent_annual"],
        )
        self.assertIs(metadata["parameter_bound_warnings"], metadata["bound_warnings"])
        self.assertEqual(len(metadata["secular_multiplier_series"]), 4)
        json.dumps(metadata, allow_nan=False)

    def test_two_epoch_secular_defaults_relative_recent_start_without_projection_year(
        self,
    ) -> None:
        age_counts = {"all": 100000.0}
        ages = list(age_counts)
        inc_hist = {offset: 12.0 for offset in range(-11, 1)}
        risk_inputs = {
            "smoker_pct": 0.0,
            "alcohol_pct": 0.0,
            "diabetes_pct": 0.0,
            "renal_pct": 0.0,
            "HIV_treated_pct": 0.0,
            "HIV_untreated_pct": 0.0,
        }

        def fake_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
            return (
                {age: 0.10 for age in ages},
                {age: 0.02 for age in ages},
                {"ari_adjustment": ari_adjustment, "shift_years": shift_years},
            )

        seen_model_params = []

        def fake_run_dynamic_model(params, years, intervention=True):
            seen_model_params.append(dict(params))
            beta_series = np.asarray(params["beta_series"], dtype=float)
            rate = float(params.get("secular_decline_rate_annual", 0.0))
            start_year = int(params["secular_decline_start_year"])
            sim_start = int(params["simulation_start_year"])
            multipliers = np.array(
                [
                    np.exp(-rate * max(0, sim_start + idx - start_year))
                    for idx in range(years)
                ],
                dtype=float,
            )
            return {"annual_incidence": beta_series * multipliers}

        with patch(
            "ui.dynamic_ui.compute_ltbi_from_inc_hist",
            side_effect=fake_ltbi_from_inc_hist,
        ), patch(
            "ui.dynamic_ui.run_dynamic_model",
            side_effect=fake_run_dynamic_model,
        ):
            _, _, _, _, metadata = calibrate_beta_two_epoch_with_secular_decline(
                age_counts=age_counts,
                ages=ages,
                inc_hist=inc_hist,
                calib_years=12,
                risk_inputs=risk_inputs,
                pre_det_months=12.0,
                delta_pre=0.0,
                beta_bounds=(1.0, 1.0),
                adj_bounds=(1.0, 1.0),
                adj_grid_points=1,
                recent_years=10,
                projection_start_year=None,
                secular_decline_rate_bounds=(0.02, 0.02),
            )

        self.assertIsNone(metadata["projection_start_year"])
        self.assertEqual(metadata["simulation_start_year"], -11)
        self.assertEqual(metadata["secular_projection_boundary_year"], 1)
        self.assertTrue(seen_model_params)
        self.assertTrue(
            all(params["projection_start_year"] is None for params in seen_model_params)
        )
        self.assertEqual(
            {params["secular_projection_boundary_year"] for params in seen_model_params},
            {1},
        )
        self.assertEqual(metadata["calibration_years"], list(range(-11, 1)))
        self.assertEqual(metadata["beta_recent_years"], list(range(-9, 1)))
        self.assertEqual(metadata["secular_decline_start_year"], -9)
        self.assertEqual(len(metadata["secular_multiplier_series"]), 12)
        np.testing.assert_allclose(
            metadata["secular_multiplier_series"][:3],
            np.ones(3),
        )
        self.assertLess(metadata["secular_multiplier_series"][3], 1.0)
        self.assertLess(
            metadata["secular_multiplier_series"][-1],
            metadata["secular_multiplier_series"][3],
        )

    def test_secular_projection_run_years_use_relative_boundary_without_projection_year(
        self,
    ) -> None:
        projection_year, boundary_year, simulation_start_year = (
            _secular_projection_run_years(
                {
                    "projection_start_year": None,
                    "secular_projection_boundary_year": 1,
                }
            )
        )

        self.assertIsNone(projection_year)
        self.assertEqual(boundary_year, 1)
        self.assertEqual(simulation_start_year, 1)

    def test_zero_secular_rate_reproduces_two_epoch_constant_target(self) -> None:
        age_counts = {"all": 100000.0}
        ages = list(age_counts)
        inc_hist = {offset: 12.0 for offset in range(-4, 1)}
        risk_inputs = {
            "smoker_pct": 0.0,
            "alcohol_pct": 0.0,
            "diabetes_pct": 0.0,
            "renal_pct": 0.0,
            "HIV_treated_pct": 0.0,
            "HIV_untreated_pct": 0.0,
        }

        def fake_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
            return (
                {age: 0.10 for age in ages},
                {age: 0.02 for age in ages},
                {"ari_adjustment": ari_adjustment, "shift_years": shift_years},
            )

        def fake_run_dynamic_model(params, years, intervention=True):
            self.assertEqual(params["simulation_start_year"], 2021)
            return {"annual_incidence": np.asarray(params["beta_series"], dtype=float)}

        with patch(
            "ui.dynamic_ui.compute_ltbi_from_inc_hist",
            side_effect=fake_ltbi_from_inc_hist,
        ), patch(
            "ui.dynamic_ui.run_dynamic_model",
            side_effect=fake_run_dynamic_model,
        ):
            beta_forward, adj, rmse, obs, metadata = (
                calibrate_beta_two_epoch_with_secular_decline(
                    age_counts=age_counts,
                    ages=ages,
                    inc_hist=inc_hist,
                    calib_years=4,
                    risk_inputs=risk_inputs,
                    pre_det_months=12.0,
                    delta_pre=0.0,
                    beta_bounds=(1.0, 20.0),
                    adj_bounds=(1.0, 1.0),
                    adj_grid_points=1,
                    recent_years=2,
                    projection_start_year=2025,
                    secular_decline_rate_bounds=(0.0, 0.0),
                )
            )

        self.assertEqual(adj, 1.0)
        self.assertLess(rmse, 1e-5)
        self.assertEqual(beta_forward, metadata["beta_recent"])
        self.assertEqual(metadata["secular_decline_rate_annual"], 0.0)
        np.testing.assert_allclose(metadata["secular_multiplier_series"], np.ones(4))
        np.testing.assert_allclose(metadata["beta_series"], np.full(4, 12.0), atol=1e-5)
        np.testing.assert_array_equal(obs, np.full(4, 12.0))

    def test_falling_target_estimates_positive_secular_decline_when_beta_lower_bound(
        self,
    ) -> None:
        age_counts = {"all": 100000.0}
        ages = list(age_counts)
        inc_hist = {-4: 10.0, -3: 10.0, -2: 8.0, -1: 6.5, 0: 5.5}
        risk_inputs = {
            "smoker_pct": 0.0,
            "alcohol_pct": 0.0,
            "diabetes_pct": 0.0,
            "renal_pct": 0.0,
            "HIV_treated_pct": 0.0,
            "HIV_untreated_pct": 0.0,
        }

        def fake_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
            return (
                {age: 0.10 for age in ages},
                {age: 0.02 for age in ages},
                {"ari_adjustment": ari_adjustment, "shift_years": shift_years},
            )

        def fake_run_dynamic_model(params, years, intervention=True):
            beta_series = np.asarray(params["beta_series"], dtype=float)
            rate = float(params.get("secular_decline_rate_annual", 0.0))
            start_year = int(params["secular_decline_start_year"])
            sim_start = int(params["simulation_start_year"])
            multipliers = np.array(
                [
                    np.exp(-rate * max(0, sim_start + idx - start_year))
                    for idx in range(years)
                ],
                dtype=float,
            )
            return {"annual_incidence": beta_series * multipliers}

        with patch(
            "ui.dynamic_ui.compute_ltbi_from_inc_hist",
            side_effect=fake_ltbi_from_inc_hist,
        ), patch(
            "ui.dynamic_ui.run_dynamic_model",
            side_effect=fake_run_dynamic_model,
        ):
            _, _, _, _, metadata = calibrate_beta_two_epoch_with_secular_decline(
                age_counts=age_counts,
                ages=ages,
                inc_hist=inc_hist,
                calib_years=4,
                risk_inputs=risk_inputs,
                pre_det_months=12.0,
                delta_pre=0.0,
                beta_bounds=(10.0, 50.0),
                adj_bounds=(1.0, 1.0),
                adj_grid_points=1,
                recent_years=4,
                projection_start_year=2025,
                secular_decline_start_year=2021,
                secular_decline_rate_bounds=(0.0, 0.20),
            )

        self.assertGreater(metadata["secular_decline_rate_annual"], 0.01)
        self.assertTrue(
            metadata["beta_recent_lower_bound_hit"]
            or metadata["beta_recent_near_lower_bound"]
        )
        self.assertLess(metadata["secular_multiplier_series"][-1], 1.0)
        json.dumps(metadata, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
