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
    _scalar_metadata,
    calibrate_beta_two_epoch,
    two_epoch_diagnostics_display_tables,
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
        json.dumps(scalar, allow_nan=False)

    def test_two_epoch_display_tables_include_requested_diagnostics(self) -> None:
        metadata = {
            "calibration_mode": "Two-epoch beta: historical + recent 10 years",
            "beta_historical": 2.0,
            "beta_recent": 1.0,
            "beta_ratio_recent_to_historical": 0.5,
            "beta_recent_start_year": 2023,
            "recent_years": 2,
            "rmse_overall": 3.0,
            "rmse_recent": 1.5,
            "beta_historical_years": [2020, 2021, 2022],
            "beta_recent_years": [2023, 2024],
        }

        summary_df, preview_df = two_epoch_diagnostics_display_tables(metadata)

        summary = dict(zip(summary_df["Metric"], summary_df["Value"]))
        self.assertEqual(summary["beta_historical"], 2.0)
        self.assertEqual(summary["beta_recent"], 1.0)
        self.assertEqual(summary["beta_recent/beta_historical"], 0.5)
        self.assertEqual(summary["beta_recent_start_year"], 2023)
        self.assertEqual(summary["recent_years"], 2)
        self.assertEqual(summary["rmse_overall"], 3.0)
        self.assertEqual(summary["rmse_recent"], 1.5)
        self.assertEqual(preview_df["Epoch"].tolist(), ["historical", "recent"])
        self.assertEqual(preview_df["Years"].tolist(), ["2020-2022", "2023-2024"])

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
            "beta_historical",
            "beta_recent",
            "beta_ratio_recent_to_historical",
            "beta_ratio",
            "recent_years",
            "projection_start_year",
            "beta_recent_start_year",
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
        self.assertEqual(metadata["beta_forward"], metadata["beta_recent"])
        self.assertEqual(len(metadata["fitted_incidence"]), 4)
        self.assertEqual(len(metadata["target_incidence"]), 4)
        self.assertEqual(len(metadata["residuals"]), 4)
        json.dumps(metadata)


if __name__ == "__main__":
    unittest.main()
