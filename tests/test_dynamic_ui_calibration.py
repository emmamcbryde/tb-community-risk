from __future__ import annotations

import sys
import types
import unittest
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

from ui.dynamic_ui import _scalar_metadata, calibrate_beta_two_epoch


class CalibrateBetaTwoEpochTests(unittest.TestCase):
    def test_scalar_metadata_excludes_beta_series_and_keeps_forward_beta(self) -> None:
        metadata = {
            "beta_historical": np.float64(2.5),
            "beta_recent": np.float64(1.25),
            "beta_forward": np.float64(1.25),
            "beta_series": np.array([2.5, 2.5, 1.25, 1.25], dtype=float),
            "recent_years": np.int64(2),
        }

        scalar = _scalar_metadata(metadata)

        self.assertNotIn("beta_series", scalar)
        self.assertEqual(scalar["beta_forward"], scalar["beta_recent"])
        self.assertEqual(scalar["beta_forward"], 1.25)
        self.assertIsInstance(scalar["beta_forward"], float)
        self.assertEqual(scalar["recent_years"], 2)
        self.assertIsInstance(scalar["recent_years"], int)

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


if __name__ == "__main__":
    unittest.main()
