from __future__ import annotations

import json
import unittest
from pathlib import Path

from engine.apy.calibration import (
    age_cumulative_infection_hazard,
    bern_prob,
    calibrate_early_hazard,
    calibrate_from_config,
    disease_multiplier_from_flags,
    expected_active_within_window,
    expected_age_or_exact,
    expected_infection_prevalence_exact,
    infection_probability_from_eta,
    odds_from_prob,
)
from engine.apy.reference_loader import load_reference_scenario_config


FIXTURE_DIR = (
    Path("validation") / "matlab_reference" / "default_random_igra_3hp_N1500_seed1"
)


class ApyCalibrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        snapshot_path = FIXTURE_DIR / "matlab_parameter_snapshot.json"
        if not snapshot_path.is_file():
            raise unittest.SkipTest("MATLAB parameter snapshot fixture is unavailable.")
        cls.snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
        cls.pars = cls.snapshot["parameters"]
        cls.calibration = cls.snapshot["calibration"]
        cls.scenario_config = load_reference_scenario_config(
            FIXTURE_DIR / "scenario_config.json"
        )

    def test_basic_probability_helpers(self) -> None:
        self.assertAlmostEqual(bern_prob(True, 0.2), 0.2)
        self.assertAlmostEqual(bern_prob(False, 0.2), 0.8)
        self.assertAlmostEqual(odds_from_prob(0.25), 1 / 3)
        self.assertAlmostEqual(
            float(age_cumulative_infection_hazard(0, -8, 1)),
            0.5 * 0.00033546262790251185,
        )
        self.assertAlmostEqual(
            float(infection_probability_from_eta(0)),
            1 - 2.718281828459045 ** -1,
        )

    def test_expected_infection_matches_matlab_snapshot(self) -> None:
        expected = expected_infection_prevalence_exact(
            self.calibration["ageInfLogLambda"],
            self.calibration["ageInfGamma"],
            self.pars,
        )
        age_or = expected_age_or_exact(
            self.calibration["ageInfLogLambda"],
            self.calibration["ageInfGamma"],
            self.pars,
        )

        self.assertAlmostEqual(expected, self.calibration["expectedInfPrev"], places=10)
        self.assertAlmostEqual(age_or, self.calibration["expectedAgeOR"], places=9)

    def test_disease_multiplier_uses_flags(self) -> None:
        mult = disease_multiplier_from_flags(
            self.pars, mj=True, contact=True, renal=True, diabetes=True
        )

        self.assertAlmostEqual(mult, 3 * 5 * 3.6 * 3)

    def test_expected_active_within_window_matches_matlab_snapshot(self) -> None:
        expected = expected_active_within_window(
            self.calibration["lambdaEarly"],
            self.pars,
            self.calibration["ageInfLogLambda"],
            self.calibration["ageInfGamma"],
            2,
        )

        self.assertAlmostEqual(expected, self.calibration["expectedActive2y"], places=10)

    def test_legacy_single_state_calibrate_early_hazard_matches_matlab_snapshot(self) -> None:
        hazard = calibrate_early_hazard(
            self.pars,
            self.calibration["ageInfLogLambda"],
            self.calibration["ageInfGamma"],
            self.calibration["targetActive2y"],
            2,
            5,
            baseline_recent_proportion=1.0,
            recent_to_remote_rate=0.0,
        )

        self.assertAlmostEqual(hazard["lambdaEarly"], self.calibration["lambdaEarly"], places=10)
        self.assertAlmostEqual(hazard["lambdaLate"], self.calibration["lambdaLate"], places=10)

    def test_calibrate_early_hazard_rejects_invalid_ratio(self) -> None:
        with self.assertRaisesRegex(ValueError, "earlyLateRatio"):
            calibrate_early_hazard(
                self.pars,
                self.calibration["ageInfLogLambda"],
                self.calibration["ageInfGamma"],
                self.calibration["targetActive2y"],
                2,
                0.5,
            )

    def test_calibrate_from_config_keeps_age_calibration_and_records_ltbi_state_assumptions(self) -> None:
        calibrated = calibrate_from_config(self.scenario_config)

        self.assertAlmostEqual(
            calibrated["ageInfLogLambda"],
            self.calibration["ageInfLogLambda"],
            places=9,
        )
        self.assertAlmostEqual(
            calibrated["ageInfGamma"], self.calibration["ageInfGamma"], places=9
        )
        self.assertNotAlmostEqual(calibrated["lambdaEarly"], self.calibration["lambdaEarly"], places=10)
        self.assertAlmostEqual(
            calibrated["expectedInfPrev"],
            self.calibration["expectedInfPrev"],
            places=10,
        )
        self.assertAlmostEqual(
            calibrated["expectedAgeOR"], self.calibration["expectedAgeOR"], places=8
        )
        self.assertEqual(
            calibrated["ltbiStateAssumptionStatus"],
            "unresolved_compatibility_placeholder",
        )
        self.assertEqual(calibrated["baselineRecentLTBIProportion"], 0.0)


if __name__ == "__main__":
    unittest.main()
