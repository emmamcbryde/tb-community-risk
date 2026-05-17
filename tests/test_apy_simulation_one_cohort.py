from __future__ import annotations

import math
import unittest

import numpy as np

from engine.apy.simulation import (
    COHORT_COLUMNS,
    RAW_FIELDS,
    safe_divide,
    safe_fraction,
    simulate_one_cohort,
    simulate_one_cohort_from_config,
)
from engine.apy.calibration import calibrate_from_config
from engine.apy.config import normalise_config
from engine.apy.cohort import make_rng
from engine.apy.regimen import (
    apply_regimen_overrides,
    default_regimen_library,
    get_regimen_from_library,
)


class ApySingleCohortSimulationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = normalise_config({"N": 100, "seed": 1})
        cls.calibration = calibrate_from_config(cls.config)
        library = default_regimen_library(
            cls.config.get("partialShortCourseMode") or "threshold80"
        )
        reg = get_regimen_from_library(cls.config["regimen"], library)
        cls.regimen = apply_regimen_overrides(reg, cls.config)

    def run_cached(self, overrides=None, n=100, seed=1, return_cohort=True):
        cfg = dict(self.config)
        if overrides:
            cfg.update(overrides)
        return simulate_one_cohort(
            n,
            self.calibration["parameters"],
            self.regimen,
            self.calibration,
            cfg,
            make_rng(seed),
            return_cohort=return_cohort,
        )

    def test_default_config_runs_small_single_cohort(self) -> None:
        result = simulate_one_cohort_from_config(n=100, seed=1)

        self.assertIsInstance(result["raw"], dict)
        self.assertEqual(len(result["cohort"]), 100)

    def test_raw_output_contains_required_fields(self) -> None:
        raw = self.run_cached(return_cohort=False)["raw"]

        for field in RAW_FIELDS:
            self.assertIn(field, raw)

    def test_cohort_dataframe_contains_required_columns(self) -> None:
        cohort = self.run_cached()["cohort"]

        self.assertEqual(list(cohort.columns), COHORT_COLUMNS)

    def test_counts_are_internally_consistent(self) -> None:
        result = self.run_cached(n=150, seed=2)
        raw = result["raw"]
        cohort = result["cohort"]

        self.assertEqual(raw["nScreened"], int(cohort["screened"].sum()))
        self.assertEqual(raw["nInfected"], int(cohort["infected"].sum()))
        self.assertEqual(raw["nStartTPT"], int(cohort["startedTPT"].sum()))
        self.assertEqual(raw["nCompleteTPT"], int(cohort["completedTPT"].sum()))
        self.assertEqual(
            raw["nPreventedActiveTB"], int(cohort["preventedActiveTB"].sum())
        )

    def test_summary_count_bounds(self) -> None:
        raw = self.run_cached(n=150, seed=3)["raw"]

        self.assertLessEqual(raw["nTestPositiveNonActive"], raw["nTestPositive"])
        self.assertLessEqual(raw["nFalsePositiveTreated"], raw["nTotalCoursesStarted"])
        self.assertLessEqual(raw["nTotalCoursesCompleted"], raw["nTotalCoursesStarted"])
        self.assertEqual(raw["nPartialCourses"], raw["nADRstop"] + raw["nStoppedOther"])
        self.assertLessEqual(raw["nActiveBy2y"], raw["nActiveBy20y"])

    def test_active_and_latent_at_screen_are_mutually_exclusive(self) -> None:
        cohort = self.run_cached(n=150, seed=4)["cohort"]

        self.assertFalse((cohort["activeAtScreen"] & cohort["latentAtScreen"]).any())

    def test_fixed_seed_is_reproducible(self) -> None:
        first = self.run_cached(n=120, seed=5)["raw"]
        second = self.run_cached(n=120, seed=5)["raw"]

        for field in RAW_FIELDS:
            left = first[field]
            right = second[field]
            if isinstance(left, float) and math.isnan(left):
                self.assertTrue(math.isnan(right), field)
            else:
                self.assertEqual(left, right, field)

    def test_different_seed_changes_stochastic_output(self) -> None:
        first = self.run_cached(n=200, seed=6)["raw"]
        second = self.run_cached(n=200, seed=7)["raw"]

        stochastic_fields = [
            "nInfected",
            "nTestPositive",
            "nStartTPT",
            "nActiveBy20y",
        ]
        self.assertTrue(any(first[field] != second[field] for field in stochastic_fields))

    def test_igra_and_tst_both_run(self) -> None:
        igra = self.run_cached({"testType": "IGRA"}, n=100, seed=8)
        tst = self.run_cached({"testType": "TST"}, n=100, seed=8)

        self.assertEqual(igra["raw"]["testType"], "IGRA")
        self.assertEqual(tst["raw"]["testType"], "TST")

    def test_all_screening_strategies_run(self) -> None:
        for strategy in ["random", "ltbi", "cure", "prevent"]:
            with self.subTest(strategy=strategy):
                result = self.run_cached({"screeningStrategy": strategy}, n=100, seed=9)
                self.assertEqual(result["raw"]["screeningStrategy"], strategy)

    def test_zero_screen_coverage_screens_no_one(self) -> None:
        result = self.run_cached({"screenCoverage": 0}, n=100, seed=10)

        self.assertEqual(result["raw"]["nScreened"], 0)
        self.assertEqual(result["raw"]["nStartTPT"], 0)
        self.assertFalse(result["cohort"]["screened"].any())

    def test_full_screen_coverage_screens_everyone(self) -> None:
        result = self.run_cached({"screenCoverage": 1}, n=100, seed=11)

        self.assertEqual(result["raw"]["nScreened"], 100)
        self.assertTrue(result["cohort"]["screened"].all())

    def test_bcg_attributable_false_positive_counts_are_bounded(self) -> None:
        raw = self.run_cached({"testType": "TST"}, n=200, seed=12)["raw"]

        self.assertLessEqual(
            raw["nExcessFalsePositiveTestsDueToBCG"], raw["nFalsePositiveTests"]
        )
        self.assertLessEqual(
            raw["nExcessCoursesStartedDueToBCG"], raw["nFalsePositiveTreated"]
        )

    def test_nns_and_nnt_fields_do_not_crash_zero_denominators(self) -> None:
        raw = self.run_cached({"screenCoverage": 0}, n=50, seed=13)["raw"]

        for field in [
            "NNS_cureInfection",
            "NNS_preventActiveTB",
            "NNS_falsePositiveTreated",
            "NNT_started_cureInfection",
            "NNT_started_preventActiveTB",
        ]:
            self.assertTrue(
                math.isfinite(raw[field]) or math.isinf(raw[field]) or math.isnan(raw[field])
            )

    def test_safe_fraction_and_safe_divide_match_matlab_conventions(self) -> None:
        self.assertTrue(math.isnan(safe_fraction(1, 0)))
        self.assertTrue(math.isinf(safe_divide(1, 0)))
        self.assertEqual(safe_fraction(1, 2), 0.5)
        self.assertEqual(safe_divide(4, 2), 2)

    def test_return_cohort_false_skips_dataframe(self) -> None:
        result = self.run_cached(n=50, seed=14, return_cohort=False)

        self.assertIsNone(result["cohort"])


if __name__ == "__main__":
    unittest.main()
