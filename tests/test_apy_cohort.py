from __future__ import annotations

import unittest

import numpy as np

from engine.apy.calibration import calibrate_from_config
from engine.apy.ltbi_state import enable_development_compatibility_mode
from engine.apy.cohort import (
    add_targeting_scores,
    age_lookup,
    average_survival_to_random_screen,
    disease_multiplier,
    discrete_draw_values,
    draw_base_population,
    draw_untreated_active_times,
    exprnd_local,
    get_counterfactual_no_bcg_specificity,
    get_test_performance,
    infection_probability,
    make_rng,
    preventable_active_risk,
    select_screened_people,
)
from engine.apy.reference_loader import load_reference_scenario_config


class ApyCohortPrimitiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cfg = load_reference_scenario_config(
            "validation/matlab_reference/default_random_igra_3hp_N1500_seed1/scenario_config.json"
        )
        cfg = enable_development_compatibility_mode(cfg)
        cls.calibration = calibrate_from_config(cfg)
        cls.pars = cls.calibration["parameters"]

    def test_rng_and_discrete_draws_are_reproducible(self) -> None:
        rng1 = make_rng(123)
        rng2 = make_rng(123)

        draws1 = discrete_draw_values([10, 20], [0.25, 0.75], 20, rng1)
        draws2 = discrete_draw_values([10, 20], [0.25, 0.75], 20, rng2)

        np.testing.assert_array_equal(draws1, draws2)
        self.assertTrue(set(draws1).issubset({10, 20}))

    def test_age_lookup_uses_one_based_age_groups(self) -> None:
        values = age_lookup([0.1, 0.2, 0.3], [1, 2, 3, 3])

        np.testing.assert_allclose(values, [0.1, 0.2, 0.3, 0.3])

    def test_exprnd_local_uses_mean_scale(self) -> None:
        draws = exprnd_local([1.0, 2.0, 3.0], make_rng(1))

        self.assertEqual(draws.shape, (3,))
        self.assertTrue(np.all(draws >= 0))

    def test_infection_probability_and_disease_multiplier_shapes(self) -> None:
        age_years = np.array([5, 25, 60])
        flags = np.array([False, True, True])

        p_inf = infection_probability(
            age_years,
            self.pars,
            self.calibration["ageInfLogLambda"],
            self.calibration["ageInfGamma"],
            flags,
            flags,
            flags,
        )
        mult = disease_multiplier(
            self.pars, flags, flags, flags, flags, flags, flags, flags
        )

        self.assertEqual(p_inf.shape, (3,))
        self.assertTrue(np.all((p_inf >= 0) & (p_inf <= 1)))
        self.assertEqual(mult.shape, (3,))
        self.assertEqual(mult[0], 1)
        self.assertGreater(mult[1], 1)

    def test_survival_and_preventable_risk_helpers(self) -> None:
        mult = np.array([1.0, 2.0, 3.0])
        avg_surv = average_survival_to_random_screen(
            mult, self.calibration["lambdaEarly"], 2
        )
        preventable = preventable_active_risk(
            mult,
            self.calibration["lambdaEarly"],
            self.calibration["lambdaLate"],
            2,
            20,
        )

        self.assertTrue(np.all((avg_surv >= 0) & (avg_surv <= 1)))
        self.assertTrue(np.all(preventable >= 0))

    def test_select_screened_people_counts_and_ranks(self) -> None:
        rng = make_rng(42)
        ltbi = np.array([0.1, 0.9, 0.3, 0.2])
        screened, score, rank = select_screened_people(
            ltbi, ltbi, ltbi, 0.5, "ltbi", rng
        )

        self.assertEqual(int(screened.sum()), 2)
        self.assertTrue(screened[1])
        self.assertEqual(score.tolist(), ltbi.tolist())
        self.assertEqual(sorted(rank.tolist()), [1, 2, 3, 4])

    def test_test_performance_helpers(self) -> None:
        bcg = np.array([True, False])
        igra_sens, igra_spec = get_test_performance(
            bcg, {"testType": "IGRA", "testSensitivity": 0.95, "testSpecificity": 0.98}
        )
        tst_sens, tst_spec = get_test_performance(
            bcg,
            {
                "testType": "TST",
                "tstSensitivity": 0.8,
                "tstSpecificityBCG": 0.55,
                "tstSpecificityNoBCG": 0.97,
            },
        )
        no_bcg_spec = get_counterfactual_no_bcg_specificity(
            bcg, {"testType": "TST", "tstSpecificityNoBCG": 0.97}
        )

        np.testing.assert_allclose(igra_sens, [0.95, 0.95])
        np.testing.assert_allclose(igra_spec, [0.98, 0.98])
        np.testing.assert_allclose(tst_sens, [0.8, 0.8])
        np.testing.assert_allclose(tst_spec, [0.55, 0.97])
        np.testing.assert_allclose(no_bcg_spec, [0.97, 0.97])

    def test_draw_base_population_outputs_expected_fields_and_shapes(self) -> None:
        population = draw_base_population(
            50,
            self.pars,
            self.calibration["ageInfLogLambda"],
            self.calibration["ageInfGamma"],
            make_rng(99),
        )

        expected = {
            "ageYears",
            "ageGroup",
            "female",
            "BCG",
            "MJ",
            "contact",
            "renal",
            "diabetes",
            "smoking",
            "cld",
            "alcohol",
            "pInfection",
            "infected",
            "diseaseMultiplier",
        }
        self.assertEqual(set(population), expected)
        for value in population.values():
            self.assertEqual(len(value), 50)
        self.assertTrue(np.all((population["pInfection"] >= 0) & (population["pInfection"] <= 1)))

    def test_targeting_scores_and_active_times(self) -> None:
        population = draw_base_population(
            30,
            self.pars,
            self.calibration["ageInfLogLambda"],
            self.calibration["ageInfGamma"],
            make_rng(100),
        )
        scored = add_targeting_scores(
            population,
            self.calibration["lambdaEarly"],
            self.calibration["lambdaLate"],
            2,
            20,
        )
        t_active = draw_untreated_active_times(
            scored["infected"],
            scored["diseaseMultiplier"],
            self.calibration["lambdaEarly"],
            self.calibration["lambdaLate"],
            2,
            make_rng(101),
        )

        self.assertIn("ltbiRiskScore", scored)
        self.assertIn("cureTargetScore", scored)
        self.assertIn("preventTargetScore", scored)
        self.assertEqual(t_active.shape, (30,))
        self.assertTrue(np.all(np.isinf(t_active[~scored["infected"]])))


if __name__ == "__main__":
    unittest.main()
