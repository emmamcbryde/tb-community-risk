from __future__ import annotations

import unittest

from engine.apy.health_economics import (
    build_default_health_economic_assumptions,
    calculate_health_outcomes,
    calculate_icers,
    calculate_post_tb_sequelae,
)


def mock_results(cases_prevented: float = 10.0) -> dict:
    return {
        "summary": [
            {"Metric": "nPreventedActiveTB", "Median": cases_prevented},
            {"Metric": "nTotalCoursesStarted", "Median": 100.0},
            {"Metric": "nADRstop", "Median": 5.0},
        ],
        "strategy": {
            "regimenMonths": 3,
        },
    }


class ApyHealthEconomicsTests(unittest.TestCase):
    def test_default_health_economic_assumptions_load(self) -> None:
        assumptions = build_default_health_economic_assumptions()

        self.assertEqual(assumptions["metadata"]["currencyCode"], "AUD")
        self.assertEqual(assumptions["metadata"]["priceYear"], 2019)
        self.assertEqual(assumptions["daly"]["activeTBDisabilityWeight"], 0.333)
        self.assertTrue(assumptions["qaly"]["includeMortalityInQaly"])
        self.assertEqual(assumptions["qaly"]["qalyMortalityMethod"], "yll_times_healthy_utility")
        self.assertTrue(assumptions["qaly"]["gbdAlignedQalyMorbiditySensitivity"])
        self.assertEqual(assumptions["postTB"]["postTBDalysPerCase"], 5.8)
        self.assertEqual(assumptions["postTB"]["postTBQalysLostPerCase"], 0.93)
        self.assertFalse(assumptions["postTB"]["includePostTBSequelae"])
        self.assertEqual(assumptions["thresholds"]["wtpLow"], 45000)

    def test_daly_calculation_works_for_known_tb_cases_prevented(self) -> None:
        assumptions = {
            "daly": {
                "tbCaseFatalityRisk": 0.1,
                "yllPerTBDeath": 20.0,
            }
        }

        out = calculate_health_outcomes(mock_results(10), assumptions)

        self.assertAlmostEqual(out["healthOutcomes"]["tbDeathsPrevented"], 1.0)
        self.assertAlmostEqual(out["healthOutcomes"]["yldAverted"], 1.665)
        self.assertAlmostEqual(out["healthOutcomes"]["yllAverted"], 20.0)
        self.assertAlmostEqual(out["healthOutcomes"]["dalysAverted"], 21.665)

    def test_yld_only_calculation_is_disability_weight_times_duration(self) -> None:
        out = calculate_health_outcomes(
            mock_results(1),
            {"daly": {"includeYLL": False}},
        )

        self.assertAlmostEqual(out["healthOutcomes"]["yldAverted"], 0.333 * 0.5)
        self.assertAlmostEqual(out["healthOutcomes"]["yllAverted"], 0.0)
        self.assertAlmostEqual(out["healthOutcomes"]["dalysAverted"], 0.333 * 0.5)

    def test_yll_calculation_uses_deaths_prevented_times_yll(self) -> None:
        out = calculate_health_outcomes(
            mock_results(4),
            {"daly": {"tbCaseFatalityRisk": 0.25, "yllPerTBDeath": 10}},
        )

        self.assertAlmostEqual(out["healthOutcomes"]["tbDeathsPrevented"], 1.0)
        self.assertAlmostEqual(out["healthOutcomes"]["yllAverted"], 10.0)

    def test_qaly_active_tb_loss_calculation_works(self) -> None:
        out = calculate_health_outcomes(mock_results(10))
        expected = 10 * (0.8733 - 0.8182) * 0.5

        self.assertAlmostEqual(out["healthOutcomes"]["qalyLossActiveTBAverted"], expected)
        self.assertAlmostEqual(
            out["healthOutcomes"]["activeTbMorbidityQalysGained_Dale"],
            expected,
        )

    def test_mortality_qalys_are_added_when_enabled(self) -> None:
        out = calculate_health_outcomes(
            mock_results(10),
            {"daly": {"tbCaseFatalityRisk": 0.1, "yllPerTBDeath": 20.0}},
        )
        health = out["healthOutcomes"]

        self.assertAlmostEqual(health["tbDeathsPrevented"], 1.0)
        self.assertAlmostEqual(
            health["qualityAdjustedLifeExpectancyPerTBDeath"],
            20.0 * 0.8733,
        )
        self.assertAlmostEqual(
            health["mortalityQalysGained"],
            health["tbDeathsPrevented"] * health["qualityAdjustedLifeExpectancyPerTBDeath"],
        )

    def test_dale_morbidity_qaly_per_case_uses_utility_decrement(self) -> None:
        out = calculate_health_outcomes(mock_results(1))

        self.assertAlmostEqual(
            out["healthOutcomes"]["activeTbMorbidityQalysGained_Dale"],
            (0.8733 - 0.8182) * 0.5,
        )

    def test_gbd_aligned_morbidity_qaly_per_case_uses_disability_weight(self) -> None:
        out = calculate_health_outcomes(mock_results(1))

        self.assertAlmostEqual(
            out["healthOutcomes"]["activeTbMorbidityQalysGained_GBD"],
            0.333 * 0.5,
        )

    def test_primary_qalys_exceed_morbidity_only_when_mortality_included(self) -> None:
        out = calculate_health_outcomes(mock_results(10))
        health = out["healthOutcomes"]

        self.assertGreater(
            health["qalysGained"],
            health["morbidityOnlyQalysGained"],
        )
        self.assertAlmostEqual(
            health["qalysGained"],
            health["qalysGained_DaleMortalityAdjusted"],
        )

    def test_dalys_and_mortality_adjusted_qalys_are_closer_than_morbidity_only_qalys(self) -> None:
        out = calculate_health_outcomes(mock_results(10))
        health = out["healthOutcomes"]

        dalys = health["dalysAverted"]
        mortality_adjusted_qalys = health["qalysGained_DaleMortalityAdjusted"]
        morbidity_only_qalys = health["morbidityOnlyQalysGained"]

        self.assertLess(
            abs(dalys - mortality_adjusted_qalys),
            abs(dalys - morbidity_only_qalys),
        )

    def test_ltbi_treatment_decrement_reduces_qalys_gained(self) -> None:
        base = calculate_health_outcomes(mock_results(10))
        decrement = calculate_health_outcomes(
            mock_results(10),
            {"qaly": {"ltbiTreatmentDecrement": 0.0133}},
        )

        self.assertLess(
            decrement["healthOutcomes"]["qalysGained"],
            base["healthOutcomes"]["qalysGained"],
        )

    def test_zero_deaths_prevented_is_handled_safely(self) -> None:
        out = calculate_health_outcomes(
            mock_results(10),
            {"daly": {"tbCaseFatalityRisk": 0.0}},
        )

        self.assertAlmostEqual(out["healthOutcomes"]["tbDeathsPrevented"], 0.0)
        self.assertAlmostEqual(out["healthOutcomes"]["mortalityQalysGained"], 0.0)

    def test_post_tb_dalys_per_case_are_added_correctly(self) -> None:
        out = calculate_post_tb_sequelae(
            tb_cases_prevented=4,
            acute_dalys_averted=6.586,
            acute_qalys_gained=5.28,
            acute_net_cost=-19956.1696,
            assumptions={"postTB": {"includePostTBSequelae": True}},
        )

        self.assertAlmostEqual(out["postTBScenarios"]["postTBDALYsAverted"], 4 * 5.8)
        self.assertAlmostEqual(
            out["postTBScenarios"]["totalDALYsAvertedIncludingPostTB"],
            6.586 + 4 * 5.8,
        )

    def test_post_tb_qalys_per_case_are_added_correctly(self) -> None:
        out = calculate_post_tb_sequelae(
            tb_cases_prevented=4,
            acute_dalys_averted=6.586,
            acute_qalys_gained=5.28,
            acute_net_cost=-19956.1696,
            assumptions={"postTB": {"includePostTBSequelae": True}},
        )

        self.assertAlmostEqual(out["postTBScenarios"]["postTBQALYsGained"], 4 * 0.93)
        self.assertAlmostEqual(
            out["postTBScenarios"]["totalQALYsGainedIncludingPostTB"],
            5.28 + 4 * 0.93,
        )

    def test_acute_only_results_are_unchanged_when_post_tb_disabled(self) -> None:
        out = calculate_post_tb_sequelae(
            tb_cases_prevented=4,
            acute_dalys_averted=6.586,
            acute_qalys_gained=5.28,
            acute_net_cost=-19956.1696,
        )

        self.assertEqual(out["postTBScenarios"]["postTBDALYsAverted"], 0.0)
        self.assertEqual(out["postTBScenarios"]["postTBQALYsGained"], 0.0)
        self.assertAlmostEqual(out["postTBScenarios"]["totalDALYsAvertedIncludingPostTB"], 6.586)
        self.assertAlmostEqual(out["postTBScenarios"]["totalQALYsGainedIncludingPostTB"], 5.28)
        self.assertAlmostEqual(out["postTBScenarios"]["netCostIncludingPostTB"], -19956.1696)

    def test_post_tb_icers_update_when_tail_is_included(self) -> None:
        acute_cost_per_daly = -167627.7536 / 19.758
        acute_cost_per_qaly = -167627.7536 / 15.8401318356
        out = calculate_post_tb_sequelae(
            tb_cases_prevented=12,
            acute_dalys_averted=19.758,
            acute_qalys_gained=15.8401318356,
            acute_net_cost=-167627.7536,
            assumptions={"postTB": {"includePostTBSequelae": True}},
        )

        self.assertNotEqual(
            out["postTBScenarios"]["costPerDALYIncludingPostTB"],
            acute_cost_per_daly,
        )
        self.assertNotEqual(
            out["postTBScenarios"]["costPerQALYIncludingPostTB"],
            acute_cost_per_qaly,
        )

    def test_zero_post_tb_cost_is_handled(self) -> None:
        out = calculate_post_tb_sequelae(
            tb_cases_prevented=11,
            acute_dalys_averted=18.1115,
            acute_qalys_gained=14.5200978356,
            acute_net_cost=-147720.6176,
            assumptions={"postTB": {"includePostTBSequelae": True, "postTBAnnualCareCost": 0}},
        )

        self.assertEqual(out["postTBScenarios"]["postTBCostsAverted"], 0.0)
        self.assertAlmostEqual(out["postTBScenarios"]["netCostIncludingPostTB"], -147720.6176)

    def test_editable_post_tb_annual_care_cost_is_included(self) -> None:
        out = calculate_post_tb_sequelae(
            tb_cases_prevented=2,
            acute_dalys_averted=3.0,
            acute_qalys_gained=2.0,
            acute_net_cost=1000.0,
            assumptions={
                "postTB": {
                    "includePostTBSequelae": True,
                    "postTBAnnualCareCost": 250.0,
                    "postTBDurationYears": 20,
                }
            },
        )

        self.assertAlmostEqual(out["postTBScenarios"]["postTBCostsAverted"], 2 * 250.0 * 20)
        self.assertAlmostEqual(out["postTBScenarios"]["netCostIncludingPostTB"], -9000.0)

    def test_icer_per_daly_uses_incremental_cost_over_dalys(self) -> None:
        health = {"healthOutcomes": {"dalysAverted": 10, "qalysGained": 5, "tbCasesPrevented": 2}}
        economics = {"costs": {"netCostVsBaseline": 1000}}

        out = calculate_icers(economics, health)

        self.assertEqual(out["icers"]["costPerDALYAverted"], 100)

    def test_icer_per_qaly_uses_incremental_cost_over_qalys(self) -> None:
        health = {"healthOutcomes": {"dalysAverted": 10, "qalysGained": 5, "tbCasesPrevented": 2}}
        economics = {"costs": {"netCostVsBaseline": 1000}}

        out = calculate_icers(economics, health)

        self.assertEqual(out["icers"]["costPerQALYGained"], 200)

    def test_icer_includes_gbd_aligned_qaly_sensitivity(self) -> None:
        health = {
            "healthOutcomes": {
                "dalysAverted": 10,
                "qalysGained": 5,
                "qalysGained_DaleMortalityAdjusted": 5,
                "qalysGained_GBDAlignedMortalityAdjusted": 8,
                "tbCasesPrevented": 2,
            }
        }
        economics = {"costs": {"netCostVsBaseline": 1000}}

        out = calculate_icers(economics, health)

        self.assertEqual(out["icers"]["costPerQALYGained"], 200)
        self.assertEqual(out["icers"]["costPerQALYGained_GBDAligned"], 125)
        self.assertEqual(out["nmb"]["netMonetaryBenefitQALY_GBDAligned_low"], 8 * 45000 - 1000)

    def test_zero_daly_qaly_denominators_are_safe(self) -> None:
        health = {"healthOutcomes": {"dalysAverted": 0, "qalysGained": 0, "tbCasesPrevented": 0}}
        economics = {"costs": {"netCostVsBaseline": 1000}}

        out = calculate_icers(economics, health)

        self.assertIsNone(out["icers"]["costPerDALYAverted"])
        self.assertIsNone(out["icers"]["costPerQALYGained"])
        self.assertIsNone(out["icers"]["costPerTBCasePrevented"])


if __name__ == "__main__":
    unittest.main()
