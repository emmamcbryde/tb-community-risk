from __future__ import annotations

import unittest

from engine.apy.health_economics import (
    build_default_health_economic_assumptions,
    calculate_health_outcomes,
    calculate_icers,
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

    def test_zero_daly_qaly_denominators_are_safe(self) -> None:
        health = {"healthOutcomes": {"dalysAverted": 0, "qalysGained": 0, "tbCasesPrevented": 0}}
        economics = {"costs": {"netCostVsBaseline": 1000}}

        out = calculate_icers(economics, health)

        self.assertIsNone(out["icers"]["costPerDALYAverted"])
        self.assertIsNone(out["icers"]["costPerQALYGained"])
        self.assertIsNone(out["icers"]["costPerTBCasePrevented"])


if __name__ == "__main__":
    unittest.main()
