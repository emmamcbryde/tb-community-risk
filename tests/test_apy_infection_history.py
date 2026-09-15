from __future__ import annotations

import unittest
from pathlib import Path

from engine.apy.config import build_default_config, normalise_config
from engine.apy.data import load_parameters_from_config
from engine.apy.expected_value import run_expected_value
from engine.apy.infection_history import (
    EXPERIMENTAL_ECONOMICS_GUARD_MESSAGE,
    EXPERIMENTAL_READINESS_ITEMS,
    EXPLICIT_EXPERIMENTAL_BASIS,
    TRAJECTORY_PRESETS,
    calibrate_infection_history,
    configure_compatibility_reference_assumptions,
    configure_infection_history_assumptions,
    has_explicit_infection_history,
    infection_history_readiness_status,
    is_experimental_infection_history_results,
)
from engine.apy.results_bundle import build_results_bundle
from engine.apy.runner import run_replicates


class ApyInfectionHistoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = normalise_config(build_default_config())
        cls.parameters = load_parameters_from_config(cls.config)

    def test_trajectory_presets_reproduce_prevalence_and_age_or(self) -> None:
        target_prev = 47 / 624
        target_or = self.config["targetAgeOR"]
        fractions = {}
        for trajectory in ("rising", "steady", "falling"):
            derived = calibrate_infection_history(
                self.parameters,
                target_prevalence=target_prev,
                target_age_or=target_or,
                trajectory=trajectory,
            )
            self.assertAlmostEqual(derived["expectedPrevalence"], target_prev, places=8)
            self.assertAlmostEqual(derived["expectedAgeOR"], target_or, places=4)
            self.assertGreaterEqual(derived["recentFraction"], 0.0)
            self.assertLessEqual(derived["recentFraction"], 1.0)
            self.assertAlmostEqual(
                derived["recentFraction"] + derived["remoteFraction"],
                1.0,
                places=10,
            )
            for row in derived["ageRows"]:
                self.assertGreaterEqual(row["ltbiPrevalence"], 0.0)
                self.assertLessEqual(row["ltbiPrevalence"], 1.0)
                self.assertGreaterEqual(row["recentFractionAmongPrevalent"], 0.0)
                self.assertLessEqual(row["recentFractionAmongPrevalent"], 1.0)
            fractions[trajectory] = derived["recentFraction"]
        self.assertGreater(fractions["rising"], fractions["steady"])
        self.assertGreater(fractions["steady"], fractions["falling"])

    def test_trajectory_definitions_are_scenario_assumptions(self) -> None:
        self.assertEqual(TRAJECTORY_PRESETS["rising"]["trendRatePerYear"], 0.01)
        self.assertEqual(TRAJECTORY_PRESETS["steady"]["trendRatePerYear"], 0.0)
        self.assertEqual(TRAJECTORY_PRESETS["falling"]["trendRatePerYear"], -0.01)

    def test_config_records_explicit_scientific_scenario_without_changing_reference(self) -> None:
        explicit = configure_infection_history_assumptions(self.config, "rising")
        self.assertTrue(has_explicit_infection_history(explicit))
        self.assertEqual(explicit["analysisBasis"], EXPLICIT_EXPERIMENTAL_BASIS)
        self.assertEqual(
            explicit["ltbiStateAssumptions"]["baselineRecentLTBIDerivationMethod"],
            "infection_history_trajectory",
        )
        self.assertIn("preceding two years", explicit["ltbiStateAssumptions"]["stateDefinition"])
        self.assertIn("five-year mean residence", explicit["ltbiStateAssumptions"]["stateDefinition"])
        restored = configure_compatibility_reference_assumptions(explicit)
        self.assertFalse(has_explicit_infection_history(restored))
        self.assertEqual(restored["analysisBasis"], "sa_health_matlab_v9_compatibility_reference")
        self.assertIsNone(restored["ltbiStateAssumptions"]["baselineRecentLTBIProportion"])

    def test_fresh_default_config_is_not_experimental(self) -> None:
        self.assertFalse(has_explicit_infection_history(self.config))
        self.assertNotEqual(self.config.get("analysisBasis"), EXPLICIT_EXPERIMENTAL_BASIS)

    def test_expected_value_uses_selected_explicit_infection_history(self) -> None:
        config = configure_infection_history_assumptions(self.config, "steady")
        config["analysisMethod"] = "expected_value"
        config["N"] = 20
        result = run_expected_value(config)
        bundle = build_results_bundle(
            {
                **result,
                "raw": __import__("pandas").DataFrame([{"nScreened": 0}]),
                "summary": __import__("pandas").DataFrame(
                    [{"Metric": "nScreened", "Median": 0, "Low95": 0, "High95": 0}]
                ),
            }
        )
        self.assertEqual(bundle["metadata"]["analysisBasis"], EXPLICIT_EXPERIMENTAL_BASIS)
        self.assertTrue(is_experimental_infection_history_results(bundle))
        derived = result["calibration"]["infectionHistory"]
        totals = result["eventLedger"]["replicateTotals"]
        intervention = totals[totals["arm"] == "intervention"]
        recent = float(intervention[intervention["eventName"] == "recent_ltbi_at_baseline"]["value"].iloc[0])
        remote = float(intervention[intervention["eventName"] == "remote_ltbi_at_baseline"]["value"].iloc[0])
        infected = float(intervention[intervention["eventName"] == "infected_at_baseline"]["value"].iloc[0])
        self.assertGreater(recent, 0.0)
        self.assertGreater(remote, 0.0)
        self.assertAlmostEqual(recent + remote, infected, places=8)
        self.assertAlmostEqual(recent / infected, derived["recentFraction"], places=3)

    def test_experimental_readiness_remains_not_validated(self) -> None:
        readiness = infection_history_readiness_status()
        self.assertEqual(readiness["status"], "not_validated_for_decision_use")
        self.assertEqual(len(readiness["items"]), len(EXPERIMENTAL_READINESS_ITEMS))
        unresolved = {item["status"] for item in readiness["items"]}
        self.assertEqual(unresolved, {"unresolved"})

    def test_standard_pages_guard_experimental_economics_and_activation(self) -> None:
        start_page = Path("pages/0_Start.py").read_text(encoding="utf-8")
        economics_page = Path("pages/4_Economics.py").read_text(encoding="utf-8")
        results_page = Path("pages/3_Results.py").read_text(encoding="utf-8")
        self.assertIn("Experimental infection-history analysis", start_page)
        self.assertIn("Enable experimental infection-history analysis", start_page)
        self.assertIn("Restore SA Health report reference", start_page)
        self.assertNotIn("Epidemiological basis\",\\n        basis_options", start_page)
        self.assertIn("is_experimental_infection_history_results", economics_page)
        self.assertIn("st.stop()", economics_page)
        self.assertIn("is_experimental_infection_history_results", results_page)
        self.assertIn("EXPERIMENTAL_ECONOMICS_GUARD_MESSAGE", economics_page)
        self.assertTrue(EXPERIMENTAL_ECONOMICS_GUARD_MESSAGE.startswith("Health-economic conclusions"))

    def test_stochastic_runner_records_selected_explicit_infection_history(self) -> None:
        config = configure_infection_history_assumptions(self.config, "falling")
        config["N"] = 50
        config["nReps"] = 1
        config["seed"] = 123
        result = run_replicates(config, keep_example_cohort=False)
        self.assertEqual(result["calibration"]["infectionHistory"]["trajectory"], "falling")
        self.assertEqual(
            result["strategy"]["baselineRecentLTBIDerivationMethod"],
            "infection_history_trajectory",
        )
        raw = result["raw"].iloc[0]
        self.assertEqual(
            raw["nRecentLTBIAtBaseline"] + raw["nRemoteLTBIAtBaseline"],
            raw["nInfected"],
        )


if __name__ == "__main__":
    unittest.main()
