from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from streamlit.testing.v1 import AppTest

from app.state import is_supported_sa_health_analysis_basis
from app.run_analysis_controls import TECHNICAL_DEMONSTRATION_ROUTE, prepare_run_config_for_recent_ltbi_route
from engine.apy.config import build_default_config, normalise_config
from engine.apy.data import load_parameters_from_config
from engine.apy.expected_value import run_expected_value
from engine.apy.infection_history import (
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
from engine.apy.working_defaults import build_unified_working_default_preset


class ApyInfectionHistoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = normalise_config(build_default_config())
        cls.parameters = load_parameters_from_config(cls.config)
        cls.economics_config = build_unified_working_default_preset()["economicsConfig"]

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

    def test_standard_pages_do_not_expose_retired_infection_history_pathway(self) -> None:
        start_page = Path("pages/0_Start.py").read_text(encoding="utf-8")
        economics_source = Path("pages/4_Economics.py").read_text(encoding="utf-8")
        economics_page = economics_source.split("st.stop()", 1)[0]
        results_page = Path("pages/3_Results.py").read_text(encoding="utf-8")
        decision_page = Path("pages/5_Decision_Analysis.py").read_text(encoding="utf-8")
        evidence_page = Path("pages/6_Evidence_Assumptions.py").read_text(encoding="utf-8")
        combined = "\n".join([start_page, results_page, decision_page, evidence_page])
        for forbidden in [
            "Enable experimental infection-history analysis",
            "Historical TB infection pressure",
            "Rising",
            "Steady",
            "Falling",
            "recent fraction",
        ]:
            self.assertNotIn(forbidden, combined)
        self.assertNotIn("EXPERIMENTAL_ECONOMICS_GUARD_MESSAGE", economics_page)
        self.assertIn("sanitize_reference_only_state", start_page)
        self.assertIn("sanitize_reference_only_state", results_page)
        self.assertIn("sanitize_reference_only_state", economics_page)

    def test_rendered_reference_only_setup_has_no_trajectory_controls(self) -> None:
        app = AppTest.from_file(str(Path("pages/0_Start.py")))
        app.run(timeout=90)

        self.assertFalse(app.exception)
        self.assertNotIn("Historical TB infection pressure", [item.label for item in app.selectbox])
        self.assertNotIn("Enable experimental infection-history analysis", [item.label for item in app.button])
        self.assertEqual([button.label for button in app.button].count("Restore APY defaults"), 1)

    def test_rendered_stale_trajectory_config_is_rejected_and_defaults_restored(self) -> None:
        stale_config = configure_infection_history_assumptions(self.config, "falling")
        app = AppTest.from_file(str(Path("pages/0_Start.py")))
        app.session_state["config"] = stale_config
        app.session_state["economics_config"] = self.economics_config
        app.session_state["experimental_infection_history_enabled"] = False
        app.session_state["infection_history_trajectory_label"] = "Falling"

        app.run(timeout=90)

        self.assertFalse(app.exception)
        self.assertFalse(has_explicit_infection_history(app.session_state["config"]))
        self.assertNotIn("experimental_infection_history_enabled", app.session_state)
        self.assertNotIn("Historical TB infection pressure", [item.label for item in app.selectbox])
        self.assertIn(
            "Future TB outcomes use the assumptions applied in the SA Health report.",
            [item.value for item in app.success],
        )
        self.assertTrue(
            any("APY defaults have been restored" in item.value for item in app.warning)
        )

    def test_rendered_stale_experimental_results_are_invalidated(self) -> None:
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
        self.assertTrue(is_experimental_infection_history_results(bundle))

        results_app = AppTest.from_file(str(Path("pages/3_Results.py")))
        results_app.session_state["config"] = self.config
        results_app.session_state["economics_config"] = self.economics_config
        results_app.session_state["results_bundle"] = bundle
        results_app.run(timeout=90)
        self.assertFalse(results_app.exception)
        self.assertIsNone(results_app.session_state["results_bundle"])
        self.assertTrue(any("not available in this SA Health version" in item.value for item in results_app.warning))

        econ_app = AppTest.from_file(str(Path("pages/4_Economics.py")))
        econ_app.session_state["config"] = self.config
        econ_app.session_state["economics_config"] = self.economics_config
        econ_app.session_state["results_bundle"] = bundle
        econ_app.run(timeout=90)
        self.assertFalse(econ_app.exception)
        self.assertIsNone(econ_app.session_state["results_bundle"])
        self.assertTrue(any("not available in this SA Health version" in item.value for item in econ_app.warning))

    def test_recent_remote_expected_value_ledger_is_rejected_by_positive_whitelist(self) -> None:
        bundle = {
            "metadata": {
                "scenarioLabel": "screenshot_recent_remote_fixture",
                "modelType": "expected_value",
                "analysisMethod": "expected_value",
                "naturalHistorySemantics": "continuous_markov_recent_remote",
                "analysisBasis": "",
            },
            "technical": {
                "interfaceConfig": {
                    "analysisMethod": "expected_value",
                    "naturalHistorySemantics": "continuous_markov_recent_remote",
                    "analysisBasis": "",
                },
                "eventLedger": {
                    "metadata": {
                        "modelType": "expected_value",
                        "naturalHistorySemantics": "continuous_markov_recent_remote",
                        "ltbiStateModel": "continuous_markov_recent_remote",
                        "analysisBasis": "",
                    },
                    "replicateTotals": pd.DataFrame(
                        [
                            {
                                "modelType": "expected_value",
                                "replicateId": 0,
                                "pairedReplicateId": 0,
                                "arm": "intervention",
                                "eventName": "active_tb_cases_prevented",
                                "value": 25.4,
                            }
                        ]
                    ),
                    "annualEvents": pd.DataFrame(),
                    "definitions": pd.DataFrame(),
                },
            },
        }
        economics = {
            "summaryRows": [
                {"metric": "incrementalCost", "discountProfile": "primary", "mean": -304182.0},
                {"metric": "dalysAverted", "discountProfile": "primary", "mean": 34.4955},
                {"metric": "activeTBCasesPrevented", "discountProfile": "primary", "mean": 25.4},
            ]
        }

        self.assertFalse(is_supported_sa_health_analysis_basis(bundle))

        app = AppTest.from_file(str(Path("pages/4_Economics.py")))
        app.session_state["config"] = self.config
        app.session_state["economics_config"] = self.economics_config
        app.session_state["results_bundle"] = bundle
        app.session_state["economics_results"] = economics
        app.session_state["results_stale"] = False
        app.session_state["economic_scenario_comparison"] = [{"Scenario": "Current assumptions"}]
        app.session_state["decision_scenario_comparison"] = [{"Scenario": "Current assumptions"}]

        app.run(timeout=90)

        self.assertFalse(app.exception)
        self.assertIsNone(app.session_state["results_bundle"])
        self.assertIsNone(app.session_state["economics_results"])
        self.assertIsNone(app.session_state["economic_scenario_comparison"])
        self.assertIsNone(app.session_state["decision_scenario_comparison"])
        self.assertTrue(any("not available in this SA Health version" in item.value for item in app.warning))
        rendered = "\n".join(
            str(getattr(item, "value", ""))
            for collection in [app.markdown, app.caption, app.warning, app.info]
            for item in collection
        )
        self.assertNotIn("34.4955", rendered)
        self.assertNotIn("25.4", rendered)
        self.assertNotIn("-304182", rendered)
        self.assertNotIn("Current assumptions", rendered)

        decision_app = AppTest.from_file(str(Path("pages/5_Decision_Analysis.py")))
        decision_app.session_state["config"] = self.config
        decision_app.session_state["economics_config"] = self.economics_config
        decision_app.session_state["results_bundle"] = bundle
        decision_app.session_state["economics_results"] = economics
        decision_app.session_state["results_stale"] = False
        decision_app.session_state["decision_scenario_comparison"] = [{"Strategy": "Current assumptions"}]
        decision_app.run(timeout=90)

        self.assertFalse(decision_app.exception)
        self.assertIsNone(decision_app.session_state["results_bundle"])
        self.assertIsNone(decision_app.session_state["economics_results"])
        self.assertIsNone(decision_app.session_state["decision_scenario_comparison"])

    def test_supported_deterministic_and_stochastic_compatibility_metadata_are_permitted(self) -> None:
        def bundle(model_type: str) -> dict[str, object]:
            return {
                "metadata": {
                    "modelType": model_type,
                    "analysisBasis": "sa_health_matlab_v9_compatibility_reference",
                    "naturalHistorySemantics": "matlab_v9_implicit_early_late",
                },
                "technical": {
                    "interfaceConfig": {
                        "analysisBasis": "sa_health_matlab_v9_compatibility_reference",
                        "naturalHistorySemantics": "matlab_v9_implicit_early_late",
                    },
                    "eventLedger": {
                        "metadata": {
                            "modelType": model_type,
                            "analysisBasis": "sa_health_matlab_v9_compatibility_reference",
                            "naturalHistorySemantics": "matlab_v9_implicit_early_late",
                        }
                    },
                },
            }

        self.assertTrue(is_supported_sa_health_analysis_basis(bundle("expected_value")))
        self.assertTrue(is_supported_sa_health_analysis_basis(bundle("agent_based")))
        self.assertFalse(is_supported_sa_health_analysis_basis({"metadata": {"modelType": "expected_value"}}))

    def test_supported_deterministic_run_survives_results_and_economics_navigation(self) -> None:
        config = build_unified_working_default_preset()["config"]
        config["analysisMethod"] = "expected_value"
        config["N"] = 50

        run_app = AppTest.from_file(str(Path("pages/2_Run_Model.py")), default_timeout=90)
        run_app.session_state["config"] = config
        run_app.session_state["economics_config"] = self.economics_config
        run_app.session_state["recent_ltbi_run_route"] = TECHNICAL_DEMONSTRATION_ROUTE
        run_app.run(timeout=90)
        next(item for item in run_app.button if item.label == "Run analysis").click().run(timeout=90)

        self.assertFalse(run_app.exception)
        bundle = run_app.session_state["results_bundle"]
        self.assertTrue(is_supported_sa_health_analysis_basis(bundle))
        self.assertFalse(run_app.session_state["results_stale"])
        self.assertEqual(bundle["metadata"]["modelType"], "expected_value")
        self.assertEqual(bundle["metadata"]["naturalHistorySemantics"], "matlab_v9_implicit_early_late")
        self.assertEqual(
            bundle["technical"]["eventLedger"]["metadata"]["naturalHistorySemantics"],
            "matlab_v9_implicit_early_late",
        )

        results_app = AppTest.from_file(str(Path("pages/3_Results.py")), default_timeout=90)
        results_app.session_state["config"] = run_app.session_state["config"]
        results_app.session_state["economics_config"] = self.economics_config
        results_app.session_state["results_bundle"] = bundle
        results_app.session_state["results_stale"] = False
        results_app.run(timeout=90)
        self.assertFalse(results_app.exception)
        self.assertIsNotNone(results_app.session_state["results_bundle"])
        self.assertFalse(any("not available in this SA Health version" in item.value for item in results_app.warning))

        econ_app = AppTest.from_file(str(Path("pages/4_Economics.py")), default_timeout=90)
        econ_app.session_state["config"] = run_app.session_state["config"]
        econ_app.session_state["economics_config"] = self.economics_config
        econ_app.session_state["results_bundle"] = bundle
        econ_app.session_state["results_stale"] = False
        econ_app.run(timeout=90)
        self.assertFalse(econ_app.exception)
        self.assertIsNotNone(econ_app.session_state["results_bundle"])
        self.assertFalse(any("not available in this SA Health version" in item.value for item in econ_app.warning))

    def test_setup_deterministic_selection_is_the_configuration_executed_by_run_analysis(self) -> None:
        setup_app = AppTest.from_file(str(Path("pages/0_Start.py")), default_timeout=90)
        setup_app.run(timeout=90)
        analysis_type = next(item for item in setup_app.radio if item.label == "Analysis type")
        analysis_type.set_value(
            "Quick deterministic preview - single expected-value calculation"
        ).run(timeout=90)

        selected_config = setup_app.session_state["config"]
        selected_config["N"] = 50
        selected_config.setdefault("scenario", {}).setdefault("eligible", {})["number"] = 50
        self.assertEqual(selected_config["analysisMethod"], "expected_value")
        self.assertNotIn("Repetitions", [item.label for item in setup_app.selectbox])
        self.assertNotIn("Random seed", [item.label for item in setup_app.number_input])
        self.assertTrue(
            any(
                "Quick deterministic preview - single expected-value calculation" in item.value
                for item in setup_app.success
            )
        )

        run_app = AppTest.from_file(str(Path("pages/2_Run_Model.py")), default_timeout=90)
        run_app.session_state["config"] = selected_config
        run_app.session_state["economics_config"] = setup_app.session_state["economics_config"]
        run_app.session_state["recent_ltbi_run_route"] = TECHNICAL_DEMONSTRATION_ROUTE
        run_app.run(timeout=90)
        self.assertNotIn("Preview analyses", [item.value for item in run_app.warning])
        expected_bundle = {
            "metadata": {
                "modelType": "expected_value",
                "analysisMethod": "expected_value",
                "analysisBasis": "sa_health_matlab_v9_compatibility_reference",
                "naturalHistorySemantics": "matlab_v9_implicit_early_late",
                "nReps": None,
                "seed": None,
            },
            "technical": {
                "interfaceConfig": selected_config,
                "eventLedger": {
                    "metadata": {
                        "modelType": "expected_value",
                        "analysisMethod": "expected_value",
                        "analysisBasis": "sa_health_matlab_v9_compatibility_reference",
                        "naturalHistorySemantics": "matlab_v9_implicit_early_late",
                    },
                    "annualEvents": [],
                },
            },
        }
        with patch(
            "adapters.python_apy_backend.PythonApyBackend.run_scenario_bundle",
            return_value=expected_bundle,
        ) as run_spy:
            next(item for item in run_app.button if item.label == "Run analysis").click().run(timeout=90)
        self.assertEqual(run_spy.call_args.args[0]["analysisMethod"], "expected_value")

        bundle = run_app.session_state["results_bundle"]
        self.assertEqual(bundle["metadata"]["modelType"], "expected_value")
        self.assertEqual(bundle["metadata"]["analysisMethod"], "expected_value")
        self.assertIsNone(bundle["metadata"].get("nReps"))
        self.assertIsNone(bundle["metadata"].get("seed"))

    def test_supported_small_stochastic_run_survives_results_and_economics_navigation(self) -> None:
        config = build_unified_working_default_preset()["config"]
        config.update({"analysisMethod": "agent_based", "N": 50, "nReps": 2, "seed": 1})

        run_app = AppTest.from_file(str(Path("pages/2_Run_Model.py")), default_timeout=90)
        run_app.session_state["config"] = config
        run_app.session_state["economics_config"] = self.economics_config
        run_app.session_state["recent_ltbi_run_route"] = TECHNICAL_DEMONSTRATION_ROUTE
        run_app.run(timeout=90)
        next(item for item in run_app.button if item.label == "Run analysis").click().run(timeout=90)

        self.assertFalse(run_app.exception)
        bundle = run_app.session_state["results_bundle"]
        self.assertTrue(is_supported_sa_health_analysis_basis(bundle))
        self.assertFalse(run_app.session_state["results_stale"])
        self.assertEqual(bundle["metadata"]["modelType"], "agent_based")
        self.assertEqual(bundle["metadata"]["nReps"], 2)
        self.assertEqual(bundle["metadata"]["seed"], 1)
        self.assertEqual(bundle["metadata"]["naturalHistorySemantics"], "matlab_v9_implicit_early_late")

        results_app = AppTest.from_file(str(Path("pages/3_Results.py")), default_timeout=90)
        results_app.session_state["config"] = run_app.session_state["config"]
        results_app.session_state["economics_config"] = self.economics_config
        results_app.session_state["results_bundle"] = bundle
        results_app.session_state["results_stale"] = False
        results_app.run(timeout=90)
        self.assertFalse(results_app.exception)
        self.assertIsNotNone(results_app.session_state["results_bundle"])

        econ_app = AppTest.from_file(str(Path("pages/4_Economics.py")), default_timeout=90)
        econ_app.session_state["config"] = run_app.session_state["config"]
        econ_app.session_state["economics_config"] = self.economics_config
        econ_app.session_state["results_bundle"] = bundle
        econ_app.session_state["results_stale"] = False
        econ_app.run(timeout=90)
        self.assertFalse(econ_app.exception)
        self.assertIsNotNone(econ_app.session_state["results_bundle"])

    def test_health_economics_widgets_apply_cost_edits_without_changing_health_outcomes(self) -> None:
        config = build_unified_working_default_preset()["config"]
        config.update({"analysisMethod": "agent_based", "N": 20, "nReps": 1, "seed": 1})
        config.setdefault("scenario", {}).setdefault("eligible", {})["number"] = 20
        run_config = prepare_run_config_for_recent_ltbi_route(
            config,
            selected_route=TECHNICAL_DEMONSTRATION_ROUTE,
        )
        bundle = build_results_bundle(run_replicates(run_config, keep_example_cohort=False))

        econ_app = AppTest.from_file(str(Path("pages/4_Economics.py")), default_timeout=90)
        econ_app.session_state["config"] = run_config
        econ_app.session_state["economics_config"] = self.economics_config
        econ_app.session_state["results_bundle"] = bundle
        econ_app.session_state["results_stale"] = False
        econ_app.run(timeout=90)
        before = econ_app.session_state["economics_results"]
        def summary_mean(payload: dict, metric: str) -> float:
            return next(
                float(row["mean"])
                for row in payload["summaryRows"]
                if row.get("metric") == metric and row.get("discountProfile") == "primary"
            )

        before_cost = summary_mean(before, "incrementalCost")
        before_dalys = summary_mean(before, "dalysAverted")
        before_tb = summary_mean(before, "activeTBCasesPrevented")

        setup = next(
            item for item in econ_app.number_input
            if item.label == "One-off programme setup cost (AUD)"
        )
        annual = next(
            item for item in econ_app.number_input
            if item.label == "Annual programme running cost (AUD per year)"
        )
        igra = next(
            item for item in econ_app.number_input
            if item.label == "IGRA screening test per person - Value used by model"
        )
        setup.set_value(100000.0)
        annual.set_value(10000.0)
        next(
            item for item in econ_app.number_input
            if item.label == "First year annual cost is incurred"
        ).set_value(1)
        igra.set_value(120.0)
        next(item for item in econ_app.button if item.label == "Recalculate economics").click().run(timeout=90)

        after = econ_app.session_state["economics_results"]
        self.assertGreater(summary_mean(after, "incrementalCost"), before_cost)
        self.assertEqual(summary_mean(after, "dalysAverted"), before_dalys)
        self.assertEqual(summary_mean(after, "activeTBCasesPrevented"), before_tb)
        controls = econ_app.session_state["health_econ_delivery_scenarios"]
        self.assertTrue(controls["includeAdditionalProgramCosts"])
        self.assertEqual(controls["standaloneSetupCost"], 100000.0)
        self.assertEqual(controls["standaloneAnnualRunningCost"], 10000.0)
        self.assertEqual(controls["annualCostFirstYear"], 1)
        igra_row = next(
            row
            for row in econ_app.session_state["health_econ_workspace"]["rows"]
            if row["assumptionId"] == "cost.test_igra"
        )
        self.assertEqual(igra_row["currentValue"], 120.0)
        self.assertEqual(igra_row["sourceCitation"], "User-defined")

        next(
            item for item in econ_app.button if item.label == "Restore SA Health economic defaults"
        ).click().run(timeout=90)
        restored = econ_app.session_state["economics_results"]
        self.assertAlmostEqual(summary_mean(restored, "incrementalCost"), before_cost, places=8)
        self.assertFalse(econ_app.session_state["health_econ_delivery_scenarios"]["includeAdditionalProgramCosts"])

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
