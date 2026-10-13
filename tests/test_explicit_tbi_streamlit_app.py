from __future__ import annotations

import subprocess
import unittest
from pathlib import Path

import pandas as pd
from streamlit.testing.v1 import AppTest

from app.explicit_tbi_workflow import (
    DETERMINISTIC_LABEL,
    STOCHASTIC_LABEL,
    active_tb_validation_for_state,
    build_explicit_workflow_config,
    default_workflow_state,
    discard_incompatible_result,
    health_signature,
    percent_to_proportion,
    proportion_to_percent,
    run_explicit_workflow_analysis,
    validate_workflow_state,
)
from engine.apy.calibration_policy import resolve_calibration_for_config
from engine.apy.explicit_recent_remote_runner import (
    EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
    TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION,
    explicit_recent_remote_event_ledger_unavailable,
    reconcile_deterministic_stochastic_explicit_intervention,
)
from engine.apy.regimen import (
    apply_regimen_overrides,
    default_regimen_library,
    get_regimen_from_library,
)


ROOT = Path(__file__).resolve().parents[1]
FROZEN_COMMIT = "03cc16e4a52a55e10019dab16bd4f599571c4b90"
APP_STATE_KEY = "explicit_tbi_workflow_state"
APP_RESULT_KEY = "explicit_tbi_results"


class ExplicitTBIInterventionReconciliationTests(unittest.TestCase):
    def _calibration_config_regimen(self, **overrides):
        state = default_workflow_state()
        state.update(
            {
                "population": 1200,
                "recentTBITarget": 0.20,
                "remoteOnlyTBITarget": 0.30,
                "screeningCoverage": 0.60,
                "testSensitivity": 1.0,
                "testSpecificity": 0.98,
                "pStartTPT": 1.0,
                "regimenPComplete": 0.90,
                "regimenADRstop": 0.0,
                "regimenEffFull": 0.90,
            }
        )
        state.update(overrides)
        config = build_explicit_workflow_config(state)
        calibration = resolve_calibration_for_config(config)
        library = default_regimen_library(config.get("partialShortCourseMode") or "threshold80")
        regimen = get_regimen_from_library(config["regimen"], library)
        regimen = apply_regimen_overrides(regimen, config)
        return state, config, calibration, regimen

    def test_intervention_reconciliation_covers_key_outputs(self) -> None:
        _state, config, calibration, regimen = self._calibration_config_regimen()
        reconciliation = reconcile_deterministic_stochastic_explicit_intervention(
            calibration,
            config=config,
            regimen=regimen,
            population_size=1200,
            seeds=(3101, 3102),
            horizon_years=20.0,
            absolute_floor_count=45.0,
        )

        required = {
            "eligiblePopulation",
            "screened",
            "truePositiveLatent",
            "falsePositive",
            "activeTBDetectedDuringWorkup",
            "tptStarts",
            "tptCompletions",
            "effectiveTPT",
            "comparatorActiveTB",
            "interventionActiveTB",
            "directActiveTBCasesAverted",
            "incidentActiveTBBeforeScreening",
            "incidentActiveTBAfterScreenBeforeEffectiveTPT",
            "subsequentIncidentActiveTB",
            "recentDirectActiveTBCasesAverted",
            "remoteOnlyDirectActiveTBCasesAverted",
        }
        self.assertTrue(reconciliation["allWithinTolerance"], reconciliation["rows"])
        self.assertTrue(required.issubset({row["metric"] for row in reconciliation["rows"]}))

    def test_intervention_reconciliation_scenarios_and_invariants(self) -> None:
        scenarios = [
            {"screeningCoverage": 0.0},
            {"screeningCoverage": 0.60},
            {"reinfectionPolicy": "recent_reinfection_does_not_reset_remote_clock_v1"},
            {"regimenEffFull": 1.0, "regimenPComplete": 1.0},
            {"regimenEffFull": 0.0},
        ]
        for overrides in scenarios:
            state, config, calibration, regimen = self._calibration_config_regimen(**overrides)
            result = run_explicit_workflow_analysis(state)
            row = result["deterministicResult"]["byHorizon"]["20"]
            self.assertLessEqual(
                row["directActiveTBCasesAvertedProportion"],
                row["comparatorActiveTBProportion"],
            )
            if state["screeningCoverage"] == 0.0:
                self.assertEqual(row["directActiveTBCasesAvertedProportion"], 0.0)
                self.assertEqual(
                    row["comparatorActiveTBProportion"],
                    row["interventionActiveTBProportion"],
                )
            if state["regimenEffFull"] == 0.0:
                self.assertEqual(row["directActiveTBCasesAvertedProportion"], 0.0)
            if state["regimenEffFull"] == 1.0 and state["regimenPComplete"] == 1.0:
                self.assertGreater(row["incidentActiveTBBeforeScreeningProportion"], 0.0)
                self.assertLessEqual(
                    row["directActiveTBCasesAvertedProportion"],
                    row["comparatorActiveTBProportion"]
                    - row["incidentActiveTBBeforeScreeningProportion"]
                    + 1e-12,
                )
            self.assertEqual(
                config["explicitRecentRemoteTBI"]["tptEffectPolicyIdentifier"],
                TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION,
            )
            self.assertFalse(
                explicit_recent_remote_event_ledger_unavailable(
                    config["explicitRecentRemoteTBI"]
                )["validation"]["isValid"]
            )
        _state, config, calibration, regimen = self._calibration_config_regimen()
        mortality = reconcile_deterministic_stochastic_explicit_intervention(
            calibration,
            config=config,
            regimen=regimen,
            population_size=900,
            seeds=(4101, 4102),
            horizon_years=20.0,
            competing_mortality_hazard=0.02,
            absolute_floor_count=45.0,
        )
        self.assertTrue(mortality["allWithinTolerance"], mortality["rows"])

    def test_common_random_pairing_and_false_positive_no_benefit(self) -> None:
        state = default_workflow_state()
        state.update(
            {
                "population": 150,
                "analysisMethod": STOCHASTIC_LABEL,
                "stochasticRepetitions": 1,
                "recentTBITarget": 0.20,
                "remoteOnlyTBITarget": 0.25,
                "screeningCoverage": 1.0,
                "testSensitivity": 1.0,
                "testSpecificity": 0.50,
                "pStartTPT": 1.0,
                "regimenPComplete": 1.0,
                "regimenADRstop": 0.0,
                "regimenEffFull": 1.0,
            }
        )
        result = run_explicit_workflow_analysis(state)["runnerResult"]
        replicate = result["explicitEpidemiology"]["replicates"][0]
        comparator_cases = replicate["comparator"]["totalActiveTBOverFollowUp"]
        intervention_cases = replicate["intervention"]["totalActiveTBOverFollowUp"]
        averted = replicate["intervention"][
            "directActiveTBCasesAvertedAmongModelledPopulation"
        ]

        self.assertEqual(comparator_cases, intervention_cases + averted)
        self.assertEqual(
            replicate["intervention"]["uninfectedFalsePositiveEpidemiologicalBenefit"],
            0,
        )
        self.assertIn("share baseline states", replicate["diagnosticMessages"][0])


class ExplicitTBIWorkflowTests(unittest.TestCase):
    def test_new_app_requires_explicit_pathway_and_rejects_other_results(self) -> None:
        state = default_workflow_state()
        config = build_explicit_workflow_config(state)

        self.assertEqual(config["analysisPathway"], EXPLICIT_RECENT_REMOTE_PATHWAY_ID)
        self.assertEqual(
            config["explicitRecentRemoteTBI"]["pathwayIdentifier"],
            EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
        )
        self.assertIsNone(discard_incompatible_result({"analysisPathway": "matlab_v9"}))
        self.assertIsNone(
            discard_incompatible_result(
                {"analysisPathway": "continuous_markov_recent_remote"}
            )
        )

    def test_frozen_app_cannot_navigate_to_explicit_workflow(self) -> None:
        text = (ROOT / "streamlit_app.py").read_text(encoding="utf-8")

        self.assertNotIn("explicit_tbi_app.py", text)
        self.assertNotIn("Explicit recent/remote TBI development workflow", text)
        self.assertNotIn("continuous_markov_recent_remote", (ROOT / "explicit_tbi_app.py").read_text(encoding="utf-8"))

    def test_percent_roundtrip_invalid_sums_and_user_defined_sources(self) -> None:
        value = 2.345678
        self.assertAlmostEqual(
            proportion_to_percent(percent_to_proportion(value)),
            value,
            places=12,
        )
        state = default_workflow_state()
        state["recentTBITarget"] = 0.95
        state["remoteOnlyTBITarget"] = 0.10
        errors = validate_workflow_state(state)

        self.assertTrue(errors)
        self.assertEqual(state["recentTBITarget"], 0.95)
        self.assertEqual(state["remoteOnlyTBITarget"], 0.10)

        edited = default_workflow_state()
        edited["recentTBITarget"] = 0.03
        from app.explicit_tbi_workflow import source_for_value

        self.assertEqual(source_for_value(edited, "recentTBITarget"), "User-defined")
        self.assertNotEqual(source_for_value(edited, "remoteOnlyTBITarget"), "User-defined")

    def test_active_tb_denominator_changes_validation_only(self) -> None:
        state = default_workflow_state()
        state.update(
            {
                "population": 200,
                "activeTBObservationEnabled": True,
                "activeTBObservedCases": 10.0,
                "activeTBDenominator": 1500.0,
                "activeTBObservationDurationYears": 2.0,
            }
        )
        result = run_explicit_workflow_analysis(state)
        summary_before = pd.DataFrame(result["summaryRows"]).copy()
        config_hash = result["configHash"]
        signature = health_signature(state)
        validation_before = active_tb_validation_for_state(state)

        state["activeTBDenominator"] = 3000.0
        validation_after = active_tb_validation_for_state(state)

        self.assertEqual(signature, health_signature(state))
        self.assertEqual(config_hash, result["configHash"])
        pd.testing.assert_frame_equal(summary_before, pd.DataFrame(result["summaryRows"]))
        self.assertNotEqual(
            validation_before["expectedCases"],
            validation_after["expectedCases"],
        )

    def test_age_specific_prevalence_reconciles_to_targets(self) -> None:
        state = default_workflow_state()
        result = run_explicit_workflow_analysis(state)
        rows = pd.DataFrame(result["ageRows"])
        recent = float(
            (rows["Population share"] * rows["Recent infection prevalence"]).sum()
        )
        remote = float(
            (rows["Population share"] * rows["Remote-only infection prevalence"]).sum()
        )

        self.assertAlmostEqual(recent, state["recentTBITarget"], places=7)
        self.assertAlmostEqual(remote, state["remoteOnlyTBITarget"], places=7)

    def test_frozen_release_files_and_headline_pointer_unchanged(self) -> None:
        changed = subprocess.check_output(
            [
                "git",
                "diff",
                "--name-only",
                "--",
                "reports/sa_health_apy_health_economic_working_report_v1",
                "validation/matlab_reference",
                "abm",
            ],
            cwd=ROOT,
            text=True,
        ).splitlines()
        branch = subprocess.check_output(
            ["git", "rev-parse", "release/sa-health-apy-he-v1.0.0"],
            cwd=ROOT,
            text=True,
        ).strip()
        tag = subprocess.check_output(
            ["git", "rev-parse", "sa-health-apy-he-v1.0.0^{}"],
            cwd=ROOT,
            text=True,
        ).strip()

        self.assertEqual(changed, [])
        self.assertEqual(branch, FROZEN_COMMIT)
        self.assertEqual(tag, FROZEN_COMMIT)


class ExplicitTBIStreamlitRenderedTests(unittest.TestCase):
    def _app(self, state: dict | None = None) -> AppTest:
        app = AppTest.from_file(str(ROOT / "explicit_tbi_app.py"), default_timeout=180)
        if state is not None:
            app.session_state[APP_STATE_KEY] = state
        return app

    def test_rendered_setup_inputs_and_guardrails(self) -> None:
        state = default_workflow_state()
        state["population"] = 100
        app = self._app(state)
        app.run(timeout=180)

        self.assertFalse(app.exception)
        labels = [item.label for item in app.number_input]
        self.assertIn("Recently infected within the past 5 years", labels)
        self.assertIn(
            "Remote infection only — infected 5 or more years ago and not reinfected within the past 5 years",
            labels,
        )
        text = _rendered_text(app)
        self.assertIn("development workflow", text)
        self.assertNotIn("Health Economics", text)
        self.assertNotIn("ICER", text)
        self.assertNotIn("Rising", text)
        self.assertEqual(
            [button.label for button in app.button].count("Restore working defaults"),
            1,
        )

        next(
            item
            for item in app.number_input
            if item.label == "Recently infected within the past 5 years"
        ).set_value(95.0).run(timeout=180)
        next(
            item
            for item in app.number_input
            if item.label.startswith("Remote infection only")
        ).set_value(10.0).run(timeout=180)
        self.assertTrue(app.error)
        self.assertIn("must not exceed 100%", app.error[0].value)

        next(button for button in app.button if button.label == "Restore working defaults").click().run(timeout=180)
        self.assertEqual(
            app.session_state[APP_STATE_KEY],
            default_workflow_state(),
        )

    def test_rendered_deterministic_and_stochastic_modes(self) -> None:
        state = default_workflow_state()
        state["population"] = 80
        app = self._app(state)
        app.run(timeout=180)

        self.assertNotIn("Repetitions", [item.label for item in app.number_input])
        next(button for button in app.button if button.label == "Continue to Run analysis").click().run(timeout=180)
        next(button for button in app.button if button.label == "Run explicit TBI analysis").click().run(timeout=180)
        self.assertFalse(app.exception)
        self.assertIn(APP_RESULT_KEY, app.session_state)
        result = app.session_state[APP_RESULT_KEY]
        self.assertEqual(result["analysisPathway"], EXPLICIT_RECENT_REMOTE_PATHWAY_ID)
        self.assertEqual(result["mode"], "deterministic")

        next(button for button in app.button if button.label == "Open Results").click().run(timeout=180)
        result_text = _rendered_text(app)
        self.assertIn("Comparator active TB", result_text)
        self.assertIn("Direct active TB cases averted among the modelled population", result_text)
        self.assertIn("N/A", result_text)

        state = default_workflow_state()
        state.update(
            {
                "population": 60,
                "analysisMethod": STOCHASTIC_LABEL,
                "stochasticRepetitions": 2,
                "seed": 99,
            }
        )
        app = self._app(state)
        app.run(timeout=180)
        self.assertIn("Repetitions", [item.label for item in app.number_input])
        next(button for button in app.button if button.label == "Continue to Run analysis").click().run(timeout=180)
        next(button for button in app.button if button.label == "Run explicit TBI analysis").click().run(timeout=180)
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state[APP_RESULT_KEY]["mode"], "stochastic")
        next(button for button in app.button if button.label == "Open Results").click().run(timeout=180)
        self.assertIn("Simulation intervals use", _rendered_text(app))

    def test_rendered_material_edit_marks_results_stale(self) -> None:
        state = default_workflow_state()
        state["population"] = 80
        result = run_explicit_workflow_analysis(state)
        app = self._app(state)
        app.session_state[APP_RESULT_KEY] = result
        app.run(timeout=180)

        next(
            item
            for item in app.number_input
            if item.label == "Screening coverage"
        ).set_value(40.0).run(timeout=180)
        next(
            radio for radio in app.radio if radio.label == "Workflow"
        ).set_value("Results").run(timeout=180)

        self.assertIn("stale", _rendered_text(app).lower())


def _rendered_text(app: AppTest) -> str:
    parts = []
    for collection_name in (
        "title",
        "header",
        "subheader",
        "markdown",
        "caption",
        "info",
        "warning",
        "error",
        "success",
    ):
        for item in getattr(app, collection_name):
            parts.append(str(getattr(item, "value", "")))
    for frame in app.dataframe:
        parts.append(frame.value.to_string())
    return "\n".join(parts)


if __name__ == "__main__":
    unittest.main()

