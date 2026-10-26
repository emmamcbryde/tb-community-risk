from __future__ import annotations

from copy import deepcopy
from decimal import Decimal
from functools import lru_cache
import json
import math
from pathlib import Path
import subprocess
import types
import unittest
from unittest.mock import patch

import numpy as np
import engine.apy.explicit_recent_remote_runner as runner_module
from streamlit.testing.v1 import AppTest

from app.explicit_tbi_formatting import format_number, format_reconciled_components, format_multiplier
from app.explicit_tbi_session import DISPATCH_COUNT_KEY, PAGE_KEY, RESULT_KEY, SCHEMA, SCHEMA_KEY, STATE_KEY, STALE_KEY
from app.explicit_tbi_workflow import (
    PROGRESSION_RISK_LABELS, WORKFLOW_PAGES, active_tb_validation_for_state,
    build_explicit_workflow_config, build_regimen_for_state, default_workflow_state,
    health_signature, progression_risk_for_state, run_explicit_workflow_analysis,
    technical_results_json, validate_projection_accounting,
)
from engine.apy.calibration_policy import resolve_calibration_for_config
from engine.apy.config import normalise_config
from engine.apy.explicit_recent_remote_progression import (
    REINFECTION_POLICY_NO_RESET_CLOCK, REINFECTION_POLICY_RESET_CLOCK,
    STATE_BASELINE_ACTIVE_TB,
    build_central_time_since_infection_progression_curve,
    competing_risk_progression_probability_for_curve,
    time_since_curve_conditional_future_progression_probability,
    time_since_curve_incremental_cumulative_hazard,
)
from engine.apy.explicit_recent_remote_runner import (
    apply_explicit_recent_remote_population_assignment,
    deterministic_explicit_recent_remote_expectations,
    deterministic_explicit_recent_remote_intervention_expectations,
    explicit_recent_remote_runner_config_hash,
    reconcile_deterministic_stochastic_explicit_intervention,
    sample_future_active_tb_time_from_curve,
    stochastic_explicit_recent_remote_natural_history,
)
from engine.apy.population_tbi_progression_risk import (
    ADJUSTMENT_FIELD, CONTRACT_VERSION, build_progression_risk_adjustment,
    validate_progression_risk_adjustment,
)

ROOT = Path(__file__).resolve().parents[1]
START = "136de0556d0ce6ce0ce9cb76b04e22efa5e2cffc"
OUTPUT = ROOT / "abm/output/milestone3g"


def scenario(mode="risk_enriched", **updates):
    state = default_workflow_state()
    state.update(progressionRiskMode=mode, **updates)
    config = build_explicit_workflow_config(state)
    calibration = resolve_calibration_for_config(config)
    return state, config, calibration, build_regimen_for_state(state, normalise_config(config))


@lru_cache(maxsize=1)
def working_result():
    result = run_explicit_workflow_analysis(default_workflow_state())
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "working-default-full-results.json").write_text(technical_results_json(result), encoding="utf-8")
    return result


def app_for(*, state=None, result=None, page="Set up"):
    app = AppTest.from_file(str(ROOT / "explicit_tbi_ui/app.py"), default_timeout=240)
    app.session_state[SCHEMA_KEY] = deepcopy(SCHEMA)
    app.session_state[STATE_KEY] = deepcopy(state or default_workflow_state())
    app.session_state[RESULT_KEY] = deepcopy(result)
    app.session_state[PAGE_KEY] = page
    return app


def radio(app, label):
    return next(item for item in app.radio if item.label == label)


def button(app, label):
    return next(item for item in app.button if item.label == label)


def text(app):
    return "\n".join(str(item.value) for kind in ("title", "subheader", "caption", "info", "warning", "error", "markdown") for item in getattr(app, kind))


class ProgressionRiskContractTests(unittest.TestCase):
    def test_standard_resolves_exactly_one(self):
        self.assertEqual(build_progression_risk_adjustment("standard")["multiplier"], 1.0)

    def test_risk_enriched_resolves_exactly_eight_provisionally(self):
        contract = build_progression_risk_adjustment("risk_enriched")
        self.assertEqual(contract["multiplier"], 8.0)
        self.assertEqual(contract["contractVersion"], CONTRACT_VERSION)
        self.assertEqual(contract["source"], "Working scenario assumption")
        self.assertIn("provisional", contract["reviewStatus"])
        self.assertEqual(progression_risk_for_state(default_workflow_state())["multiplier"], 8.0)

    def test_custom_retained_exactly_and_modes_are_validated(self):
        value = 3.141592653589793
        self.assertEqual(build_progression_risk_adjustment("custom", multiplier=value)["multiplier"], value)
        with self.assertRaises(ValueError):
            build_progression_risk_adjustment("standard", multiplier=8.0)
        with self.assertRaises(ValueError):
            build_progression_risk_adjustment("unknown")

    def test_invalid_custom_values_rejected_without_clipping(self):
        for value in (None, 0.0, -1.0, math.nan, math.inf, -math.inf, True, "8", [], {}, np.array([8.0])):
            with self.subTest(value=value), self.assertRaises(ValueError):
                build_progression_risk_adjustment("custom", multiplier=value)
        contract = build_progression_risk_adjustment("custom", multiplier=2.0)
        del contract["multiplier"]
        with self.assertRaises(ValueError):
            validate_progression_risk_adjustment(contract)

    def test_mode_and_resolved_multiplier_in_hash_and_staleness_signature(self):
        from engine.apy.runner import _calibration_key_payload
        signatures, hashes, cache_keys = [], [], []
        for mode, custom in (("standard", 8.0), ("risk_enriched", 8.0), ("custom", 1.0), ("custom", 1.123456789012345)):
            state, config, _cal, _reg = scenario(mode, customProgressionRiskMultiplier=custom)
            contract = config["explicitRecentRemoteTBI"][ADJUSTMENT_FIELD]
            signatures.append(health_signature(state))
            hashes.append(explicit_recent_remote_runner_config_hash(config["explicitRecentRemoteTBI"]))
            cache_keys.append(json.dumps(_calibration_key_payload(config), sort_keys=True))
            self.assertIn(mode, signatures[-1])
            self.assertIn("multiplier", signatures[-1])
            if mode == "custom":
                self.assertEqual(contract["multiplier"], custom)
        self.assertEqual(len(set(signatures)), 4)
        self.assertEqual(len(set(hashes)), 4)
        self.assertEqual(len(set(cache_keys)), 4)


class ProgressionHazardMathematicsTests(unittest.TestCase):
    def test_curve_and_anchors_unchanged_from_starting_commit(self):
        source = subprocess.check_output(["git", "show", f"{START}:engine/apy/explicit_recent_remote_progression.py"], cwd=ROOT, text=True, encoding="utf-8")
        baseline = types.ModuleType("milestone3f_progression_baseline")
        exec(compile(source, "<milestone3f-baseline>", "exec"), baseline.__dict__)
        for policy in (REINFECTION_POLICY_RESET_CLOCK, REINFECTION_POLICY_NO_RESET_CLOCK):
            self.assertEqual(baseline.build_central_time_since_infection_progression_curve(reinfection_policy=policy), build_central_time_since_infection_progression_curve(reinfection_policy=policy))

    def test_probability_is_hazard_scaled_not_probability_multiplied(self):
        curve = build_central_time_since_infection_progression_curve()
        for clock in (0.0, 0.5, 2.5, 4.9, 5.0, 10.0, 24.9, 25.0, 70.0):
            for horizon in (0.0, 1.0, 5.0, 20.0):
                increment = time_since_curve_incremental_cumulative_hazard(curve, time_since_infection_at_baseline=clock, horizon_years=horizon)
                probability = time_since_curve_conditional_future_progression_probability(curve, time_since_infection_at_baseline=clock, horizon_years=horizon, progression_hazard_multiplier=8.0)
                self.assertEqual(probability, 1.0 - math.exp(-8.0 * increment))
                if increment > 0:
                    self.assertLess(probability, 8.0 * (1.0 - math.exp(-increment)))

    def test_stochastic_inverse_uses_same_scaled_increment_and_uniform_draw(self):
        curve = build_central_time_since_infection_progression_curve()
        class FixedUniform:
            def random(self):
                return 0.4
        for clock in (0.0, 2.0, 5.0, 24.0, 40.0):
            times = []
            for multiplier in (1.0, 8.0):
                event = sample_future_active_tb_time_from_curve(curve, time_since_infection_at_baseline=clock, progression_hazard_multiplier=multiplier, rng=FixedUniform())
                increment = time_since_curve_incremental_cumulative_hazard(curve, time_since_infection_at_baseline=clock, horizon_years=event)
                self.assertAlmostEqual(multiplier * increment, -math.log1p(-0.4), delta=9e-12)
                times.append(event)
            self.assertLess(times[1], times[0])

    def test_mortality_hazard_is_not_scaled_and_survival_route_agrees(self):
        curve = build_central_time_since_infection_progression_curve()
        base = curve["postFinalAnchorHazard"]
        death, horizon = 0.02, 20.0
        for multiplier in (1.0, 8.0):
            expected = multiplier * base / (multiplier * base + death) * -math.expm1(-(multiplier * base + death) * horizon)
            analytic = competing_risk_progression_probability_for_curve(curve, time_since_infection_at_baseline=30.0, horizon_years=horizon, competing_mortality_hazard=death, progression_hazard_multiplier=multiplier)
            external = competing_risk_progression_probability_for_curve(curve, time_since_infection_at_baseline=30.0, horizon_years=horizon, competing_survival_curve=lambda time: math.exp(-death * time), integration_steps=2048, progression_hazard_multiplier=multiplier)
            self.assertAlmostEqual(analytic["cumulativeIncidence"], expected, places=13)
            self.assertAlmostEqual(external["cumulativeIncidence"], expected, delta=1e-7)
            callable_route = runner_module._progression_probability_for_clock(
                curve, 30.0, horizon, 40.0, competing_survival_curve=None,
                competing_mortality_hazard=lambda time, age: death,
                integration_steps=512, progression_hazard_multiplier=multiplier,
            )
            self.assertAlmostEqual(callable_route, expected, delta=1e-7)

    def test_tbi_age_calibration_and_sampled_infection_history_are_identical(self):
        _s1, _c1, cal1, _r1 = scenario("standard", baselineActiveTBProportion=0.01)
        _s8, _c8, cal8, _r8 = scenario("risk_enriched", baselineActiveTBProportion=0.01)
        self.assertEqual(cal1["explicitRecentRemoteCalibration"], cal8["explicitRecentRemoteCalibration"])
        self.assertEqual(cal1["progressionCurve"], cal8["progressionCurve"])
        for key in ("exactAgeValues", "exactAgeProb"):
            np.testing.assert_array_equal(cal1["parameters"][key], cal8["parameters"][key])
        ages = np.tile([10.0, 35.0, 70.0], 500)
        outputs = [apply_explicit_recent_remote_population_assignment({"ageYears": ages}, cal, np.random.default_rng(223344)) for cal in (cal1, cal8)]
        for key in ("explicitTBIState", "priorRemoteExposure", "priorRemotePlusRecent", "baselineActiveTBAtBaseline", "timeSinceMostRecentInfection", "timeSinceRemoteInfectionNoReset", "progressionClockYears", "pInfection", "diseaseMultiplier"):
            np.testing.assert_array_equal(outputs[0][0][key], outputs[1][0][key])
        for output in outputs:
            exempt = np.isin(output[0]["explicitTBIState"], ["uninfected", STATE_BASELINE_ACTIVE_TB])
            self.assertTrue(np.all(np.isinf(output[4][exempt])))
            np.testing.assert_array_equal(output[0]["diseaseMultiplier"], np.ones(len(ages)))
        self.assertTrue(np.all(outputs[1][4] <= outputs[0][4]))

    def test_reset_and_no_reset_intervention_match_natural_history(self):
        for policy in (REINFECTION_POLICY_RESET_CLOCK, REINFECTION_POLICY_NO_RESET_CLOCK):
            for mode in ("standard", "risk_enriched"):
                _state, cfg, cal, reg = scenario(mode, reinfectionPolicy=policy, recentTBITarget=0.2, remoteOnlyTBITarget=0.3, baselineActiveTBProportion=0.05)
                det = deterministic_explicit_recent_remote_intervention_expectations(cal, config=cfg, regimen=reg, infection_quadrature_points=4, screen_quadrature_points=4)
                natural = deterministic_explicit_recent_remote_expectations(cal, horizons_years=(20.0,), quadrature_points=4)
                self.assertAlmostEqual(det["byHorizon"]["20"]["comparatorActiveTBProportion"], natural["incidentActiveTBByHorizon"]["20"]["incidentActiveTBProportion"], places=12)
                self.assertTrue(det["accountingValidation"]["isValid"])

    def test_zero_coverage_and_zero_efficacy_prevent_no_cases(self):
        for mode in ("standard", "risk_enriched"):
            for changes in ({"screeningCoverage": 0.0}, {"regimenEffFull": 0.0}):
                _state, cfg, cal, reg = scenario(mode, **changes)
                det = deterministic_explicit_recent_remote_intervention_expectations(cal, config=cfg, regimen=reg, infection_quadrature_points=4, screen_quadrature_points=4)
                horizon = det["byHorizon"]["20"]
                self.assertEqual(horizon["directActiveTBCasesAvertedProportion"], 0.0)
                self.assertEqual(horizon["interventionActiveTBProportion"], horizon["comparatorActiveTBProportion"])

    def test_working_output_non_linear_enrichment_and_group_accounting(self):
        result = working_result()
        standard = result["progressionRiskSensitivity"]["standardProjection"]["byHorizon"]["20"]
        enriched = result["analyticalProjection"]["byHorizon"]["20"]
        self.assertGreater(enriched["comparatorActiveTBProportion"], standard["comparatorActiveTBProportion"])
        self.assertLess(enriched["comparatorActiveTBProportion"], 8.0 * standard["comparatorActiveTBProportion"])
        self.assertNotEqual(enriched["comparatorActiveTBProportion"] * 1500, 20.0)
        rows = result["comparatorDecomposition"]
        self.assertAlmostEqual(rows[0]["numberOfPeople"] + rows[1]["numberOfPeople"], 150.0, places=6)
        self.assertAlmostEqual(rows[0]["expectedComparatorActiveTBCases"] + rows[1]["expectedComparatorActiveTBCases"], rows[-1]["expectedComparatorActiveTBCases"], places=12)
        self.assertEqual(rows[2]["expectedComparatorActiveTBCases"], 0.0)
        self.assertIsNone(rows[3]["expectedComparatorActiveTBCases"])
        self.assertTrue(result["projectionAccountingValidation"]["isValid"])
        for horizon in (standard, enriched):
            self.assertAlmostEqual(horizon["interventionActiveTBProportion"] + horizon["directActiveTBCasesAvertedProportion"], horizon["comparatorActiveTBProportion"], places=14)
        selected_config = deepcopy(result["config"])
        standard_config = deepcopy(result["progressionRiskSensitivity"]["standardConfiguration"])
        selected_config["explicitRecentRemoteTBI"].pop(ADJUSTMENT_FIELD)
        standard_config["explicitRecentRemoteTBI"].pop(ADJUSTMENT_FIELD)
        self.assertEqual(selected_config, standard_config)

    def test_accounting_fails_on_corrupt_components(self):
        result = working_result()
        for index, field in ((0, "expectedComparatorActiveTBCases"), (2, "expectedComparatorActiveTBCases"), (3, "expectedComparatorActiveTBCases"), (4, "numberOfPeople")):
            rows = deepcopy(result["comparatorDecomposition"])
            rows[index][field] = 999.0
            self.assertFalse(validate_projection_accounting(result["config"], result["analyticalProjection"], rows)["isValid"])

    def test_large_valid_custom_multiplier_is_not_clipped_or_rejected_for_roundoff(self):
        state = default_workflow_state()
        state.update(progressionRiskMode="custom", customProgressionRiskMultiplier=1_000_000.0)
        result = run_explicit_workflow_analysis(state)
        self.assertEqual(result[ADJUSTMENT_FIELD]["multiplier"], 1_000_000.0)
        self.assertTrue(result["projectionAccountingValidation"]["isValid"])
        self.assertLessEqual(result["comparatorDecomposition"][-1]["expectedComparatorActiveTBCases"], 150.0 + 1e-6)
        for row in result["comparatorDecomposition"][:2]:
            self.assertLessEqual(row["averageProgressionProbability"], 1.0 + 1e-12)

    def test_observations_cannot_change_multiplier_curve_or_projection(self):
        state, cfg, cal, _reg = scenario()
        before = deterministic_explicit_recent_remote_expectations(cal, horizons_years=(20.0,), quadrature_points=4)
        signature = health_signature(state)
        for cases, denominator, duration in ((20, 1500, 20), (10, 3000, 2), (40, 2000, 5)):
            changed = {**state, "activeTBObservationEnabled": True, "activeTBObservedCases": cases, "activeTBDenominator": denominator, "activeTBObservationDurationYears": duration}
            self.assertEqual(build_explicit_workflow_config(changed), cfg)
            self.assertEqual(health_signature(changed), signature)
            validation = active_tb_validation_for_state(changed)
            self.assertEqual(validation["status"], "validation_only")
            self.assertEqual(validation["observedCases"], cases)
            self.assertEqual(progression_risk_for_state(changed)["multiplier"], 8.0)
        after = deterministic_explicit_recent_remote_expectations(cal, horizons_years=(20.0,), quadrature_points=4)
        self.assertEqual(before, after)

    def test_predeclared_monte_carlo_reconciliation_at_one_and_eight(self):
        # Declared before seeing results: 20 seeds, N=1500, 4 MCSE, 0.5 case floor.
        reports = []
        for mode in ("standard", "risk_enriched"):
            _state, cfg, cal, reg = scenario(mode)
            report = reconcile_deterministic_stochastic_explicit_intervention(cal, config=cfg, regimen=reg, population_size=1500, seeds=tuple(range(8800, 8820)), infection_quadrature_points=16, screen_quadrature_points=12, absolute_floor_count=0.5, z_score=4.0)
            reports.append({"mode": mode, "rows": report["rows"], "allWithinTolerance": report["allWithinTolerance"], "seedCount": 20})
            self.assertTrue(report["allWithinTolerance"], report["rows"])
        OUTPUT.mkdir(parents=True, exist_ok=True)
        (OUTPUT / "progression-risk-mcse-reconciliation.json").write_text(json.dumps(reports, indent=2), encoding="utf-8")

    def test_no_reset_and_mortality_routes_reconcile_without_changing_death_draws(self):
        _state, cfg, cal, reg = scenario("risk_enriched", reinfectionPolicy=REINFECTION_POLICY_NO_RESET_CLOCK)
        # Same predeclared rule; mortality enters independently of progression.
        report = reconcile_deterministic_stochastic_explicit_intervention(cal, config=cfg, regimen=reg, population_size=1500, seeds=tuple(range(8900, 8920)), infection_quadrature_points=4, screen_quadrature_points=4, competing_mortality_hazard=0.02, absolute_floor_count=0.5, z_score=4.0)
        # This is an epidemiological progression check. The pre-existing helper's
        # ancillary cascade approximations under mortality are outside this UI,
        # which still supplies no mortality data. Do not claim all such counts.
        core = {"comparatorActiveTB", "interventionActiveTB", "directActiveTBCasesAverted"}
        self.assertTrue(all(row["withinTolerance"] for row in report["rows"] if row["metric"] in core), report["rows"])
        OUTPUT.mkdir(parents=True, exist_ok=True)
        (OUTPUT / "no-reset-mortality-mcse-reconciliation.json").write_text(json.dumps(report["rows"], indent=2), encoding="utf-8")
        captured = []
        actual_draw = runner_module._draw_death_times_from_hazard
        def capture(*args, **kwargs):
            result = actual_draw(*args, **kwargs)
            captured.append(result.copy())
            return result
        with patch.object(runner_module, "_draw_death_times_from_hazard", side_effect=capture):
            for mode in ("standard", "risk_enriched"):
                _state, _cfg, calibration, _reg = scenario(mode, reinfectionPolicy=REINFECTION_POLICY_NO_RESET_CLOCK)
                stochastic_explicit_recent_remote_natural_history(calibration, population_size=400, horizons_years=(20.0,), seed=76123, competing_mortality_hazard=0.02)
        np.testing.assert_array_equal(captured[0], captured[1])

    def test_standard_projection_matches_starting_runner(self):
        source = subprocess.check_output(["git", "show", f"{START}:engine/apy/explicit_recent_remote_runner.py"], cwd=ROOT, text=True, encoding="utf-8")
        baseline = types.ModuleType("milestone3f_runner_baseline")
        exec(compile(source, "<milestone3f-runner-baseline>", "exec"), baseline.__dict__)
        _state, cfg, cal, reg = scenario("standard")
        arguments = {"config": cfg, "regimen": reg, "infection_quadrature_points": 4, "screen_quadrature_points": 4}
        old = baseline.deterministic_explicit_recent_remote_intervention_expectations(cal, **arguments)
        current = deterministic_explicit_recent_remote_intervention_expectations(cal, **arguments)
        self.assertEqual(old["byHorizon"], current["byHorizon"])
        self.assertEqual(old["stateProportions"], current["stateProportions"])

    def test_nonzero_baseline_active_tb_is_excluded_under_both_adjustments(self):
        for mode in ("standard", "risk_enriched"):
            _state, cfg, cal, reg = scenario(mode, baselineActiveTBProportion=0.2)
            projection = deterministic_explicit_recent_remote_intervention_expectations(cal, config=cfg, regimen=reg, infection_quadrature_points=4, screen_quadrature_points=4)
            from app.explicit_tbi_workflow import comparator_progression_decomposition
            rows = comparator_progression_decomposition(projection, 1500, 20.0)
            self.assertAlmostEqual(rows[3]["numberOfPeople"], 300.0)
            self.assertIsNone(rows[3]["expectedComparatorActiveTBCases"])
            self.assertEqual(rows[2]["expectedComparatorActiveTBCases"], 0.0)
            self.assertTrue(validate_projection_accounting(cfg, projection, rows)["isValid"])


class ProgressionRiskRenderedTests(unittest.TestCase):
    def test_legacy_economics_and_costable_ledger_are_unavailable(self):
        app = AppTest.from_file(str(ROOT / "pages/4_Economics.py"), default_timeout=180)
        app.session_state[SCHEMA_KEY] = deepcopy(SCHEMA)
        app.session_state[RESULT_KEY] = deepcopy(working_result())
        app.session_state["economics_results"] = {"stale": 12345}
        app.run()
        self.assertFalse(app.exception)
        self.assertIn("no validated costable event ledger", text(app))
        self.assertFalse(app.metric)
        self.assertFalse(app.dataframe)
        self.assertNotIn("economics_results", app.session_state)
        from engine.apy.event_ledger_economics import run_event_ledger_health_economics
        with self.assertRaisesRegex(ValueError, "valid APY event ledger"):
            run_event_ledger_health_economics(working_result()["eventLedger"], {})

    def test_setup_default_custom_validation_persistence_and_single_restore(self):
        app = app_for().run()
        self.assertFalse(app.exception)
        self.assertIn("Progression risk among people with TBI", text(app))
        self.assertEqual(radio(app, "Progression-risk scenario").value, PROGRESSION_RISK_LABELS["risk_enriched"])
        self.assertIn("not a calibrated or independently estimated", text(app))
        self.assertNotIn(CONTRACT_VERSION, text(app))
        radio(app, "Progression-risk scenario").set_value(PROGRESSION_RISK_LABELS["custom"]).run()
        custom = next(item for item in app.number_input if item.label.startswith("Custom progression-hazard"))
        custom.set_value(0.0).run()
        self.assertIn("positive", text(app))
        radio(app, "Workflow").set_value("Run analysis").run()
        self.assertNotIn("Run explicit TBI analysis", [item.label for item in app.button])
        radio(app, "Workflow").set_value("Set up").run()
        value = 3.141592653589793
        next(item for item in app.number_input if item.label.startswith("Custom progression-hazard")).set_value(value).run()
        app.run()
        radio(app, "Workflow").set_value("Run analysis").run()
        self.assertEqual(app.session_state[STATE_KEY]["customProgressionRiskMultiplier"], value)
        self.assertEqual(progression_risk_for_state(app.session_state[STATE_KEY])["multiplier"], value)
        radio(app, "Workflow").set_value("Set up").run()
        self.assertEqual([item.label for item in app.button].count("Restore working defaults"), 1)
        button(app, "Restore working defaults").click().run()
        self.assertEqual(app.session_state[STATE_KEY], default_workflow_state())
        self.assertEqual(radio(app, "Progression-risk scenario").value, PROGRESSION_RISK_LABELS["risk_enriched"])

    def test_standard_selection_marks_completed_result_stale(self):
        app = app_for(result=working_result()).run()
        radio(app, "Progression-risk scenario").set_value(PROGRESSION_RISK_LABELS["standard"]).run()
        self.assertTrue(app.session_state[STALE_KEY])
        radio(app, "Workflow").set_value("Results").run()
        self.assertFalse(app.exception)
        self.assertIn("stale", text(app).lower())
        self.assertEqual(app.session_state[STATE_KEY]["progressionRiskMode"], "standard")
        # The stored run is labelled by its own selected multiplier, not relabelled.
        self.assertIn("8×", text(app))

    def test_one_dispatch_and_visible_reconciled_components_and_sensitivity(self):
        result = working_result()
        app = app_for(page="Run analysis").run()
        with patch("app.explicit_tbi_workflow.run_explicit_workflow_analysis", return_value=deepcopy(result)) as dispatch:
            button(app, "Run explicit TBI analysis").click().run()
            self.assertEqual(dispatch.call_count, 1)
            self.assertEqual(app.session_state[DISPATCH_COUNT_KEY], 1)
            button(app, "Open Results").click().run()
            self.assertFalse(app.exception)
            self.assertEqual(dispatch.call_count, 1)
        self.assertIn("Expected future active TB cases over 20 years among the baseline population", [item.label for item in app.metric])
        self.assertIn("No future infections or endogenous transmission", text(app))
        table = next(item.value for item in app.dataframe if "Baseline group" in item.value.columns)
        self.assertEqual(table["Baseline group"].tolist(), ["Recently infected within 5 years", "Remote-only infection", "Uninfected", "Baseline active TB", "Total"])
        cases = table["Expected comparator active-TB cases"].tolist()
        self.assertEqual(sum(Decimal(value.replace(",", "")) for value in cases[:3]), Decimal(cases[-1].replace(",", "")))
        self.assertEqual(cases[3], "Excluded from incident projection")
        sensitivity = next(item.value for item in app.dataframe if "Scenario" in item.value.columns)
        self.assertEqual(sensitivity["Multiplier"].tolist(), ["1×", "8×"])
        decoded = json.loads(technical_results_json(app.session_state[RESULT_KEY]))
        self.assertEqual(decoded[ADJUSTMENT_FIELD]["multiplier"], 8.0)
        self.assertEqual(decoded["comparatorDecomposition"][0]["expectedComparatorActiveTBCases"], result["comparatorDecomposition"][0]["expectedComparatorActiveTBCases"])
        self.assertEqual(radio(app, "Workflow").options, list(WORKFLOW_PAGES))
        self.assertNotIn("Economics", radio(app, "Workflow").options)

    def test_contamination_rejected_and_old_or_machinery_remains_inactive(self):
        app = app_for(page="Results")
        app.session_state["results_bundle"] = {"analysisPathway": "frozen_sa_health"}
        app.session_state["economics_results"] = {"ICER": 99}
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(radio(app, "Workflow").value, "Set up")
        self.assertIsNone(app.session_state[RESULT_KEY])
        self.assertNotIn("economics_results", app.session_state)
        _s, cfg, _cal, _reg = scenario()
        self.assertEqual(cfg["explicitRecentRemoteTBI"]["riskFactorProgressionPolicy"], "none")
        self.assertEqual(cfg["explicitRecentRemoteTBI"]["progressionCurve"]["riskFactorProgressionPolicy"], "none")
        for word in ("Scenario", "Economics", "Decision Analysis", "Dynamic Model", "Compare"):
            self.assertNotIn(word, radio(app, "Workflow").options)
        self.assertFalse((ROOT / "explicit_tbi_app.py").exists())

    def test_display_rounding_reconciles_without_changing_source(self):
        for values in ([0.335, 0.335, 0.0], [1.005, 0.0, 0.0], [1.234567, 8.901234, 0.0]):
            before = list(values)
            displayed = format_reconciled_components(values, sum(values))
            self.assertEqual(sum(Decimal(value) for value in displayed), Decimal(format_number(sum(values))))
            self.assertEqual(values, before)
        self.assertEqual(format_multiplier(8.0), "8×")
        self.assertEqual(format_multiplier(0.001), "1e-3×")
        self.assertEqual(format_multiplier(1e12), "1e12×")


if __name__ == "__main__":
    unittest.main()
