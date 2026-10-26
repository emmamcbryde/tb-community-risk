from __future__ import annotations

from datetime import datetime, timezone
from copy import deepcopy
import json
import math
from time import perf_counter, process_time
from typing import Any, Mapping

import numpy as np
import pandas as pd

from engine.apy.calibration_policy import resolve_calibration_for_config
from engine.apy.config import normalise_config
from engine.apy.explicit_recent_remote_progression import (
    REINFECTION_POLICY_NO_RESET_CLOCK,
    REINFECTION_POLICY_RESET_CLOCK,
    STATE_BASELINE_ACTIVE_TB,
    STATE_RECENT,
    STATE_REMOTE_ONLY,
    STATE_UNINFECTED,
)
from engine.apy.explicit_recent_remote_runner import (
    ACTIVE_TB_OBSERVATION_POLICY_VALIDATION_ONLY,
    EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
    SCREENING_COVERAGE_DEFINITION,
    SCREENING_COVERAGE_POLICY_IDENTIFIER,
    TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION,
    deterministic_explicit_recent_remote_intervention_expectations,
    deterministic_explicit_recent_remote_expectations,
    explicit_recent_remote_analysis_config,
    explicit_recent_remote_runner_config_hash,
    explicit_recent_remote_event_ledger_unavailable,
    validate_deterministic_intervention_accounting,
)
from engine.apy.population_tbi_progression_risk import (
    ADJUSTMENT_FIELD, build_progression_risk_adjustment,
    validate_progression_risk_adjustment,
)
from engine.apy.regimen import (
    apply_regimen_overrides,
    default_regimen_library,
    get_regimen_from_library,
)
from engine.apy.runner import run_replicates


DETERMINISTIC_LABEL = "Quick deterministic preview — single expected-value calculation"
STOCHASTIC_LABEL = "Stochastic analysis — multiple simulated populations"
WORKFLOW_PAGES = (
    "Set up",
    "Run analysis",
    "Results",
    "Evidence and technical information",
)
DETERMINISTIC_PREVIEW_INFECTION_QUADRATURE_POINTS = 16
DETERMINISTIC_PREVIEW_SCREEN_QUADRATURE_POINTS = 12
PROGRESSION_RISK_LABELS = {
    "standard": "Standard progression risk — 1×",
    "risk_enriched": "Risk-enriched population — 8×",
    "custom": "Custom multiplier",
}


def default_workflow_state() -> dict[str, Any]:
    return {
        "population": 1500,
        "ageDistributionLabel": "Inherited APY exact-age distribution",
        "baselineActiveTBProportion": 0.0,
        "baselineActiveTBSource": "Working default: no baseline active TB assigned",
        "baselineActiveTBReviewStatus": "working_default",
        "recentTBITarget": 0.02,
        "remoteOnlyTBITarget": 0.08,
        "progressionRiskMode": "risk_enriched",
        "customProgressionRiskMultiplier": 8.0,
        "recentTBISource": "Working demonstration value",
        "remoteOnlyTBISource": "Working demonstration value",
        "targetReferenceYear": 2026,
        "targetReviewStatus": "working_demo_not_approved_default",
        "analysisHorizonYears": 20.0,
        "screeningCoverage": 0.30,
        "screeningWindowYears": 3.0,
        "testType": "IGRA",
        "testSensitivity": 0.95,
        "testSpecificity": 0.98,
        "tstSensitivity": 0.80,
        "tstSpecificityBCG": 0.55,
        "tstSpecificityNoBCG": 0.97,
        "regimen": "3HP",
        "pStartTPT": 0.85,
        "regimenPComplete": 0.80,
        "regimenADRstop": 0.05,
        "regimenEffFull": 0.85,
        "analysisMethod": DETERMINISTIC_LABEL,
        "stochasticRepetitions": 10,
        "seed": 20261013,
        "reinfectionPolicy": REINFECTION_POLICY_RESET_CLOCK,
        "activeTBObservationEnabled": False,
        "activeTBObservedCases": 0.0,
        "activeTBDenominator": 1500.0,
        "activeTBObservationDurationYears": 2.0,
        "activeTBObservationClassification": "prospective incident active TB",
        "activeTBAscertainment": "complete ascertainment explicit assumption",
        "activeTBAscertainmentProbability": 1.0,
        "activeTBObservationSource": "Not supplied",
        "activeTBObservationReviewStatus": "not_reviewed",
        "activeTBObservationNotes": "",
    }


def proportion_to_percent(value: float) -> float:
    return 100.0 * float(value)


def percent_to_proportion(value: float) -> float:
    return float(value) / 100.0


def source_for_value(state: Mapping[str, Any], key: str, default_key: str | None = None) -> str:
    defaults = default_workflow_state()
    reference_key = default_key or key
    if abs(float(state.get(key, 0.0)) - float(defaults.get(reference_key, 0.0))) > 1e-12:
        return "User-defined"
    if key == "recentTBITarget":
        return str(state.get("recentTBISource") or defaults["recentTBISource"])
    if key == "remoteOnlyTBITarget":
        return str(state.get("remoteOnlyTBISource") or defaults["remoteOnlyTBISource"])
    if key == "baselineActiveTBProportion":
        return str(
            state.get("baselineActiveTBSource")
            or defaults["baselineActiveTBSource"]
        )
    return "Working demonstration value"


def build_explicit_workflow_config(state: Mapping[str, Any]) -> dict[str, Any]:
    n_reps = (
        int(state["stochasticRepetitions"])
        if state.get("analysisMethod") == STOCHASTIC_LABEL
        else 1
    )
    target_source = (
        f"Recent: {source_for_value(state, 'recentTBITarget')}; "
        f"remote-only: {source_for_value(state, 'remoteOnlyTBITarget')}"
    )
    cfg = explicit_recent_remote_analysis_config(
        N=int(state["population"]),
        nReps=n_reps,
        seed=int(state["seed"]),
        screeningStrategy="random",
        screenCoverage=float(state["screeningCoverage"]),
        screeningWindowYears=float(state["screeningWindowYears"]),
        screenWindow=float(state["screeningWindowYears"]),
        followUpHorizonYears=float(state["analysisHorizonYears"]),
        followHorizon=float(state["analysisHorizonYears"]),
        testType=str(state["testType"]),
        testSensitivity=float(state["testSensitivity"]),
        testSpecificity=float(state["testSpecificity"]),
        tstSensitivity=float(state["tstSensitivity"]),
        tstSpecificityBCG=float(state["tstSpecificityBCG"]),
        tstSpecificityNoBCG=float(state["tstSpecificityNoBCG"]),
        regimen=str(state["regimen"]),
        pStartTPT=float(state["pStartTPT"]),
        regimenPComplete=float(state["regimenPComplete"]),
        regimenADRstop=float(state["regimenADRstop"]),
        regimenEffFull=float(state["regimenEffFull"]),
        baselineActiveTBProportion=float(state["baselineActiveTBProportion"]),
        baselineActiveTBSource=str(state["baselineActiveTBSource"]),
        baselineActiveTBReviewStatus=str(state["baselineActiveTBReviewStatus"]),
        recent_tbi_target=float(state["recentTBITarget"]),
        remote_only_tbi_target=float(state["remoteOnlyTBITarget"]),
        target_source=target_source,
        target_reference_year=int(state["targetReferenceYear"]),
        review_status=str(state["targetReviewStatus"]),
        notes="Explicit recent/remote TBI development workflow.",
        reinfection_policy=str(state["reinfectionPolicy"]),
        progression_risk_adjustment=progression_risk_for_state(state),
    )
    cfg["analysisMethod"] = (
        "agent_based" if state.get("analysisMethod") == STOCHASTIC_LABEL else "expected_value"
    )
    cfg["screeningCoveragePolicyIdentifier"] = SCREENING_COVERAGE_POLICY_IDENTIFIER
    cfg["screeningCoverageDefinition"] = SCREENING_COVERAGE_DEFINITION
    return cfg


def progression_risk_for_state(state: Mapping[str, Any]) -> dict[str, Any]:
    return build_progression_risk_adjustment(
        state.get("progressionRiskMode"),
        multiplier=state.get("customProgressionRiskMultiplier") if state.get("progressionRiskMode") == "custom" else None,
        notes="Population risk enrichment scenario. Not calibrated to an active-TB target or independently estimated. Inherited risk-factor ORs remain inactive.",
    )


def build_regimen_for_state(state: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, Any]:
    library = default_regimen_library(config.get("partialShortCourseMode") or "threshold80")
    regimen = get_regimen_from_library(str(state["regimen"]), library)
    return apply_regimen_overrides(regimen, dict(config))


def validate_workflow_state(state: Mapping[str, Any], *, check_calibration: bool = True) -> list[str]:
    try:
        progression_risk_for_state(state)
    except ValueError as exc:
        return [str(exc)]
    checks = [
        ("total modelled population", float(state.get("population", 0.0)), 1.0, None),
        ("baseline active TB", float(state.get("baselineActiveTBProportion", -1.0)), 0.0, 1.0),
        ("recent TBI target", float(state.get("recentTBITarget", -1.0)), 0.0, 1.0),
        ("remote-only TBI target", float(state.get("remoteOnlyTBITarget", -1.0)), 0.0, 1.0),
        ("screening coverage", float(state.get("screeningCoverage", -1.0)), 0.0, 1.0),
        ("analysis horizon", float(state.get("analysisHorizonYears", 0.0)), 0.01, None),
        ("screening period", float(state.get("screeningWindowYears", 0.0)), 0.01, None),
        ("treatment initiation", float(state.get("pStartTPT", -1.0)), 0.0, 1.0),
        ("completion", float(state.get("regimenPComplete", -1.0)), 0.0, 1.0),
        ("adverse-event stop", float(state.get("regimenADRstop", -1.0)), 0.0, 1.0),
        ("full-course efficacy", float(state.get("regimenEffFull", -1.0)), 0.0, 1.0),
    ]
    for label, value, lower, upper in checks:
        if value < lower or (upper is not None and value > upper):
            return [f"{label} must be between {lower:g} and {upper:g}." if upper is not None else f"{label} must be at least {lower:g}."]
    occupied = (
        float(state["baselineActiveTBProportion"])
        + float(state["recentTBITarget"])
        + float(state["remoteOnlyTBITarget"])
    )
    if occupied > 1.0 + 1e-12:
        return [
            "Baseline active TB plus recent and remote-only TBI targets must not exceed 100%."
        ]
    if float(state["regimenPComplete"]) + float(state["regimenADRstop"]) > 1.0 + 1e-12:
        return ["Completion plus adverse-event discontinuation must not exceed 100%."]
    if state.get("analysisMethod") == STOCHASTIC_LABEL:
        if int(state.get("stochasticRepetitions", 0)) <= 0:
            return ["Stochastic repetitions must be at least one."]
    if bool(state.get("activeTBObservationEnabled")):
        if float(state.get("activeTBDenominator", 0.0)) <= 0.0:
            return ["Active-TB observation denominator must be positive."]
        if float(state.get("activeTBObservedCases", -1.0)) < 0.0:
            return ["Observed active-TB cases must be non-negative."]
        if float(state.get("activeTBObservationDurationYears", 0.0)) <= 0.0:
            return ["Active-TB observation duration must be positive."]
        q = float(state.get("activeTBAscertainmentProbability", 0.0))
        if q <= 0.0 or q > 1.0:
            return ["Ascertainment probability must be in (0, 1]."]
    if not check_calibration:
        return []
    try:
        cfg = build_explicit_workflow_config(state)
        resolve_calibration_for_config(cfg)
    except Exception as exc:
        return [f"Recent/remote calibration is infeasible: {exc}"]
    return []


def health_signature(state: Mapping[str, Any]) -> str:
    material_keys = [
        "population",
        "baselineActiveTBProportion",
        "recentTBITarget",
        "remoteOnlyTBITarget",
        "analysisHorizonYears",
        "screeningCoverage",
        "screeningWindowYears",
        "testType",
        "testSensitivity",
        "testSpecificity",
        "tstSensitivity",
        "tstSpecificityBCG",
        "tstSpecificityNoBCG",
        "regimen",
        "pStartTPT",
        "regimenPComplete",
        "regimenADRstop",
        "regimenEffFull",
        "analysisMethod",
        "stochasticRepetitions",
        "seed",
        "reinfectionPolicy",
    ]
    payload = {key: state.get(key) for key in material_keys}
    try:
        payload[ADJUSTMENT_FIELD] = progression_risk_for_state(state)
    except ValueError:
        payload[ADJUSTMENT_FIELD] = {"mode": state.get("progressionRiskMode"), "invalidCustomMultiplier": state.get("customProgressionRiskMultiplier")}
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def run_explicit_workflow_analysis(state: Mapping[str, Any]) -> dict[str, Any]:
    clock_start = perf_counter()
    cpu_start = process_time()
    errors = validate_workflow_state(state, check_calibration=False)
    if errors:
        raise ValueError(errors[0])
    validation_finished = perf_counter()
    config = build_explicit_workflow_config(state)
    calibration = resolve_calibration_for_config(config)
    calibration_finished = perf_counter()
    regimen = build_regimen_for_state(state, normalise_config(config))
    horizon = float(state["analysisHorizonYears"])
    population = int(state["population"])
    analysis_started = perf_counter()
    if state.get("analysisMethod") == STOCHASTIC_LABEL:
        runner_result = run_replicates(config=config, keep_example_cohort=False)
        analysis_finished = perf_counter()
        summary_rows = _summary_rows_from_stochastic(runner_result)
        raw_result = {
            "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
            "mode": "stochastic",
            "runnerResult": runner_result,
        }
    else:
        deterministic = deterministic_explicit_recent_remote_intervention_expectations(
            calibration,
            config=normalise_config(config),
            regimen=regimen,
            horizons_years=(horizon,),
            infection_quadrature_points=DETERMINISTIC_PREVIEW_INFECTION_QUADRATURE_POINTS,
            screen_quadrature_points=DETERMINISTIC_PREVIEW_SCREEN_QUADRATURE_POINTS,
        )
        analysis_finished = perf_counter()
        validation = deterministic.get("accountingValidation") or {}
        if not validation.get("isValid", False):
            raise ValueError(
                "Explicit deterministic accounting failed: "
                + "; ".join(validation.get("errors") or [])
            )
        summary_rows = _summary_rows_from_deterministic(
            deterministic,
            horizon_key=f"{horizon:g}",
            population=population,
        )
        raw_result = {
            "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
            "mode": "deterministic",
            "deterministicResult": deterministic,
            "accountingValidation": validation,
        }
    if raw_result.get("mode") == "stochastic":
        replicate_validations = [
            replicate.get("accountingValidation") or {}
            for replicate in raw_result["runnerResult"]["explicitEpidemiology"][
                "replicates"
            ]
        ]
        failed = [
            validation
            for validation in replicate_validations
            if not validation.get("isValid", False)
        ]
        if failed:
            raise ValueError(
                "Explicit stochastic accounting failed: "
                + "; ".join(failed[0].get("errors") or [])
            )
        raw_result["accountingValidation"] = {
            "isValid": True,
            "replicateCount": len(replicate_validations),
        }
    sensitivity_started = perf_counter()
    selected_projection = raw_result.get("deterministicResult")
    if selected_projection is None:
        selected_projection = deterministic_explicit_recent_remote_intervention_expectations(
            calibration, config=normalise_config(config), regimen=regimen,
            horizons_years=(horizon,),
            infection_quadrature_points=DETERMINISTIC_PREVIEW_INFECTION_QUADRATURE_POINTS,
            screen_quadrature_points=DETERMINISTIC_PREVIEW_SCREEN_QUADRATURE_POINTS,
        )
    adjustment = config["explicitRecentRemoteTBI"][ADJUSTMENT_FIELD]
    standard_adjustment = build_progression_risk_adjustment("standard")
    standard_config = deepcopy(config)
    standard_config["explicitRecentRemoteTBI"][ADJUSTMENT_FIELD] = standard_adjustment
    if adjustment["multiplier"] == 1.0:
        standard_projection = deepcopy(selected_projection)
        standard_projection[ADJUSTMENT_FIELD] = standard_adjustment
    else:
        standard_calibration = resolve_calibration_for_config(standard_config)
        standard_projection = deterministic_explicit_recent_remote_intervention_expectations(
            standard_calibration, config=normalise_config(standard_config), regimen=regimen,
            horizons_years=(horizon,),
            infection_quadrature_points=DETERMINISTIC_PREVIEW_INFECTION_QUADRATURE_POINTS,
            screen_quadrature_points=DETERMINISTIC_PREVIEW_SCREEN_QUADRATURE_POINTS,
        )
    decomposition = comparator_progression_decomposition(selected_projection, population, horizon)
    projection_validation = validate_projection_accounting(config, selected_projection, decomposition)
    standard_validation = validate_deterministic_intervention_accounting(standard_projection)
    if not projection_validation["isValid"] or not standard_validation["isValid"]:
        raise ValueError("Progression projection accounting failed: " + "; ".join(projection_validation["errors"] + standard_validation["errors"]))
    sensitivity = {
        "method": "analytical_progression_hazard_sensitivity_v1",
        "standardAdjustment": standard_adjustment,
        "selectedAdjustment": adjustment,
        "standardProjection": standard_projection,
        "standardConfiguration": standard_config,
        "standardConfigHash": explicit_recent_remote_runner_config_hash(standard_config["explicitRecentRemoteTBI"]),
        "rows": [
            {
                "Scenario": label, "Mode": contract["mode"], "Multiplier": contract["multiplier"],
                "Comparator active-TB cases": projection["byHorizon"][f"{horizon:g}"]["comparatorActiveTBProportion"] * population,
                "Intervention active-TB cases": projection["byHorizon"][f"{horizon:g}"]["interventionActiveTBProportion"] * population,
                "Direct active-TB cases averted": projection["byHorizon"][f"{horizon:g}"]["directActiveTBCasesAvertedProportion"] * population,
            }
            for label, contract, projection in (
                ("Standard progression risk", standard_adjustment, standard_projection),
                ("Selected progression risk", adjustment, selected_projection),
            )
        ],
    }
    sensitivity_finished = perf_counter()
    natural_history = deterministic_explicit_recent_remote_expectations(
        calibration,
        horizons_years=(horizon,),
        population_size=float(population),
    )
    age_rows = age_specific_prevalence_table(
        calibration,
        horizon_years=horizon,
        population_size=float(population),
        natural_history=natural_history,
    )
    active_tb_validation = active_tb_validation_for_state(state, natural_history)
    ledger_marker = explicit_recent_remote_event_ledger_unavailable(config["explicitRecentRemoteTBI"])
    config_hash = explicit_recent_remote_runner_config_hash(config["explicitRecentRemoteTBI"])
    finished = datetime.now(timezone.utc)
    return {
        **raw_result,
        "createdAt": finished.isoformat(),
        "elapsedSeconds": perf_counter() - clock_start,
        "timings": {
            "processCpuSeconds": process_time() - cpu_start,
            "inputValidationSeconds": validation_finished - clock_start,
            "calibrationSeconds": calibration_finished - validation_finished,
            "analysisPreparationSeconds": analysis_started - calibration_finished,
            "deterministicQuadratureSeconds": analysis_finished - analysis_started if state.get("analysisMethod") != STOCHASTIC_LABEL else 0.0,
            "stochasticSimulationSeconds": analysis_finished - analysis_started if state.get("analysisMethod") == STOCHASTIC_LABEL else 0.0,
            "analyticalSensitivitySeconds": sensitivity_finished - sensitivity_started,
            "resultConstructionSeconds": perf_counter() - sensitivity_finished,
        },
        "healthSignature": health_signature(state),
        "config": config,
        "configHash": config_hash,
        "eventLedger": ledger_marker,
        ADJUSTMENT_FIELD: adjustment,
        "analyticalProjection": selected_projection,
        "comparatorDecomposition": decomposition,
        "projectionAccountingValidation": projection_validation,
        "progressionRiskSensitivity": sensitivity,
        "calibration": calibration,
        "summaryRows": summary_rows,
        "ageRows": age_rows,
        "stateProportions": natural_history["stateProportions"],
        "activeTBValidation": active_tb_validation,
        "limitations": [
            "Direct active TB cases averted among the modelled population only.",
            "Active-TB observations are validation-only and do not recalibrate the progression curve.",
            "No endogenous transmission, cost, disability-weight or quality-adjusted outputs are produced.",
            "The population progression-hazard multiplier is a scenario assumption, not calibrated or independently estimated.",
            "No new infections occur during follow-up.",
        ],
    }


def comparator_progression_decomposition(projection: Mapping[str, Any], population: int, horizon: float) -> list[dict[str, Any]]:
    states = projection["stateProportions"]
    outcome = projection["byHorizon"][f"{horizon:g}"]
    contributions = {
        STATE_RECENT: outcome["recentComparatorActiveTBProportion"] * population,
        STATE_REMOTE_ONLY: outcome["remoteOnlyComparatorActiveTBProportion"] * population,
        STATE_UNINFECTED: 0.0,
        STATE_BASELINE_ACTIVE_TB: None,
    }
    rows = []
    for state, label in (
        (STATE_RECENT, "Recently infected within 5 years"),
        (STATE_REMOTE_ONLY, "Remote-only infection"),
        (STATE_UNINFECTED, "Uninfected"),
        (STATE_BASELINE_ACTIVE_TB, "Baseline active TB"),
    ):
        people = population * states[state]
        cases = contributions[state]
        rows.append({
            "state": state, "baselineGroup": label, "numberOfPeople": people,
            "averageProgressionProbability": None if cases is None else (cases / people if people > 0.0 else 0.0),
            "expectedComparatorActiveTBCases": cases,
        })
    rows.append({"state": "total", "baselineGroup": "Total", "numberOfPeople": population,
                 "averageProgressionProbability": None,
                 "expectedComparatorActiveTBCases": population * outcome["comparatorActiveTBProportion"]})
    return rows


def validate_projection_accounting(config: Mapping[str, Any], projection: Mapping[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Hard checks on unrounded projection values, before any display formatting."""
    checked = validate_deterministic_intervention_accounting(projection)
    errors = list(checked["errors"])
    population = float(config["N"])
    states = {row["state"]: row for row in rows}
    required = {STATE_RECENT, STATE_REMOTE_ONLY, STATE_UNINFECTED, STATE_BASELINE_ACTIVE_TB, "total"}
    if set(states) != required or len(rows) != 5:
        return {"isValid": False, "errors": ["Comparator decomposition must contain each baseline group and a total exactly once."]}
    def close(left, right, label, tolerance=1e-8):
        if not math.isfinite(float(left)) or not math.isfinite(float(right)) or not math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance * max(1.0, population)):
            errors.append(label)
    close(sum(states[key]["numberOfPeople"] for key in required if key != "total"), population, "Baseline groups must sum to the population.")
    close(states["total"]["numberOfPeople"], population, "Total population mismatch.")
    runner = config["explicitRecentRemoteTBI"]
    if projection.get(ADJUSTMENT_FIELD) != runner.get(ADJUSTMENT_FIELD):
        errors.append("Projection progression-risk contract does not match configuration.")
    close(states[STATE_RECENT]["numberOfPeople"], population * runner["recentTBITargetProportion"], "Recent TBI target mismatch.")
    close(states[STATE_REMOTE_ONLY]["numberOfPeople"], population * runner["remoteOnlyTBITargetProportion"], "Remote-only TBI target mismatch.")
    close(states[STATE_BASELINE_ACTIVE_TB]["numberOfPeople"], population * runner["baselineActiveTBProportion"], "Baseline active TB mismatch.")
    close(states[STATE_UNINFECTED]["expectedComparatorActiveTBCases"], 0.0, "Uninfected progression must be zero.")
    if states[STATE_BASELINE_ACTIVE_TB]["expectedComparatorActiveTBCases"] is not None or states[STATE_BASELINE_ACTIVE_TB]["averageProgressionProbability"] is not None:
        errors.append("Baseline active TB must be excluded from incident projection.")
    for key in (STATE_RECENT, STATE_REMOTE_ONLY, STATE_UNINFECTED):
        row = states[key]
        probability = row["averageProgressionProbability"]
        # At saturated hazards, weighted quadrature/count division can exceed
        # one by a few floating-point ulps. Accept numerical error only; retain
        # the raw value and never clip the multiplier or a multiplied risk.
        if not math.isfinite(probability) or not 0.0 <= probability <= 1.0 + 1e-12:
            errors.append("Invalid group progression probability.")
        close(row["numberOfPeople"] * probability, row["expectedComparatorActiveTBCases"], "Group probability/count mismatch.")
    horizon = projection["byHorizon"][f"{float(config['followUpHorizonYears']):g}"]
    group_cases = states[STATE_RECENT]["expectedComparatorActiveTBCases"] + states[STATE_REMOTE_ONLY]["expectedComparatorActiveTBCases"]
    close(group_cases, states["total"]["expectedComparatorActiveTBCases"], "Recent plus remote-only cases must sum to total.")
    close(group_cases, population * horizon["comparatorActiveTBProportion"], "Decomposition must equal comparator incident cases.")
    close(horizon["interventionActiveTBProportion"] + horizon["directActiveTBCasesAvertedProportion"], horizon["comparatorActiveTBProportion"], "Intervention plus cases averted must equal comparator.")
    return {"isValid": not errors, "errors": errors, "relativeAndPopulationScaledAbsoluteTolerance": 1e-8, "groupProbabilityUpperBoundTolerance": 1e-12}


def active_tb_validation_for_state(
    state: Mapping[str, Any],
    natural_history: Mapping[str, Any] | None = None,
) -> dict[str, Any] | None:
    if not bool(state.get("activeTBObservationEnabled")):
        return None
    denominator = float(state["activeTBDenominator"])
    observed = float(state["activeTBObservedCases"])
    duration = float(state["activeTBObservationDurationYears"])
    if natural_history is None:
        cfg = build_explicit_workflow_config(state)
        calibration = resolve_calibration_for_config(cfg)
        natural_history = deterministic_explicit_recent_remote_expectations(
            calibration,
            horizons_years=(duration,),
            population_size=denominator,
        )
    horizon_row = natural_history["incidentActiveTBByHorizon"].get(f"{duration:g}")
    if horizon_row is None:
        cfg = build_explicit_workflow_config(state)
        calibration = resolve_calibration_for_config(cfg)
        natural_history = deterministic_explicit_recent_remote_expectations(
            calibration,
            horizons_years=(duration,),
            population_size=denominator,
        )
        horizon_row = natural_history["incidentActiveTBByHorizon"][f"{duration:g}"]
    q = float(state.get("activeTBAscertainmentProbability", 1.0))
    expected = float(horizon_row["incidentActiveTBProportion"]) * denominator * q
    observed_rate = observed / denominator
    expected_rate = expected / denominator
    ratio = None if expected <= 0.0 else observed / expected
    return {
        "status": "validation_only",
        "observationPolicy": ACTIVE_TB_OBSERVATION_POLICY_VALIDATION_ONLY,
        "observedCases": observed,
        "sourcePopulationDenominator": denominator,
        "observationDurationYears": duration,
        "classification": str(state["activeTBObservationClassification"]),
        "ascertainment": str(state["activeTBAscertainment"]),
        "ascertainmentProbability": q,
        "source": str(state["activeTBObservationSource"]),
        "reviewStatus": str(state["activeTBObservationReviewStatus"]),
        "expectedCases": expected,
        "observedRate": observed_rate,
        "expectedRate": expected_rate,
        "absoluteDifference": observed - expected,
        "ratioObservedToExpected": ratio,
        "message": (
            "These observations compare observed and model-expected active TB. "
            "They do not automatically recalibrate the progression curve."
        ),
    }


def age_specific_prevalence_table(
    calibration: Mapping[str, Any],
    *,
    horizon_years: float,
    population_size: float,
    natural_history: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    if natural_history is None:
        natural_history = deterministic_explicit_recent_remote_expectations(
            calibration,
            horizons_years=(float(horizon_years),),
            population_size=float(population_size),
        )
    age_rows = natural_history["ageSpecificStateDistributions"]
    incident_rows = [
        row
        for row in natural_history["ageSpecificIncidentActiveTBByHorizon"]
        if abs(float(row["horizonYears"]) - float(horizon_years)) < 1e-9
    ]
    incident_by_age = {float(row["ageYears"]): row for row in incident_rows}
    groups: dict[str, dict[str, float]] = {}
    for row in age_rows:
        age = float(row["ageYears"])
        label = age_group_label(age)
        weight = float(row["ageWeight"])
        target = groups.setdefault(
            label,
            {
                "populationShare": 0.0,
                "recent": 0.0,
                "remote": 0.0,
                "uninfected": 0.0,
                "risk": 0.0,
            },
        )
        target["populationShare"] += weight
        target["recent"] += weight * float(row["recentProbabilityTotalPopulation"])
        target["remote"] += weight * float(row["remoteOnlyProbabilityTotalPopulation"])
        target["uninfected"] += weight * float(row["uninfectedProbabilityTotalPopulation"])
        incident = incident_by_age.get(age)
        if incident:
            target["risk"] += float(incident["incidentActiveTBProbabilityTotalPopulation"])
    output = []
    for label in sorted(groups, key=_age_group_sort_key):
        row = groups[label]
        share = row["populationShare"]
        output.append(
            {
                "Age group": label,
                "Population share": share,
                "Recent infection prevalence": 0.0 if share <= 0 else row["recent"] / share,
                "Remote-only infection prevalence": 0.0 if share <= 0 else row["remote"] / share,
                "Uninfected prevalence": 0.0 if share <= 0 else row["uninfected"] / share,
                "Expected active-TB risk": 0.0 if share <= 0 else row["risk"] / share,
            }
        )
    return output


def age_group_label(age: float) -> str:
    if age < 5:
        return "0-4"
    if age < 15:
        return "5-14"
    if age < 25:
        return "15-24"
    if age < 35:
        return "25-34"
    if age < 45:
        return "35-44"
    if age < 55:
        return "45-54"
    if age < 65:
        return "55-64"
    if age < 75:
        return "65-74"
    if age < 85:
        return "75-84"
    return "85+"


def results_to_csv(result: Mapping[str, Any]) -> str:
    return pd.DataFrame(result.get("summaryRows") or []).to_csv(index=False)


def validation_to_csv(validation: Mapping[str, Any] | None) -> str:
    if not validation:
        return "status,message\nnot_supplied,No active-TB observation supplied\n"
    return pd.DataFrame([validation]).to_csv(index=False)


def provenance_text(result: Mapping[str, Any]) -> str:
    cfg = result["config"]["explicitRecentRemoteTBI"]
    return "\n".join(
        [
            "Explicit recent/remote TBI development workflow",
            f"Pathway: {cfg['pathwayIdentifier']}",
            f"Runner contract: {cfg['contractVersion']}",
            f"Progression curve: {cfg['progressionCurveIdentifier']}",
            f"TPT effect policy: {cfg['tptEffectPolicyIdentifier']}",
            f"Risk-factor progression policy: {cfg['riskFactorProgressionPolicy']}",
            f"Population progression-risk adjustment: {json.dumps(cfg[ADJUSTMENT_FIELD], sort_keys=True)}",
            "Population enrichment is an assumption, not a calibrated or independently estimated relative risk.",
            "No future infections or endogenous transmission are modelled.",
            f"Active-TB observation policy: {cfg['activeTBObservationPolicy']}",
            f"Screening coverage policy: {cfg['screeningCoveragePolicyIdentifier']}",
            f"Screening coverage definition: {cfg['screeningCoverageDefinition']}",
            f"Analysis method: {result['mode']}",
            "Limitation: direct active TB cases averted among the modelled population only.",
        ]
    )


def discard_incompatible_result(result: Any) -> Any:
    if not isinstance(result, Mapping):
        return None
    if result.get("analysisPathway") != EXPLICIT_RECENT_REMOTE_PATHWAY_ID:
        return None
    if result.get("mode") not in ("deterministic", "stochastic"):
        return None
    config = result.get("config") or {}
    if not isinstance(config, Mapping):
        return None
    if config.get("analysisPathway") != EXPLICIT_RECENT_REMOTE_PATHWAY_ID:
        return None
    pathway = config.get("explicitRecentRemoteTBI") or {}
    if not isinstance(pathway, Mapping):
        return None
    if pathway.get("pathwayIdentifier") != EXPLICIT_RECENT_REMOTE_PATHWAY_ID:
        return None
    try:
        adjustment = validate_progression_risk_adjustment(pathway.get(ADJUSTMENT_FIELD))
        if result.get(ADJUSTMENT_FIELD) != adjustment:
            return None
    except ValueError:
        return None
    if not all(key in result for key in ("analyticalProjection", "comparatorDecomposition", "progressionRiskSensitivity", "projectionAccountingValidation")):
        return None
    if not all(key in result for key in ("healthSignature", "summaryRows", "ageRows", "calibration", "configHash")):
        return None
    if not isinstance(result.get("summaryRows"), list) or not isinstance(result.get("ageRows"), list):
        return None
    if not isinstance(result.get("calibration"), Mapping):
        return None
    payload_key = "deterministicResult" if result["mode"] == "deterministic" else "runnerResult"
    if not isinstance(result.get(payload_key), Mapping):
        return None
    if any(key in result for key in ("economics", "economics_results", "economicsResults", "dynamicResults", "decisionAnalysis")):
        return None
    ledger = result.get("eventLedger")
    if not isinstance(ledger, Mapping) or ledger.get("available") is not False:
        return None
    return result


def technical_results_json(result: Mapping[str, Any]) -> str:
    """Full-precision machine-readable provenance, independent of display formatting."""
    def encode(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"Unsupported provenance value: {type(value).__name__}")
    return json.dumps(result, indent=2, sort_keys=True, default=encode)


def _summary_rows_from_deterministic(
    deterministic: Mapping[str, Any],
    *,
    horizon_key: str,
    population: int,
) -> list[dict[str, Any]]:
    horizon = deterministic["byHorizon"][horizon_key]
    states = deterministic["stateProportions"]
    rows = [
        ("Population", population),
        ("Screened", population * horizon["screenedProportion"]),
        ("Recently infected within the past 5 years", population * states[STATE_RECENT]),
        (
            "Remote infection only — infected 5 or more years ago and not reinfected within the past 5 years",
            population * states[STATE_REMOTE_ONLY],
        ),
        ("Baseline active TB", population * states[STATE_BASELINE_ACTIVE_TB]),
        ("TPT starts", population * horizon["tptStartsProportion"]),
        ("TPT completions", population * horizon["tptCompletionsProportion"]),
        ("Comparator active TB", population * horizon["comparatorActiveTBProportion"]),
        ("Intervention active TB", population * horizon["interventionActiveTBProportion"]),
        (
            "Direct active TB cases averted among the modelled population",
            population * horizon["directActiveTBCasesAvertedProportion"],
        ),
    ]
    return [
        {
            "Metric": label,
            "Expected value": float(value),
            "Low 95% simulation interval": "N/A",
            "High 95% simulation interval": "N/A",
        }
        for label, value in rows
    ]


def _summary_rows_from_stochastic(result: Mapping[str, Any]) -> list[dict[str, Any]]:
    summary = result["explicitEpidemiology"]["metricSummary"]
    lookup = {
        (row["arm"], row["metric"]): row
        for row in summary
    }
    metric_map = [
        ("Population", "comparator", "population"),
        ("Screened", "intervention", "screened"),
        ("Recently infected within the past 5 years", "comparator", "recentTBI"),
        (
            "Remote infection only — infected 5 or more years ago and not reinfected within the past 5 years",
            "comparator",
            "remoteOnlyTBI",
        ),
        ("Baseline active TB", "comparator", "baselineActiveTB"),
        ("TPT starts", "intervention", "tptStarts"),
        ("TPT completions", "intervention", "tptCompletions"),
        ("Comparator active TB", "comparator", "totalActiveTBOverFollowUp"),
        ("Intervention active TB", "intervention", "totalActiveTBOverFollowUp"),
        (
            "Direct active TB cases averted among the modelled population",
            "intervention",
            "directActiveTBCasesAvertedAmongModelledPopulation",
        ),
    ]
    rows = []
    for label, arm, metric in metric_map:
        row = lookup[(arm, metric)]
        rows.append(
            {
                "Metric": label,
                "Mean": float(row["mean"]),
                "Low 95% simulation interval": float(row["low95"]),
                "High 95% simulation interval": float(row["high95"]),
            }
        )
    return rows


def _age_group_sort_key(label: str) -> int:
    order = {
        "0-4": 0,
        "5-14": 1,
        "15-24": 2,
        "25-34": 3,
        "35-44": 4,
        "45-54": 5,
        "55-64": 6,
        "65-74": 7,
        "75-84": 8,
        "85+": 9,
    }
    return order.get(label, 999)
