from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import json
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
    TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION,
    deterministic_explicit_recent_remote_intervention_expectations,
    deterministic_explicit_recent_remote_expectations,
    explicit_recent_remote_analysis_config,
    explicit_recent_remote_runner_config_hash,
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


def default_workflow_state() -> dict[str, Any]:
    return {
        "population": 1500,
        "ageDistributionLabel": "Inherited APY exact-age distribution",
        "baselineActiveTBProportion": 0.0,
        "baselineActiveTBSource": "Working default: no baseline active TB assigned",
        "baselineActiveTBReviewStatus": "working_default",
        "recentTBITarget": 0.02,
        "remoteOnlyTBITarget": 0.08,
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
    )
    cfg["analysisMethod"] = (
        "agent_based" if state.get("analysisMethod") == STOCHASTIC_LABEL else "expected_value"
    )
    return cfg


def build_regimen_for_state(state: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, Any]:
    library = default_regimen_library(config.get("partialShortCourseMode") or "threshold80")
    regimen = get_regimen_from_library(str(state["regimen"]), library)
    return apply_regimen_overrides(regimen, dict(config))


def validate_workflow_state(state: Mapping[str, Any]) -> list[str]:
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
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def run_explicit_workflow_analysis(state: Mapping[str, Any]) -> dict[str, Any]:
    errors = validate_workflow_state(state)
    if errors:
        raise ValueError(errors[0])
    started = datetime.now(timezone.utc)
    config = build_explicit_workflow_config(state)
    calibration = resolve_calibration_for_config(config)
    regimen = build_regimen_for_state(state, normalise_config(config))
    horizon = float(state["analysisHorizonYears"])
    population = int(state["population"])
    if state.get("analysisMethod") == STOCHASTIC_LABEL:
        runner_result = run_replicates(config=config, keep_example_cohort=False)
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
            infection_quadrature_points=32,
            screen_quadrature_points=24,
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
        }
    natural_history = deterministic_explicit_recent_remote_expectations(
        calibration,
        horizons_years=(horizon,),
        population_size=float(population),
    )
    age_rows = age_specific_prevalence_table(
        calibration,
        horizon_years=horizon,
        population_size=float(population),
    )
    finished = datetime.now(timezone.utc)
    return {
        **raw_result,
        "createdAt": finished.isoformat(),
        "elapsedSeconds": (finished - started).total_seconds(),
        "healthSignature": health_signature(state),
        "config": config,
        "configHash": explicit_recent_remote_runner_config_hash(
            config["explicitRecentRemoteTBI"]
        ),
        "calibration": calibration,
        "summaryRows": summary_rows,
        "ageRows": age_rows,
        "stateProportions": natural_history["stateProportions"],
        "activeTBValidation": active_tb_validation_for_state(state, natural_history),
        "limitations": [
            "Direct active TB cases averted among the modelled population only.",
            "Active-TB observations are validation-only and do not recalibrate the progression curve.",
            "No endogenous transmission, cost, disability-weight or quality-adjusted outputs are produced.",
        ],
    }


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
) -> list[dict[str, Any]]:
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
            f"Active-TB observation policy: {cfg['activeTBObservationPolicy']}",
            "Limitation: direct active TB cases averted among the modelled population only.",
        ]
    )


def discard_incompatible_result(result: Any) -> Any:
    if not isinstance(result, Mapping):
        return None
    if result.get("analysisPathway") != EXPLICIT_RECENT_REMOTE_PATHWAY_ID:
        return None
    return result


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

