from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from typing import Any, Mapping

import numpy as np

from engine.apy.data import load_parameters_from_config
from engine.apy.explicit_recent_remote_progression import (
    COMPETING_MORTALITY_NOT_MODELLED,
    COMPETING_MORTALITY_CONSTANT_HAZARD,
    CURVE_CENTRAL_TIME_SINCE_INFECTION,
    POST_FINAL_ANCHOR_FROM_LAST_SEGMENT,
    REINFECTION_POLICY_NO_RESET_CLOCK,
    REINFECTION_POLICY_RESET_CLOCK,
    RISK_POLICY_NONE,
    STATE_BASELINE_ACTIVE_TB,
    build_central_time_since_infection_progression_curve,
    build_progression_curve_by_identifier,
    competing_risk_progression_probability_for_curve,
    progression_curve_hash,
    time_since_curve_cumulative_hazard,
    time_since_curve_conditional_future_progression_probability,
    validate_progression_curve_contract,
)
from engine.apy.explicit_recent_remote_tbi import (
    AGE_SUPPORT_PROVENANCE,
    AGE85_PLUS_MAX_DEFAULT,
    CONFIG_CONTRACT_VERSION as TBI_CONFIG_CONTRACT_VERSION,
    RECENT_HAZARD_SHAPE,
    RECENT_WINDOW_YEARS,
    REMOTE_HAZARD_SHAPE,
    REMOTE_HISTORY_CAP_YEARS,
    STATE_RECENT,
    STATE_REMOTE_ONLY,
    STATE_UNINFECTED,
    build_explicit_recent_remote_config,
    calibrate_recent_remote_hazards,
    draw_recent_remote_states_for_ages,
    explicit_recent_remote_config_from_dict,
    explicit_recent_remote_config_hash,
    exposure_durations_for_ages,
    mutually_exclusive_state_probabilities,
    sample_time_since_most_recent_event,
)


EXPLICIT_RECENT_REMOTE_PATHWAY_ID = "explicit_recent_remote_tbi_v2"
EXPLICIT_RECENT_REMOTE_RUNNER_CONTRACT_VERSION = (
    "explicit_recent_remote_tbi_runner_config_v1"
)
EXPLICIT_RECENT_REMOTE_ANALYSIS_METHOD = (
    "isolated_agent_based_explicit_recent_remote_tbi_v2"
)
ACTIVE_TB_OBSERVATION_POLICY_VALIDATION_ONLY = "active_tb_validation_only_v1"
RUNNER_INTEGRATION_STATUS = "milestone_3b_isolated_runner_integration"
PROGRESSION_TIME_SAMPLING_METHOD = (
    "inverse_survivor_conditioned_time_since_infection_cumulative_hazard_v1"
)
BASELINE_ACTIVE_TB_ASSIGNMENT_POLICY = (
    "baseline_active_tb_separate_before_tbi_assignment_v1"
)
EXPLICIT_EVENT_LEDGER_OUT_OF_SCOPE_SOURCE = (
    "explicit_recent_remote_tbi_v2_event_ledger_out_of_scope_milestone_3b"
)
DETERMINISTIC_INFECTION_TIME_INTEGRATION_METHOD = (
    "gauss_legendre_quadrature_over_conditional_infection_time_v1"
)
STOCHASTIC_NATURAL_HISTORY_METHOD = (
    "isolated_explicit_recent_remote_natural_history_stochastic_v1"
)


def build_explicit_recent_remote_runner_config(
    *,
    recent_tbi_target: float,
    remote_only_tbi_target: float,
    target_source: str,
    target_reference_year: int | None = None,
    review_status: str = "approved_for_runner_integration_review",
    notes: str = "",
    recent_window_years: float = RECENT_WINDOW_YEARS,
    remote_history_cap_years: float = REMOTE_HISTORY_CAP_YEARS,
    age85_plus_max: int = AGE85_PLUS_MAX_DEFAULT,
    progression_curve_identifier: str = CURVE_CENTRAL_TIME_SINCE_INFECTION,
    reinfection_policy: str = REINFECTION_POLICY_RESET_CLOCK,
    risk_factor_progression_policy: str = RISK_POLICY_NONE,
    mortality_policy: str = COMPETING_MORTALITY_NOT_MODELLED,
    active_tb_observation_policy: str = ACTIVE_TB_OBSERVATION_POLICY_VALIDATION_ONLY,
    analysis_method: str = EXPLICIT_RECENT_REMOTE_ANALYSIS_METHOD,
    stochastic_seed: int | None = None,
    repetitions: int | None = None,
    baseline_active_tb_proportion: float = 0.0,
    baseline_active_tb_source: str = "not configured; baseline active TB assigned as zero",
    baseline_active_tb_review_status: str = "not_configured",
) -> dict[str, Any]:
    tbi_config = build_explicit_recent_remote_config(
        enabled=True,
        recent_tbi_target=recent_tbi_target,
        remote_only_tbi_target=remote_only_tbi_target,
        recent_window_years=recent_window_years,
        remote_history_cap_years=remote_history_cap_years,
        age85_plus_max=age85_plus_max,
        age_support_provenance=AGE_SUPPORT_PROVENANCE,
        recent_hazard_shape=RECENT_HAZARD_SHAPE,
        remote_hazard_shape=REMOTE_HAZARD_SHAPE,
        target_source=target_source,
        target_reference_year=target_reference_year,
        review_status=review_status,
        notes=notes,
    )
    curve = _build_progression_curve(
        progression_curve_identifier,
        reinfection_policy=reinfection_policy,
        mortality_policy=mortality_policy,
    )
    payload = {
        "contractVersion": EXPLICIT_RECENT_REMOTE_RUNNER_CONTRACT_VERSION,
        "pathwayIdentifier": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
        "enabled": True,
        "integrationStatus": RUNNER_INTEGRATION_STATUS,
        "analysisMethod": analysis_method,
        "explicitRecentRemoteTBIConfig": tbi_config.as_dict(),
        "explicitRecentRemoteTBIConfigHash": explicit_recent_remote_config_hash(
            tbi_config
        ),
        "targetDefinitions": tbi_config.as_dict()["targetDefinitions"],
        "recentTBITargetProportion": tbi_config.recent_tbi_target,
        "remoteOnlyTBITargetProportion": tbi_config.remote_only_tbi_target,
        "recentWindowYears": tbi_config.recent_window_years,
        "remoteHistoryCapYears": tbi_config.remote_history_cap_years,
        "age85PlusMax": tbi_config.age85_plus_max,
        "ageSupportProvenance": tbi_config.age_support_provenance,
        "recentHazardShape": tbi_config.recent_hazard_shape,
        "remoteHazardShape": tbi_config.remote_hazard_shape,
        "progressionCurveIdentifier": curve["curveIdentifier"],
        "progressionCurveHash": progression_curve_hash(curve),
        "progressionCurve": curve,
        "progressionAnchors": {
            "timeYears": list(curve["timeAnchors"]),
            "cumulativeRisk": list(curve["cumulativeRiskAnchors"]),
            "cumulativeHazard": list(curve["cumulativeHazardAnchors"]),
            "segmentHazards": list(curve["segmentHazards"]),
        },
        "postFinalAnchorRule": curve["postFinalAnchorExtrapolationMethod"],
        "postFinalAnchorHazard": curve["postFinalAnchorHazard"],
        "postFinalAnchorIsExtrapolation": (
            curve["postFinalAnchorExtrapolationMethod"]
            == POST_FINAL_ANCHOR_FROM_LAST_SEGMENT
        ),
        "progressionTimeSamplingMethod": PROGRESSION_TIME_SAMPLING_METHOD,
        "reinfectionPolicy": curve["reinfectionPolicy"],
        "riskFactorProgressionPolicy": risk_factor_progression_policy,
        "mortalityPolicy": mortality_policy,
        "mortalitySource": "none supplied; competing mortality not modelled",
        "baselineActiveTBAssignmentPolicy": BASELINE_ACTIVE_TB_ASSIGNMENT_POLICY,
        "baselineActiveTBProportion": float(baseline_active_tb_proportion),
        "baselineActiveTBSource": str(baseline_active_tb_source),
        "baselineActiveTBReviewStatus": str(baseline_active_tb_review_status),
        "activeTBObservationPolicy": active_tb_observation_policy,
        "stochasticSeed": None if stochastic_seed is None else int(stochastic_seed),
        "repetitions": None if repetitions is None else int(repetitions),
        "scientificInputs": {
            "tbiTargets": {
                "source": tbi_config.target_source,
                "referenceYear": tbi_config.target_reference_year,
                "reviewStatus": tbi_config.review_status,
                "notes": tbi_config.notes,
            },
            "progressionCurve": {
                "sources": deepcopy(curve["sources"]),
                "reviewStatus": curve["reviewStatus"],
                "notes": curve["notes"],
                "transformation": curve["transformation"],
            },
            "riskFactors": {
                "policy": risk_factor_progression_policy,
                "source": "Milestone 3B approved working decision",
                "reviewStatus": "use_none_initially",
                "notes": "Inherited disease ORs are not applied by the central explicit pathway.",
            },
            "mortality": {
                "policy": mortality_policy,
                "source": "No external mortality table supplied in Milestone 3B.",
                "reviewStatus": "not_modelled",
            },
            "activeTBObservations": {
                "policy": active_tb_observation_policy,
                "reviewStatus": "validation_only_by_default",
            },
            "baselineActiveTB": {
                "assignmentPolicy": BASELINE_ACTIVE_TB_ASSIGNMENT_POLICY,
                "proportion": float(baseline_active_tb_proportion),
                "source": str(baseline_active_tb_source),
                "reviewStatus": str(baseline_active_tb_review_status),
                "notes": (
                    "Baseline/prevalent active TB is sequenced before TBI states "
                    "and is excluded from the latent-TB treatment cascade."
                ),
            },
        },
    }
    return validate_explicit_recent_remote_runner_config(payload)


def explicit_recent_remote_analysis_config(
    *,
    recent_tbi_target: float,
    remote_only_tbi_target: float,
    target_source: str,
    target_reference_year: int | None = None,
    review_status: str = "approved_for_runner_integration_review",
    notes: str = "",
    progression_curve_identifier: str = CURVE_CENTRAL_TIME_SINCE_INFECTION,
    reinfection_policy: str = REINFECTION_POLICY_RESET_CLOCK,
    mortality_policy: str = COMPETING_MORTALITY_NOT_MODELLED,
    **overrides: Any,
) -> dict[str, Any]:
    runner_config = build_explicit_recent_remote_runner_config(
        recent_tbi_target=recent_tbi_target,
        remote_only_tbi_target=remote_only_tbi_target,
        target_source=target_source,
        target_reference_year=target_reference_year,
        review_status=review_status,
        notes=notes,
        progression_curve_identifier=progression_curve_identifier,
        reinfection_policy=reinfection_policy,
        mortality_policy=mortality_policy,
        stochastic_seed=overrides.get("seed"),
        repetitions=overrides.get("nReps"),
        age85_plus_max=int(overrides.get("age85PlusMax") or AGE85_PLUS_MAX_DEFAULT),
        baseline_active_tb_proportion=float(
            overrides.get("baselineActiveTBProportion") or 0.0
        ),
        baseline_active_tb_source=str(
            overrides.get(
                "baselineActiveTBSource",
                "not configured; baseline active TB assigned as zero",
            )
        ),
        baseline_active_tb_review_status=str(
            overrides.get("baselineActiveTBReviewStatus", "not_configured")
        ),
    )
    return {
        **overrides,
        "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
        "age85PlusMax": runner_config["age85PlusMax"],
        "explicitRecentRemoteTBI": runner_config,
    }


def validate_explicit_recent_remote_runner_config(
    payload: Mapping[str, Any],
) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("explicitRecentRemoteTBI must be a mapping.")
    contract = str(
        payload.get("contractVersion", EXPLICIT_RECENT_REMOTE_RUNNER_CONTRACT_VERSION)
    )
    if contract != EXPLICIT_RECENT_REMOTE_RUNNER_CONTRACT_VERSION:
        raise ValueError(
            "explicitRecentRemoteTBI.contractVersion must be "
            f"{EXPLICIT_RECENT_REMOTE_RUNNER_CONTRACT_VERSION!r}."
        )
    pathway = str(payload.get("pathwayIdentifier") or "")
    if pathway != EXPLICIT_RECENT_REMOTE_PATHWAY_ID:
        raise ValueError(
            "explicitRecentRemoteTBI.pathwayIdentifier must be "
            f"{EXPLICIT_RECENT_REMOTE_PATHWAY_ID!r}."
        )
    if not bool(payload.get("enabled", False)):
        raise ValueError("explicitRecentRemoteTBI.enabled must be true when selected.")
    tbi_payload = payload.get("explicitRecentRemoteTBIConfig")
    if not isinstance(tbi_payload, Mapping):
        raise ValueError("explicitRecentRemoteTBIConfig must be supplied.")
    tbi_config = explicit_recent_remote_config_from_dict(tbi_payload)
    if not tbi_config.enabled:
        raise ValueError("explicitRecentRemoteTBIConfig.enabled must be true.")
    if str(tbi_payload.get("configContractVersion")) != TBI_CONFIG_CONTRACT_VERSION:
        raise ValueError("Nested TBI configuration contract version is inconsistent.")
    curve = validate_progression_curve_contract(payload.get("progressionCurve") or {})
    curve_identifier = str(payload.get("progressionCurveIdentifier") or "")
    if curve_identifier != curve["curveIdentifier"]:
        raise ValueError("progressionCurveIdentifier must match progressionCurve.")
    if curve_identifier != CURVE_CENTRAL_TIME_SINCE_INFECTION:
        raise ValueError(
            "Milestone 3B runner integration only supports the central progression curve."
        )
    if str(payload.get("progressionCurveHash") or "") != progression_curve_hash(curve):
        raise ValueError("progressionCurveHash does not match progressionCurve.")
    sampling_method = str(payload.get("progressionTimeSamplingMethod") or "")
    if sampling_method != PROGRESSION_TIME_SAMPLING_METHOD:
        raise ValueError("Unsupported progressionTimeSamplingMethod.")
    recent_window = _positive_float(payload.get("recentWindowYears"), "recentWindowYears")
    if abs(recent_window - tbi_config.recent_window_years) > 1e-12:
        raise ValueError("Runner recentWindowYears must match nested TBI config.")
    if abs(float(curve["recentWindowYears"]) - recent_window) > 1e-12:
        raise ValueError("Progression curve recentWindowYears must match runner config.")
    remote_cap = _positive_float(payload.get("remoteHistoryCapYears"), "remoteHistoryCapYears")
    if abs(remote_cap - tbi_config.remote_history_cap_years) > 1e-12:
        raise ValueError("Runner remoteHistoryCapYears must match nested TBI config.")
    age85_plus_max = _positive_int(payload.get("age85PlusMax"), "age85PlusMax")
    if age85_plus_max != tbi_config.age85_plus_max:
        raise ValueError("Runner age85PlusMax must match nested TBI config.")
    risk_policy = str(payload.get("riskFactorProgressionPolicy") or "")
    if risk_policy != RISK_POLICY_NONE or curve["riskFactorProgressionPolicy"] != RISK_POLICY_NONE:
        raise ValueError("Milestone 3B central runner policy requires riskFactorProgressionPolicy='none'.")
    mortality_policy = str(payload.get("mortalityPolicy") or "")
    if mortality_policy != COMPETING_MORTALITY_NOT_MODELLED:
        raise ValueError("Milestone 3B runner integration does not bundle mortality data.")
    baseline_policy = str(payload.get("baselineActiveTBAssignmentPolicy") or "")
    if baseline_policy != BASELINE_ACTIVE_TB_ASSIGNMENT_POLICY:
        raise ValueError("Unsupported baselineActiveTBAssignmentPolicy.")
    baseline_active = _probability(
        payload.get("baselineActiveTBProportion", 0.0),
        "baselineActiveTBProportion",
    )
    if tbi_config.recent_tbi_target + tbi_config.remote_only_tbi_target + baseline_active > 1.0 + 1e-12:
        raise ValueError(
            "baselineActiveTBProportion plus mutually exclusive recent and "
            "remote-only TBI targets cannot exceed one."
        )
    active_policy = str(payload.get("activeTBObservationPolicy") or "")
    if active_policy != ACTIVE_TB_OBSERVATION_POLICY_VALIDATION_ONLY:
        raise ValueError("Active-TB observations must remain validation-only in Milestone 3B.")
    reinfection_policy = str(payload.get("reinfectionPolicy") or "")
    if reinfection_policy != curve["reinfectionPolicy"]:
        raise ValueError("reinfectionPolicy must match progressionCurve.")
    if reinfection_policy not in {
        REINFECTION_POLICY_RESET_CLOCK,
        REINFECTION_POLICY_NO_RESET_CLOCK,
    }:
        raise ValueError("Unsupported reinfectionPolicy.")
    return deepcopy(dict(payload))


def explicit_recent_remote_runner_config_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(
        validate_explicit_recent_remote_runner_config(payload),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def explicit_recent_remote_runner_config_hash(payload: Mapping[str, Any]) -> str:
    text = explicit_recent_remote_runner_config_json(payload)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def explicit_recent_remote_pathway_fields_present(config: Mapping[str, Any]) -> bool:
    return (
        "analysisPathway" in config
        or "explicitRecentRemoteTBI" in config
        or "explicitRecentRemoteTBIConfig" in config
    )


def is_explicit_recent_remote_pathway_selected(config: Mapping[str, Any]) -> bool:
    return str(config.get("analysisPathway") or "") == EXPLICIT_RECENT_REMOTE_PATHWAY_ID


def validate_explicit_recent_remote_pathway_selection(
    config: Mapping[str, Any],
) -> dict[str, Any] | None:
    pathway = config.get("analysisPathway")
    nested = config.get("explicitRecentRemoteTBI")
    legacy_nested = config.get("explicitRecentRemoteTBIConfig")
    if pathway not in (None, "", EXPLICIT_RECENT_REMOTE_PATHWAY_ID):
        raise ValueError(f"Unknown analysisPathway: {pathway}")
    if pathway in (None, ""):
        if isinstance(nested, Mapping) and bool(nested.get("enabled", False)):
            raise ValueError(
                "explicitRecentRemoteTBI.enabled=true requires "
                f"analysisPathway={EXPLICIT_RECENT_REMOTE_PATHWAY_ID!r}."
            )
        if isinstance(legacy_nested, Mapping) and bool(legacy_nested.get("enabled", False)):
            raise ValueError(
                "explicitRecentRemoteTBIConfig.enabled=true requires "
                f"analysisPathway={EXPLICIT_RECENT_REMOTE_PATHWAY_ID!r}."
            )
        return None
    if str(config.get("analysisBasis") or "") == "sa_health_matlab_v9_compatibility_reference":
        raise ValueError("The explicit pathway must not use the frozen SA Health analysisBasis.")
    if not isinstance(nested, Mapping):
        raise ValueError(
            "analysisPathway='explicit_recent_remote_tbi_v2' requires "
            "explicitRecentRemoteTBI runner configuration."
        )
    return validate_explicit_recent_remote_runner_config(nested)


def build_explicit_recent_remote_runner_calibration(config: Mapping[str, Any]) -> dict[str, Any]:
    runner_config = validate_explicit_recent_remote_pathway_selection(config)
    if runner_config is None:
        raise ValueError("Explicit recent/remote pathway is not selected.")
    runner_config = deepcopy(runner_config)
    runner_config["screeningWindowYears"] = float(
        config.get("screeningWindowYears", config.get("screenWindow", 3.0))
    )
    runner_config["followUpHorizonYears"] = float(
        config.get("followUpHorizonYears", config.get("followHorizon", 20.0))
    )
    runner_config["stochasticSeed"] = int(config["seed"]) if "seed" in config else None
    runner_config["repetitions"] = int(config["nReps"]) if "nReps" in config else None
    cfg_age85_max = int(config.get("age85PlusMax") or AGE85_PLUS_MAX_DEFAULT)
    if cfg_age85_max != int(runner_config["age85PlusMax"]):
        raise ValueError("Top-level age85PlusMax must match explicit runner provenance.")
    pars = load_parameters_from_config(dict(config))
    tbi_config = explicit_recent_remote_config_from_dict(
        runner_config["explicitRecentRemoteTBIConfig"]
    )
    baseline_active = _probability(
        runner_config.get("baselineActiveTBProportion", 0.0),
        "baselineActiveTBProportion",
    )
    tbi_at_risk_fraction = 1.0 - baseline_active
    if tbi_at_risk_fraction <= 0.0 and (
        tbi_config.recent_tbi_target > 0.0 or tbi_config.remote_only_tbi_target > 0.0
    ):
        raise ValueError(
            "Positive TBI targets are infeasible when everyone is assigned "
            "baseline active TB."
        )
    if tbi_at_risk_fraction > 0.0:
        conditional_recent_target = (
            tbi_config.recent_tbi_target / tbi_at_risk_fraction
        )
        conditional_remote_target = (
            tbi_config.remote_only_tbi_target / tbi_at_risk_fraction
        )
    else:
        conditional_recent_target = 0.0
        conditional_remote_target = 0.0
    if conditional_recent_target + conditional_remote_target > 1.0 + 1e-12:
        raise ValueError(
            "Recent and remote-only TBI targets are infeasible after excluding "
            "baseline active TB from the at-risk population."
        )
    calibration = calibrate_recent_remote_hazards(
        pars["exactAgeValues"],
        pars["exactAgeProb"],
        conditional_recent_target,
        conditional_remote_target,
        recent_window_years=tbi_config.recent_window_years,
        remote_history_cap_years=tbi_config.remote_history_cap_years,
    )
    calibration_payload = calibration.as_dict()
    calibration_payload.update(
        {
            "requestedRecentPrevalenceTotalPopulation": tbi_config.recent_tbi_target,
            "requestedRemoteOnlyPrevalenceTotalPopulation": tbi_config.remote_only_tbi_target,
            "baselineActiveTBProportion": baseline_active,
            "tbiAtRiskPopulationProportion": tbi_at_risk_fraction,
            "conditionalRecentTargetAmongNonBaselineActiveTB": conditional_recent_target,
            "conditionalRemoteOnlyTargetAmongNonBaselineActiveTB": conditional_remote_target,
            "achievedRecentPrevalenceTotalPopulation": (
                tbi_at_risk_fraction * calibration.achieved_recent_prevalence
            ),
            "achievedRemoteOnlyPrevalenceTotalPopulation": (
                tbi_at_risk_fraction * calibration.achieved_remote_only_prevalence
            ),
            "achievedTotalTBIPrevalenceTotalPopulation": (
                tbi_at_risk_fraction * calibration.achieved_total_tbi_prevalence
            ),
            "achievedUninfectedPrevalenceTotalPopulation": (
                tbi_at_risk_fraction * calibration.achieved_uninfected_prevalence
            ),
            "baselineStateSequencing": [
                STATE_BASELINE_ACTIVE_TB,
                STATE_RECENT,
                STATE_REMOTE_ONLY,
                STATE_UNINFECTED,
            ],
            "diagnosticMessages": list(calibration_payload.get("diagnosticMessages", ()))
            + [
                "Baseline/prevalent active TB is excluded before TBI state assignment.",
                "Source denominators for observations are retained separately from the at-risk population.",
            ],
        }
    )
    curve = validate_progression_curve_contract(runner_config["progressionCurve"])
    segment_hazards = [float(value) for value in curve["segmentHazards"]]
    return {
        "parameters": pars,
        "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
        "analysisMethod": runner_config["analysisMethod"],
        "calibrationPolicy": "explicit_recent_remote_tbi_v2_runner_calibration",
        "calibrationParametersRetained": [],
        "calibrationParametersRecalibrated": [
            "recentInfectionHazard",
            "remoteOnlyInfectionHazard",
        ],
        "explicitRecentRemoteRunnerConfig": runner_config,
        "explicitRecentRemoteRunnerConfigHash": explicit_recent_remote_runner_config_hash(
            runner_config
        ),
        "explicitRecentRemoteCalibration": calibration_payload,
        "progressionCurve": curve,
        "progressionCurveHash": progression_curve_hash(curve),
        "ageInfLogLambda": -20.0,
        "ageInfGamma": 0.0,
        "expectedInfPrev": calibration_payload["achievedTotalTBIPrevalenceTotalPopulation"],
        "expectedAgeOR": None,
        "lambdaEarly": segment_hazards[0],
        "lambdaLate": float(curve["postFinalAnchorHazard"]),
        "expectedActiveAtCalibrationHorizon": None,
        "expectedActive2y": None,
        "targetInfPrev": calibration_payload["achievedTotalTBIPrevalenceTotalPopulation"],
        "targetAgeOR": None,
        "targetActiveAtCalibrationHorizon": None,
        "targetActive2y": None,
        "baselineRecentLTBIProportion": None,
        "recentToRemoteTransitionRatePerYear": None,
        "ltbiStateAssumptionStatus": RUNNER_INTEGRATION_STATUS,
        "infectionHistory": None,
        "zeroInfectionPrevalence": calibration_payload["achievedTotalTBIPrevalenceTotalPopulation"] == 0.0,
    }


def apply_explicit_recent_remote_population_assignment(
    population: dict[str, np.ndarray],
    calibration: Mapping[str, Any],
    rng: Any,
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    runner_config = validate_explicit_recent_remote_runner_config(
        calibration["explicitRecentRemoteRunnerConfig"]
    )
    tbi_calibration = calibration["explicitRecentRemoteCalibration"]
    recent_hazard = float(tbi_calibration["fittedRecentHazard"])
    remote_hazard = float(tbi_calibration["fittedRemoteHazard"])
    curve = validate_progression_curve_contract(calibration["progressionCurve"])
    ages = np.asarray(population["ageYears"], dtype=float)
    baseline_active_probability = _probability(
        runner_config.get("baselineActiveTBProportion", 0.0),
        "baselineActiveTBProportion",
    )
    baseline_active = (
        rng.random(len(ages)) < baseline_active_probability
        if baseline_active_probability > 0.0
        else np.zeros(len(ages), dtype=bool)
    )
    tbi_eligible = ~baseline_active
    states = np.full(len(ages), STATE_UNINFECTED, dtype=object)
    states[baseline_active] = STATE_BASELINE_ACTIVE_TB
    prior_remote = np.zeros(len(ages), dtype=bool)
    prior_remote_plus_recent = np.zeros(len(ages), dtype=bool)
    time_since = np.full(len(ages), np.nan, dtype=float)
    remaining_early_risk = np.zeros(len(ages), dtype=float)
    eligible_idx = np.flatnonzero(tbi_eligible)
    if len(eligible_idx) > 0:
        assignment = draw_recent_remote_states_for_ages(
            ages[eligible_idx],
            recent_hazard,
            remote_hazard,
            rng=rng,
            recent_window_years=float(runner_config["recentWindowYears"]),
            remote_history_cap_years=float(runner_config["remoteHistoryCapYears"]),
        )
        states[eligible_idx] = np.asarray(assignment["effectiveStates"], dtype=object)
        prior_remote[eligible_idx] = np.asarray(
            assignment["priorRemoteExposure"], dtype=bool
        )
        prior_remote_plus_recent[eligible_idx] = np.asarray(
            assignment["priorRemotePlusRecent"], dtype=bool
        )
        time_since[eligible_idx] = np.asarray(
            [
                np.nan if value is None else float(value)
                for value in assignment["timeSinceMostRecentInfection"]
            ],
            dtype=float,
        )
        remaining_early_risk[eligible_idx] = np.asarray(
            assignment["remainingEarlyRiskYears"], dtype=float
        )
    recent_at_baseline = states == STATE_RECENT
    remote_at_baseline = states == STATE_REMOTE_ONLY
    infected = recent_at_baseline | remote_at_baseline
    remote_clock = _remote_clock_for_no_reset_sensitivity(
        ages,
        prior_remote_plus_recent,
        remote_hazard,
        rng,
        recent_window_years=float(runner_config["recentWindowYears"]),
        remote_history_cap_years=float(runner_config["remoteHistoryCapYears"]),
    )
    progression_clock = _progression_clock_for_assignment(
        states,
        time_since,
        remote_clock,
        prior_remote_plus_recent,
        str(runner_config["reinfectionPolicy"]),
    )
    t_active = _draw_active_times_from_progression_curve(
        states,
        progression_clock,
        curve,
        rng,
    )
    probabilities = mutually_exclusive_state_probabilities(
        ages,
        recent_hazard,
        remote_hazard,
        recent_window_years=float(runner_config["recentWindowYears"]),
        remote_history_cap_years=float(runner_config["remoteHistoryCapYears"]),
    )
    p_recent = np.asarray(probabilities["recent"], dtype=float)
    p_remote_only = np.asarray(probabilities["remote_only"], dtype=float)
    p_tbi = p_recent + p_remote_only
    p_recent[baseline_active] = 0.0
    p_remote_only[baseline_active] = 0.0
    p_tbi[baseline_active] = 0.0

    out = dict(population)
    out["pInfection"] = p_tbi
    out["pRecentGivenInfected"] = np.divide(
        p_recent,
        p_tbi,
        out=np.zeros_like(p_tbi, dtype=float),
        where=p_tbi > 0.0,
    )
    out["infected"] = infected
    out["explicitTBIState"] = states
    out["baselineActiveTBAtBaseline"] = baseline_active
    out["tbiEligibleAtBaseline"] = tbi_eligible
    out["priorRemoteExposure"] = prior_remote
    out["priorRemotePlusRecent"] = prior_remote_plus_recent
    out["timeSinceMostRecentInfection"] = time_since
    out["timeSinceRemoteInfectionNoReset"] = remote_clock
    out["progressionClockYears"] = progression_clock
    out["remainingEarlyRiskYears"] = remaining_early_risk
    out["diseaseMultiplier"] = np.ones_like(ages, dtype=float)
    _attach_explicit_targeting_scores(out, curve, runner_config)
    transition_time = np.full(len(ages), np.inf, dtype=float)
    transition_time[recent_at_baseline] = out["remainingEarlyRiskYears"][
        recent_at_baseline
    ]
    return out, recent_at_baseline, remote_at_baseline, transition_time, t_active


def is_explicit_recent_remote_calibration(calibration: Mapping[str, Any]) -> bool:
    return str(calibration.get("analysisPathway") or "") == EXPLICIT_RECENT_REMOTE_PATHWAY_ID


def sample_future_active_tb_time_from_curve(
    curve: Mapping[str, Any],
    *,
    time_since_infection_at_baseline: float,
    rng: Any,
    tolerance: float = 1e-12,
) -> float:
    validated = validate_progression_curve_contract(curve)
    u = float(rng.random())
    if u <= 0.0:
        return 0.0
    if u >= 1.0:
        return math.inf
    start_time = _finite_nonnegative_float(
        time_since_infection_at_baseline,
        "time_since_infection_at_baseline",
    )
    target_hazard = time_since_curve_cumulative_hazard(
        validated,
        start_time,
    ) - math.log1p(-u)
    lo = 0.0
    hi = 1.0
    while time_since_curve_cumulative_hazard(validated, start_time + hi) < target_hazard:
        hi *= 2.0
        if hi > 1.0e9:
            return math.inf
    for _ in range(160):
        mid = (lo + hi) / 2.0
        value = time_since_curve_cumulative_hazard(validated, start_time + mid)
        if abs(value - target_hazard) <= tolerance or (hi - lo) / 2.0 <= tolerance:
            return mid
        if value < target_hazard:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def active_tb_validation_only_for_runner(
    observation: Mapping[str, Any],
    strata: list[Mapping[str, Any]],
    calibration: Mapping[str, Any],
    *,
    model_baseline_year: int,
    ascertainment_probability: float,
    ascertainment_source: str,
    ascertainment_review_status: str,
    observation_model_id: str | None = None,
) -> dict[str, Any]:
    from engine.apy.explicit_recent_remote_progression import (
        validate_active_tb_observation_against_progression_curve,
    )

    curve = validate_progression_curve_contract(calibration["progressionCurve"])
    result = validate_active_tb_observation_against_progression_curve(
        observation,
        strata,
        curve,
        model_baseline_year=model_baseline_year,
        ascertainment_probability=ascertainment_probability,
        ascertainment_source=ascertainment_source,
        ascertainment_review_status=ascertainment_review_status,
        observation_model_id=observation_model_id,
    )
    result["activeTBObservationPolicy"] = ACTIVE_TB_OBSERVATION_POLICY_VALIDATION_ONLY
    result["progressionCurveHash"] = progression_curve_hash(curve)
    return result


def explicit_recent_remote_event_ledger_unavailable(
    runner_config: Mapping[str, Any],
) -> dict[str, Any]:
    validated = validate_explicit_recent_remote_runner_config(runner_config)
    return {
        "available": False,
        "source": EXPLICIT_EVENT_LEDGER_OUT_OF_SCOPE_SOURCE,
        "metadata": {
            "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
            "contractVersion": validated["contractVersion"],
            "explicitRecentRemoteRunnerConfigHash": explicit_recent_remote_runner_config_hash(
                validated
            ),
            "progressionCurveIdentifier": validated["progressionCurveIdentifier"],
            "progressionCurveHash": validated["progressionCurveHash"],
            "integrationStatus": RUNNER_INTEGRATION_STATUS,
        },
        "definitions": [],
        "replicateTotals": [],
        "annualEvents": [],
        "summaries": {},
        "validation": {
            "isValid": True,
            "errors": [],
            "warnings": [
                "Event-ledger generation, economics and DALYs are intentionally out of scope for the explicit recent/remote pathway in Milestone 3B."
            ],
        },
    }


def deterministic_explicit_recent_remote_expectations(
    calibration: Mapping[str, Any],
    *,
    horizons_years: tuple[float, ...] = (1.0, 5.0, 10.0, 20.0),
    population_size: float = 1.0,
    quadrature_points: int = 48,
    reinfection_policy: str | None = None,
    competing_survival_curve: Any = None,
    competing_mortality_hazard: Any = None,
    integration_steps: int = 512,
) -> dict[str, Any]:
    if not is_explicit_recent_remote_calibration(calibration):
        raise ValueError("deterministic expectations require explicit recent/remote calibration.")
    runner_config = validate_explicit_recent_remote_runner_config(
        calibration["explicitRecentRemoteRunnerConfig"]
    )
    curve = validate_progression_curve_contract(calibration["progressionCurve"])
    policy = str(reinfection_policy or runner_config["reinfectionPolicy"])
    if policy not in {REINFECTION_POLICY_RESET_CLOCK, REINFECTION_POLICY_NO_RESET_CLOCK}:
        raise ValueError("Unsupported reinfection_policy.")
    if competing_survival_curve is not None and competing_mortality_hazard is not None:
        raise ValueError(
            "Provide either competing_survival_curve or competing_mortality_hazard, not both."
        )
    horizons = tuple(_finite_nonnegative_float(value, "horizons_years") for value in horizons_years)
    pop_size = _finite_nonnegative_float(population_size, "population_size")
    q_points = _positive_int(quadrature_points, "quadrature_points")
    recent_window = float(runner_config["recentWindowYears"])
    remote_cap = float(runner_config["remoteHistoryCapYears"])
    baseline_active = _probability(
        runner_config.get("baselineActiveTBProportion", 0.0),
        "baselineActiveTBProportion",
    )
    at_risk_fraction = 1.0 - baseline_active
    tbi_calibration = calibration["explicitRecentRemoteCalibration"]
    recent_hazard = float(tbi_calibration["fittedRecentHazard"])
    remote_hazard = float(tbi_calibration["fittedRemoteHazard"])
    pars = calibration["parameters"]
    ages = np.asarray(pars["exactAgeValues"], dtype=float)
    weights = np.asarray(pars["exactAgeProb"], dtype=float)
    weights = weights / weights.sum()
    probabilities = mutually_exclusive_state_probabilities(
        ages,
        recent_hazard,
        remote_hazard,
        recent_window_years=recent_window,
        remote_history_cap_years=remote_cap,
    )
    durations = exposure_durations_for_ages(
        ages,
        recent_window_years=recent_window,
        remote_history_cap_years=remote_cap,
    )
    p_recent = np.asarray(probabilities["recent"], dtype=float)
    p_remote = np.asarray(probabilities["remote_only"], dtype=float)
    p_uninfected = np.asarray(probabilities["uninfected"], dtype=float)
    p_recent_with_prior = np.asarray(
        probabilities["recent_with_prior_remote"], dtype=float
    )
    p_recent_without_prior = np.asarray(
        probabilities["recent_without_prior_remote"], dtype=float
    )
    state_proportions = {
        STATE_BASELINE_ACTIVE_TB: baseline_active,
        STATE_RECENT: at_risk_fraction * float(np.sum(weights * p_recent)),
        STATE_REMOTE_ONLY: at_risk_fraction * float(np.sum(weights * p_remote)),
        STATE_UNINFECTED: at_risk_fraction * float(np.sum(weights * p_uninfected)),
        "priorRemotePlusRecent": at_risk_fraction
        * float(np.sum(weights * p_recent_with_prior)),
    }
    state_counts = {
        key: value * pop_size for key, value in state_proportions.items()
    }
    age_rows = []
    horizon_rows: dict[str, dict[str, Any]] = {}
    for horizon in horizons:
        recent_total = 0.0
        remote_total = 0.0
        reset_total = 0.0
        for idx, age in enumerate(ages):
            recent_duration = float(durations.recent_years[idx])
            remote_duration = float(durations.remote_years[idx])
            recent_expectation = _conditional_time_expectation(
                window_start_years_before_baseline=0.0,
                window_duration_years=recent_duration,
                hazard_per_year=recent_hazard,
                quadrature_points=q_points,
                fn=lambda clock, age_value=float(age): _progression_probability_for_clock(
                    curve,
                    clock,
                    horizon,
                    age_value,
                    competing_survival_curve=competing_survival_curve,
                    competing_mortality_hazard=competing_mortality_hazard,
                    integration_steps=integration_steps,
                ),
            )
            remote_expectation = _conditional_time_expectation(
                window_start_years_before_baseline=recent_window,
                window_duration_years=remote_duration,
                hazard_per_year=remote_hazard,
                quadrature_points=q_points,
                fn=lambda clock, age_value=float(age): _progression_probability_for_clock(
                    curve,
                    clock,
                    horizon,
                    age_value,
                    competing_survival_curve=competing_survival_curve,
                    competing_mortality_hazard=competing_mortality_hazard,
                    integration_steps=integration_steps,
                ),
            )
            if policy == REINFECTION_POLICY_NO_RESET_CLOCK:
                recent_component = (
                    float(p_recent_without_prior[idx]) * recent_expectation
                    + float(p_recent_with_prior[idx]) * remote_expectation
                )
            else:
                recent_component = float(p_recent[idx]) * recent_expectation
            remote_component = float(p_remote[idx]) * remote_expectation
            recent_total += float(weights[idx]) * recent_component
            remote_total += float(weights[idx]) * remote_component
            reset_total += float(weights[idx]) * float(p_recent_with_prior[idx]) * (
                recent_expectation - remote_expectation
            )
        proportion = at_risk_fraction * (recent_total + remote_total)
        horizon_rows[f"{horizon:g}"] = {
            "horizonYears": horizon,
            "incidentActiveTBProportion": proportion,
            "incidentActiveTBCount": pop_size * proportion,
            "recentContributionProportion": at_risk_fraction * recent_total,
            "remoteOnlyContributionProportion": at_risk_fraction * remote_total,
            "resetMinusNoResetRecentPriorRemoteContribution": at_risk_fraction
            * reset_total,
        }
    remaining_recent = 0.0
    for idx, age in enumerate(ages):
        recent_duration = float(durations.recent_years[idx])
        remaining_recent += float(weights[idx]) * float(p_recent[idx]) * _conditional_time_expectation(
            window_start_years_before_baseline=0.0,
            window_duration_years=recent_duration,
            hazard_per_year=recent_hazard,
            quadrature_points=q_points,
            fn=lambda clock: max(0.0, recent_window - clock),
        )
        age_rows.append(
            {
                "ageYears": float(age),
                "ageWeight": float(weights[idx]),
                "baselineActiveTBProbability": baseline_active,
                "recentProbabilityTotalPopulation": at_risk_fraction
                * float(p_recent[idx]),
                "remoteOnlyProbabilityTotalPopulation": at_risk_fraction
                * float(p_remote[idx]),
                "uninfectedProbabilityTotalPopulation": at_risk_fraction
                * float(p_uninfected[idx]),
                "priorRemotePlusRecentProbabilityTotalPopulation": at_risk_fraction
                * float(p_recent_with_prior[idx]),
            }
        )
    recent_denominator = float(np.sum(weights * p_recent))
    expected_remaining = 0.0 if recent_denominator <= 0.0 else remaining_recent / recent_denominator
    return {
        "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
        "method": "deterministic_expected_value",
        "integrationMethod": DETERMINISTIC_INFECTION_TIME_INTEGRATION_METHOD,
        "quadraturePoints": q_points,
        "usesMeanTimeSinceInfection": False,
        "reinfectionPolicy": policy,
        "mortalityPolicy": (
            COMPETING_MORTALITY_CONSTANT_HAZARD
            if competing_mortality_hazard is not None
            else (
                "external_survival_curve"
                if competing_survival_curve is not None
                else COMPETING_MORTALITY_NOT_MODELLED
            )
        ),
        "baselineStateSequencing": [
            STATE_BASELINE_ACTIVE_TB,
            STATE_RECENT,
            STATE_REMOTE_ONLY,
            STATE_UNINFECTED,
        ],
        "stateProportions": state_proportions,
        "stateCounts": state_counts,
        "ageSpecificStateDistributions": age_rows,
        "incidentActiveTBByHorizon": horizon_rows,
        "expectedRemainingEarlyRiskYearsAmongRecentlyInfected": expected_remaining,
        "diagnosticMessages": [
            "Deterministic expectations integrate over the conditional infection-time distribution.",
            "Baseline/prevalent active TB is excluded before TBI assignment and prospective incident progression.",
        ],
    }


def stochastic_explicit_recent_remote_natural_history(
    calibration: Mapping[str, Any],
    *,
    population_size: int,
    horizons_years: tuple[float, ...] = (1.0, 5.0, 10.0, 20.0),
    seed: int | None = None,
    rng: Any = None,
    competing_mortality_hazard: Any = None,
    competing_survival_curve: Any = None,
) -> dict[str, Any]:
    if rng is None and seed is None:
        raise ValueError("A seed or explicit rng is required for stochastic explicit runner summaries.")
    if rng is not None and seed is not None:
        raise ValueError("Provide either seed or rng, not both.")
    if competing_mortality_hazard is not None and competing_survival_curve is not None:
        raise ValueError(
            "Provide either competing_mortality_hazard or competing_survival_curve, not both."
        )
    generator = rng if rng is not None else np.random.default_rng(int(seed))
    n = _positive_int(population_size, "population_size")
    horizons = tuple(_finite_nonnegative_float(value, "horizons_years") for value in horizons_years)
    pars = calibration["parameters"]
    age_values = np.asarray(pars["exactAgeValues"], dtype=float)
    age_prob = np.asarray(pars["exactAgeProb"], dtype=float)
    age_prob = age_prob / age_prob.sum()
    ages = generator.choice(age_values, size=n, p=age_prob)
    population = {"ageYears": ages}
    assigned, recent, remote, _transition, active_time = apply_explicit_recent_remote_population_assignment(
        population,
        calibration,
        generator,
    )
    baseline_active = np.asarray(assigned["baselineActiveTBAtBaseline"], dtype=bool)
    death_time = np.full(n, np.inf, dtype=float)
    mortality_draws_used = False
    if competing_mortality_hazard is not None:
        mortality_draws_used = True
        death_time = _draw_death_times_from_hazard(
            competing_mortality_hazard,
            ages,
            generator,
            max(horizons) if horizons else 0.0,
        )
    elif competing_survival_curve is not None:
        mortality_draws_used = True
        death_time = _draw_death_times_from_survival_curve(
            competing_survival_curve,
            ages,
            generator,
            max(horizons) if horizons else 0.0,
        )
    by_horizon = {}
    at_risk_for_incident = ~baseline_active
    for horizon in horizons:
        incident = (
            at_risk_for_incident
            & (active_time <= horizon)
            & (active_time < death_time)
        )
        by_horizon[f"{horizon:g}"] = {
            "horizonYears": horizon,
            "incidentActiveTBCount": int(incident.sum()),
            "incidentActiveTBProportion": float(incident.sum() / n),
        }
    return {
        "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
        "method": STOCHASTIC_NATURAL_HISTORY_METHOD,
        "populationSize": n,
        "seed": None if seed is None else int(seed),
        "mortalityDrawsUsed": mortality_draws_used,
        "stateCounts": {
            STATE_BASELINE_ACTIVE_TB: int(baseline_active.sum()),
            STATE_RECENT: int(recent.sum()),
            STATE_REMOTE_ONLY: int(remote.sum()),
            STATE_UNINFECTED: int(
                (
                    np.asarray(assigned["explicitTBIState"], dtype=object)
                    == STATE_UNINFECTED
                ).sum()
            ),
            "priorRemotePlusRecent": int(
                np.asarray(assigned["priorRemotePlusRecent"], dtype=bool).sum()
            ),
        },
        "stateProportions": {
            key: value / n
            for key, value in {
                STATE_BASELINE_ACTIVE_TB: int(baseline_active.sum()),
                STATE_RECENT: int(recent.sum()),
                STATE_REMOTE_ONLY: int(remote.sum()),
                STATE_UNINFECTED: int(
                    (
                        np.asarray(assigned["explicitTBIState"], dtype=object)
                        == STATE_UNINFECTED
                    ).sum()
                ),
                "priorRemotePlusRecent": int(
                    np.asarray(assigned["priorRemotePlusRecent"], dtype=bool).sum()
                ),
            }.items()
        },
        "incidentActiveTBByHorizon": by_horizon,
        "audit": {
            "baselineActiveTBExcludedFromIncidentRisk": True,
            "globalRandomStateUsed": False,
            "progressionTimeSamplingMethod": PROGRESSION_TIME_SAMPLING_METHOD,
        },
    }


def reconcile_deterministic_stochastic_explicit_runner(
    calibration: Mapping[str, Any],
    *,
    population_size: int = 50000,
    seed: int = 20261011,
    horizons_years: tuple[float, ...] = (1.0, 5.0, 10.0, 20.0),
    z_score: float = 4.0,
    absolute_floor: float = 0.003,
) -> dict[str, Any]:
    deterministic = deterministic_explicit_recent_remote_expectations(
        calibration,
        horizons_years=horizons_years,
        population_size=population_size,
    )
    stochastic = stochastic_explicit_recent_remote_natural_history(
        calibration,
        population_size=population_size,
        seed=seed,
        horizons_years=horizons_years,
    )
    rows = []
    n = float(population_size)
    for key, expected in deterministic["stateProportions"].items():
        if key == "priorRemotePlusRecent":
            continue
        observed = float(stochastic["stateProportions"][key])
        tolerance = max(absolute_floor, z_score * math.sqrt(max(expected * (1.0 - expected), 0.0) / n))
        rows.append(
            {
                "quantity": f"state:{key}",
                "deterministic": expected,
                "stochastic": observed,
                "absoluteDifference": abs(observed - expected),
                "tolerance": tolerance,
                "withinTolerance": abs(observed - expected) <= tolerance,
            }
        )
    for key, expected_row in deterministic["incidentActiveTBByHorizon"].items():
        expected = float(expected_row["incidentActiveTBProportion"])
        observed = float(
            stochastic["incidentActiveTBByHorizon"][key][
                "incidentActiveTBProportion"
            ]
        )
        tolerance = max(absolute_floor, z_score * math.sqrt(max(expected * (1.0 - expected), 0.0) / n))
        rows.append(
            {
                "quantity": f"incidentActiveTB:{key}y",
                "deterministic": expected,
                "stochastic": observed,
                "absoluteDifference": abs(observed - expected),
                "tolerance": tolerance,
                "withinTolerance": abs(observed - expected) <= tolerance,
            }
        )
    return {
        "method": "deterministic_stochastic_reconciliation_v1",
        "populationSize": int(population_size),
        "seed": int(seed),
        "rows": rows,
        "allWithinTolerance": all(row["withinTolerance"] for row in rows),
        "deterministic": deterministic,
        "stochastic": stochastic,
    }


def _build_progression_curve(
    curve_identifier: str,
    *,
    reinfection_policy: str,
    mortality_policy: str,
) -> dict[str, Any]:
    if str(curve_identifier) != CURVE_CENTRAL_TIME_SINCE_INFECTION:
        return build_progression_curve_by_identifier(str(curve_identifier))
    return build_central_time_since_infection_progression_curve(
        reinfection_policy=reinfection_policy,
        mortality_policy=mortality_policy,
    )


def _attach_explicit_targeting_scores(
    population: dict[str, np.ndarray],
    curve: Mapping[str, Any],
    runner_config: Mapping[str, Any],
) -> None:
    p_tbi = np.asarray(population["pInfection"], dtype=float)
    infected = np.asarray(population["infected"], dtype=bool)
    clocks = np.asarray(population["progressionClockYears"], dtype=float)
    follow_horizon = float(
        runner_config.get("followUpHorizonYears")
        or runner_config.get("followHorizon")
        or 20.0
    )
    future_risk = np.zeros_like(p_tbi, dtype=float)
    for idx in np.flatnonzero(infected):
        future_risk[idx] = time_since_curve_conditional_future_progression_probability(
            curve,
            time_since_infection_at_baseline=float(clocks[idx]),
            horizon_years=follow_horizon,
        )
    population["ltbiRiskScore"] = p_tbi
    population["cureTargetScore"] = p_tbi
    population["preventTargetScore"] = p_tbi * future_risk


def _draw_active_times_from_progression_curve(
    states: np.ndarray,
    progression_clock: np.ndarray,
    curve: Mapping[str, Any],
    rng: Any,
) -> np.ndarray:
    t_active = np.full(len(states), np.inf, dtype=float)
    infected_idx = np.flatnonzero(
        (states == STATE_RECENT) | (states == STATE_REMOTE_ONLY)
    )
    for idx in infected_idx:
        t_active[idx] = sample_future_active_tb_time_from_curve(
            curve,
            time_since_infection_at_baseline=float(progression_clock[idx]),
            rng=rng,
        )
    return t_active


def _progression_clock_for_assignment(
    states: np.ndarray,
    time_since: np.ndarray,
    remote_clock: np.ndarray,
    prior_remote_plus_recent: np.ndarray,
    reinfection_policy: str,
) -> np.ndarray:
    clock = np.full(len(states), np.nan, dtype=float)
    infected = (states == STATE_RECENT) | (states == STATE_REMOTE_ONLY)
    clock[infected] = time_since[infected]
    if reinfection_policy == REINFECTION_POLICY_NO_RESET_CLOCK:
        no_reset = prior_remote_plus_recent & np.isfinite(remote_clock)
        clock[no_reset] = remote_clock[no_reset]
    return clock


def _remote_clock_for_no_reset_sensitivity(
    ages: np.ndarray,
    prior_remote_plus_recent: np.ndarray,
    remote_hazard: float,
    rng: Any,
    *,
    recent_window_years: float,
    remote_history_cap_years: float,
) -> np.ndarray:
    clock = np.full(len(ages), np.nan, dtype=float)
    if remote_hazard <= 0.0 or not np.any(prior_remote_plus_recent):
        return clock
    durations = exposure_durations_for_ages(
        ages,
        recent_window_years=recent_window_years,
        remote_history_cap_years=remote_history_cap_years,
    )
    remote_years = np.asarray(durations.remote_years, dtype=float)
    for idx in np.flatnonzero(prior_remote_plus_recent & (remote_years > 0.0)):
        clock[idx] = sample_time_since_most_recent_event(
            window_start_years_before_baseline=recent_window_years,
            window_duration_years=float(remote_years[idx]),
            hazard_per_year=remote_hazard,
            rng=rng,
        )
    return clock


def _progression_probability_for_clock(
    curve: Mapping[str, Any],
    clock: float,
    horizon: float,
    age: float,
    *,
    competing_survival_curve: Any,
    competing_mortality_hazard: Any,
    integration_steps: int,
) -> float:
    if competing_survival_curve is None and competing_mortality_hazard is None:
        return time_since_curve_conditional_future_progression_probability(
            curve,
            time_since_infection_at_baseline=clock,
            horizon_years=horizon,
        )
    if competing_survival_curve is not None:
        result = competing_risk_progression_probability_for_curve(
            curve,
            time_since_infection_at_baseline=clock,
            horizon_years=horizon,
            competing_survival_curve=_age_wrapped_survival_curve(
                competing_survival_curve,
                age,
            ),
            integration_steps=integration_steps,
        )
        return float(result["cumulativeIncidence"])
    if callable(competing_mortality_hazard):
        result = competing_risk_progression_probability_for_curve(
            curve,
            time_since_infection_at_baseline=clock,
            horizon_years=horizon,
            competing_survival_curve=_hazard_to_survival_curve(
                competing_mortality_hazard,
                age,
            ),
            integration_steps=integration_steps,
        )
        return float(result["cumulativeIncidence"])
    result = competing_risk_progression_probability_for_curve(
        curve,
        time_since_infection_at_baseline=clock,
        horizon_years=horizon,
        competing_mortality_hazard=competing_mortality_hazard,
        integration_steps=integration_steps,
    )
    return float(result["cumulativeIncidence"])


def _conditional_time_expectation(
    *,
    window_start_years_before_baseline: float,
    window_duration_years: float,
    hazard_per_year: float,
    quadrature_points: int,
    fn,
) -> float:
    start = _finite_nonnegative_float(
        window_start_years_before_baseline,
        "window_start_years_before_baseline",
    )
    duration = _finite_nonnegative_float(
        window_duration_years,
        "window_duration_years",
    )
    hazard = _finite_nonnegative_float(hazard_per_year, "hazard_per_year")
    if duration <= 0.0:
        return 0.0
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_points)
    years_after_start = 0.5 * duration * (nodes + 1.0)
    clocks = start + years_after_start
    if hazard <= 0.0:
        densities = np.full_like(years_after_start, 1.0 / duration, dtype=float)
    else:
        normaliser = -math.expm1(-hazard * duration)
        densities = hazard * np.exp(-hazard * years_after_start) / normaliser
    values = np.asarray([float(fn(float(clock))) for clock in clocks], dtype=float)
    return float(0.5 * duration * np.sum(weights * densities * values))


def _age_wrapped_survival_curve(survival_curve: Any, age: float):
    if not callable(survival_curve):
        return survival_curve

    def survival_at(t: float) -> float:
        return _call_age_time_function(survival_curve, t, age)

    return survival_at


def _hazard_to_survival_curve(hazard_function: Any, age: float):
    def survival_at(t: float) -> float:
        horizon = _finite_nonnegative_float(t, "mortality_horizon")
        if horizon <= 0.0:
            return 1.0
        steps = max(32, int(math.ceil(horizon * 16)))
        grid = np.linspace(0.0, horizon, steps + 1)
        hazards = np.asarray(
            [
                _finite_nonnegative_float(
                    _call_age_time_function(hazard_function, float(time), age),
                    "competing_mortality_hazard",
                )
                for time in grid
            ],
            dtype=float,
        )
        cumulative = float(np.trapz(hazards, grid))
        return math.exp(-cumulative)

    return survival_at


def _call_age_time_function(fn, t: float, age: float) -> float:
    try:
        return float(fn(t, {"ageYears": age}))
    except TypeError:
        try:
            return float(fn(t, age))
        except TypeError:
            return float(fn(t))


def _draw_death_times_from_hazard(
    mortality_hazard: Any,
    ages: np.ndarray,
    rng: Any,
    max_horizon: float,
) -> np.ndarray:
    if mortality_hazard is None:
        return np.full(len(ages), np.inf, dtype=float)
    if not callable(mortality_hazard):
        hazard = _finite_nonnegative_float(
            mortality_hazard,
            "competing_mortality_hazard",
        )
        if hazard <= 0.0:
            return np.full(len(ages), np.inf, dtype=float)
        return rng.exponential(scale=1.0 / hazard, size=len(ages))
    death_times = np.full(len(ages), np.inf, dtype=float)
    for idx, age in enumerate(ages):
        target = -math.log1p(-float(rng.random()))
        hi = max(max_horizon, 1.0)
        survival = _hazard_to_survival_curve(mortality_hazard, float(age))
        while -math.log(max(survival(hi), 1e-300)) < target and hi < 1e4:
            hi *= 2.0
        if hi >= 1e4:
            continue
        lo = 0.0
        for _ in range(80):
            mid = (lo + hi) / 2.0
            cumulative = -math.log(max(survival(mid), 1e-300))
            if cumulative < target:
                lo = mid
            else:
                hi = mid
        death_times[idx] = (lo + hi) / 2.0
    return death_times


def _draw_death_times_from_survival_curve(
    survival_curve: Any,
    ages: np.ndarray,
    rng: Any,
    max_horizon: float,
) -> np.ndarray:
    death_times = np.full(len(ages), np.inf, dtype=float)
    for idx, age in enumerate(ages):
        survival = _age_wrapped_survival_curve(survival_curve, float(age))
        target = float(rng.random())
        if float(survival(max_horizon)) > target:
            continue
        lo = 0.0
        hi = max_horizon
        for _ in range(80):
            mid = (lo + hi) / 2.0
            if float(survival(mid)) > target:
                lo = mid
            else:
                hi = mid
        death_times[idx] = (lo + hi) / 2.0
    return death_times


def _probability(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0 or number > 1.0:
        raise ValueError(f"{label} must be a finite probability in [0,1].")
    return number


def _positive_float(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{label} must be finite and positive.")
    return number


def _positive_int(value: Any, label: str) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a positive integer.") from exc
    if number <= 0:
        raise ValueError(f"{label} must be a positive integer.")
    return number


def _finite_nonnegative_float(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{label} must be finite and non-negative.")
    return number
