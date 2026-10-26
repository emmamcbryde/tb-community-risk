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
from engine.apy.regimen import regimen_partial_efficacy
from engine.apy.population_tbi_progression_risk import (
    ADJUSTMENT_FIELD, adjustment_for_runner, build_progression_risk_adjustment,
    positive_multiplier, validate_progression_risk_adjustment,
)


EXPLICIT_RECENT_REMOTE_PATHWAY_ID = "explicit_recent_remote_tbi_v2"
EXPLICIT_RECENT_REMOTE_RUNNER_CONTRACT_VERSION = (
    "explicit_recent_remote_tbi_runner_config_v1"
)
EXPLICIT_RECENT_REMOTE_ANALYSIS_METHOD = (
    "isolated_agent_based_explicit_recent_remote_tbi_v2"
)
ACTIVE_TB_OBSERVATION_POLICY_VALIDATION_ONLY = "active_tb_validation_only_v1"
RUNNER_INTEGRATION_STATUS = "milestone_3c_runner_hardened_intervention_integration"
RUNNER_HARDENING_STATUS = "milestone_3c_runner_hardened_intervention_integration"
PROGRESSION_TIME_SAMPLING_METHOD = (
    "inverse_survivor_conditioned_time_since_infection_cumulative_hazard_v1"
)
BASELINE_ACTIVE_TB_ASSIGNMENT_POLICY = (
    "baseline_active_tb_separate_before_tbi_assignment_v1"
)
EXPLICIT_EVENT_LEDGER_OUT_OF_SCOPE_SOURCE = (
    "explicit_recent_remote_tbi_v2_event_ledger_out_of_scope_milestone_3c"
)
DETERMINISTIC_INFECTION_TIME_INTEGRATION_METHOD = (
    "gauss_legendre_quadrature_over_conditional_infection_time_v1"
)
STOCHASTIC_NATURAL_HISTORY_METHOD = (
    "isolated_explicit_recent_remote_natural_history_stochastic_v1"
)
STOCHASTIC_INTERVENTION_METHOD = (
    "isolated_explicit_recent_remote_intervention_stochastic_v1"
)
INTERVENTION_RECONCILIATION_METHOD = (
    "explicit_recent_remote_intervention_reconciliation_with_mcse_v1"
)
SCREENING_COVERAGE_POLICY_IDENTIFIER = "programme_coverage_at_least_once_v1"
SCREENING_COVERAGE_DEFINITION = (
    "The proportion of the eligible population screened at least once over "
    "the complete screening programme window."
)
TPT_EFFECT_POLICY_CONTRACT_VERSION = "explicit_recent_remote_tpt_effect_policy_v1"
TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION = (
    "effective_tpt_prevents_future_progression_v1"
)
COMMON_RANDOM_NUMBER_CONTRACT_VERSION = (
    "explicit_recent_remote_common_random_number_v1"
)
COMMON_RANDOM_NUMBER_PAIRED_BASELINE = (
    "paired_baseline_untreated_progression_shared_v1"
)
EPIDEMIOLOGICAL_INTERVENTION_OUTPUT_CONTRACT_VERSION = (
    "explicit_recent_remote_epidemiological_intervention_outputs_v1"
)
DETERMINISTIC_INTERVENTION_INTEGRATION_METHOD = (
    "gauss_legendre_infection_and_screening_time_quadrature_v1"
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
    tpt_effect_policy_identifier: str = TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION,
    common_random_number_contract_identifier: str = COMMON_RANDOM_NUMBER_PAIRED_BASELINE,
    progression_risk_adjustment: Mapping[str, Any] | None = None,
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
    tpt_effect_policy = build_tpt_effect_policy_contract(
        policy_identifier=tpt_effect_policy_identifier,
        review_status="milestone_3c_matches_inherited_stochastic_semantics",
    )
    common_random_number_contract = build_common_random_number_contract(
        contract_identifier=common_random_number_contract_identifier,
    )
    payload = {
        ADJUSTMENT_FIELD: build_progression_risk_adjustment() if progression_risk_adjustment is None else validate_progression_risk_adjustment(progression_risk_adjustment),
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
        "tptEffectPolicy": tpt_effect_policy,
        "tptEffectPolicyHash": tpt_effect_policy_hash(tpt_effect_policy),
        "tptEffectPolicyIdentifier": tpt_effect_policy["policyIdentifier"],
        "commonRandomNumberContract": common_random_number_contract,
        "commonRandomNumberContractHash": common_random_number_contract_hash(
            common_random_number_contract
        ),
        "commonRandomNumberContractIdentifier": common_random_number_contract[
            "contractIdentifier"
        ],
        "screeningCoveragePolicyIdentifier": SCREENING_COVERAGE_POLICY_IDENTIFIER,
        "screeningCoverageDefinition": SCREENING_COVERAGE_DEFINITION,
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
                "source": "Milestone 3C runner-hardening working decision",
                "reviewStatus": "use_none_initially",
                "notes": "Inherited disease ORs are not applied by the central explicit pathway.",
            },
            "mortality": {
                "policy": mortality_policy,
                "source": "No external mortality table supplied in Milestone 3C.",
                "reviewStatus": "not_modelled",
            },
            "activeTBObservations": {
                "policy": active_tb_observation_policy,
                "reviewStatus": "validation_only_by_default",
            },
            "tptEffect": {
                "policy": tpt_effect_policy["policyIdentifier"],
                "reviewStatus": tpt_effect_policy["reviewStatus"],
                "notes": tpt_effect_policy["notes"],
            },
            "commonRandomNumbers": {
                "contractIdentifier": common_random_number_contract[
                    "contractIdentifier"
                ],
                "reviewStatus": common_random_number_contract["reviewStatus"],
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
    progression_risk_adjustment: Mapping[str, Any] | None = None,
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
        progression_risk_adjustment=progression_risk_adjustment,
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


def build_tpt_effect_policy_contract(
    *,
    policy_identifier: str = TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION,
    review_status: str = "milestone_3c_policy_contract",
    notes: str = "",
) -> dict[str, Any]:
    if policy_identifier != TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION:
        raise ValueError("Unsupported explicit recent/remote TPT effect policy.")
    payload = {
        "contractVersion": TPT_EFFECT_POLICY_CONTRACT_VERSION,
        "policyIdentifier": policy_identifier,
        "effectTiming": "protection_begins_at_treatment_stop_time",
        "historicalInfectionStateRetained": True,
        "changesHistoricalInfectionPrevalence": False,
        "preventsOnlyFutureProgression": True,
        "falsePositiveEpidemiologicalBenefit": 0.0,
        "baselineActiveTBEligibleForOrdinaryTPT": False,
        "incompleteOrIneffectiveTreatment": "untreated_survivor_conditioned_progression_curve_retained",
        "matchesInheritedImplementation": True,
        "reviewStatus": str(review_status),
        "notes": str(
            notes
            or "Inherited stochastic APY semantics: a successful TPT efficacy draw prevents only active TB whose untreated event time occurs after the treatment stop time."
        ),
    }
    return validate_tpt_effect_policy_contract(payload)


def validate_tpt_effect_policy_contract(payload: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("tptEffectPolicy must be a mapping.")
    if str(payload.get("contractVersion") or "") != TPT_EFFECT_POLICY_CONTRACT_VERSION:
        raise ValueError("Unsupported tptEffectPolicy contractVersion.")
    if str(payload.get("policyIdentifier") or "") != TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION:
        raise ValueError("Unsupported tptEffectPolicy policyIdentifier.")
    if float(payload.get("falsePositiveEpidemiologicalBenefit", 1.0)) != 0.0:
        raise ValueError("False-positive TPT epidemiological benefit must be zero.")
    if bool(payload.get("baselineActiveTBEligibleForOrdinaryTPT", True)):
        raise ValueError("Baseline active TB must not be eligible for ordinary TPT.")
    return deepcopy(dict(payload))


def tpt_effect_policy_json(payload: Mapping[str, Any]) -> str:
    return _canonical_json(validate_tpt_effect_policy_contract(payload))


def tpt_effect_policy_hash(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(tpt_effect_policy_json(payload).encode("utf-8")).hexdigest()


def build_common_random_number_contract(
    *,
    contract_identifier: str = COMMON_RANDOM_NUMBER_PAIRED_BASELINE,
    review_status: str = "milestone_3c_common_random_number_review",
) -> dict[str, Any]:
    if contract_identifier != COMMON_RANDOM_NUMBER_PAIRED_BASELINE:
        raise ValueError("Unsupported common-random-number contract.")
    payload = {
        "contractVersion": COMMON_RANDOM_NUMBER_CONTRACT_VERSION,
        "contractIdentifier": contract_identifier,
        "sharedDraws": [
            "baseline age/sex/BCG/risk-factor population",
            "baseline active-TB assignment",
            "recent/remote/uninfected assignment",
            "infection-time history",
            "untreated active-TB progression potential",
            "mortality history when supplied to explicit helpers",
        ],
        "interventionSpecificDraws": [
            "screening selection and timing",
            "test outcome",
            "TPT initiation",
            "adverse-event discontinuation",
            "completion",
            "TPT efficacy",
        ],
        "caseAvertedDefinition": "direct active TB cases averted among the modelled population",
        "reviewStatus": str(review_status),
    }
    return validate_common_random_number_contract(payload)


def validate_common_random_number_contract(payload: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("commonRandomNumberContract must be a mapping.")
    if str(payload.get("contractVersion") or "") != COMMON_RANDOM_NUMBER_CONTRACT_VERSION:
        raise ValueError("Unsupported commonRandomNumberContract contractVersion.")
    if str(payload.get("contractIdentifier") or "") != COMMON_RANDOM_NUMBER_PAIRED_BASELINE:
        raise ValueError("Unsupported commonRandomNumberContract identifier.")
    return deepcopy(dict(payload))


def common_random_number_contract_json(payload: Mapping[str, Any]) -> str:
    return _canonical_json(validate_common_random_number_contract(payload))


def common_random_number_contract_hash(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(
        common_random_number_contract_json(payload).encode("utf-8")
    ).hexdigest()


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
            "Milestone 3C runner integration only supports the central progression curve."
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
        raise ValueError("Milestone 3C central runner policy requires riskFactorProgressionPolicy='none'.")
    mortality_policy = str(payload.get("mortalityPolicy") or "")
    if mortality_policy != COMPETING_MORTALITY_NOT_MODELLED:
        raise ValueError("Milestone 3C runner integration does not bundle mortality data.")
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
    tpt_policy = validate_tpt_effect_policy_contract(
        payload.get("tptEffectPolicy") or {}
    )
    if str(payload.get("tptEffectPolicyIdentifier") or "") != tpt_policy["policyIdentifier"]:
        raise ValueError("tptEffectPolicyIdentifier must match tptEffectPolicy.")
    if str(payload.get("tptEffectPolicyHash") or "") != tpt_effect_policy_hash(tpt_policy):
        raise ValueError("tptEffectPolicyHash does not match tptEffectPolicy.")
    crn_contract = validate_common_random_number_contract(
        payload.get("commonRandomNumberContract") or {}
    )
    if (
        str(payload.get("commonRandomNumberContractIdentifier") or "")
        != crn_contract["contractIdentifier"]
    ):
        raise ValueError(
            "commonRandomNumberContractIdentifier must match commonRandomNumberContract."
        )
    if str(payload.get("commonRandomNumberContractHash") or "") != common_random_number_contract_hash(
        crn_contract
    ):
        raise ValueError(
            "commonRandomNumberContractHash does not match commonRandomNumberContract."
        )
    active_policy = str(payload.get("activeTBObservationPolicy") or "")
    if active_policy != ACTIVE_TB_OBSERVATION_POLICY_VALIDATION_ONLY:
        raise ValueError("Active-TB observations must remain validation-only in Milestone 3C.")
    reinfection_policy = str(payload.get("reinfectionPolicy") or "")
    if reinfection_policy != curve["reinfectionPolicy"]:
        raise ValueError("reinfectionPolicy must match progressionCurve.")
    if reinfection_policy not in {
        REINFECTION_POLICY_RESET_CLOCK,
        REINFECTION_POLICY_NO_RESET_CLOCK,
    }:
        raise ValueError("Unsupported reinfectionPolicy.")
    validated = deepcopy(dict(payload))
    validated[ADJUSTMENT_FIELD] = adjustment_for_runner(payload)
    return validated


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
        ADJUSTMENT_FIELD: adjustment_for_runner(runner_config),
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
        progression_hazard_multiplier=adjustment_for_runner(runner_config)["multiplier"],
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
    progression_hazard_multiplier: float = 1.0,
) -> float:
    multiplier = positive_multiplier(progression_hazard_multiplier)
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
    ) - math.log1p(-u) / multiplier
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
        progression_hazard_multiplier=adjustment_for_runner(calibration["explicitRecentRemoteRunnerConfig"])["multiplier"],
    )
    result["activeTBObservationPolicy"] = ACTIVE_TB_OBSERVATION_POLICY_VALIDATION_ONLY
    result[ADJUSTMENT_FIELD] = adjustment_for_runner(calibration["explicitRecentRemoteRunnerConfig"])
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
            "isValid": False,
            "errors": [
                "The explicit recent/remote pathway does not produce an authoritative event ledger in Milestone 3C."
            ],
            "warnings": [
                "Event-ledger generation, economics and DALYs are intentionally out of scope for the explicit recent/remote pathway in Milestone 3C."
            ],
        },
    }


def build_explicit_epidemiology_replicate_summary(
    *,
    population: Mapping[str, np.ndarray],
    arrays: Mapping[str, np.ndarray],
    follow_horizon_years: float,
    runner_config: Mapping[str, Any],
    replicate_id: int | None = None,
    replicate_seed: int | None = None,
) -> dict[str, Any]:
    validated = validate_explicit_recent_remote_runner_config(runner_config)
    horizon = _finite_nonnegative_float(follow_horizon_years, "follow_horizon_years")
    infected = np.asarray(population["infected"], dtype=bool)
    recent = np.asarray(population["recentLTBIAtBaseline"], dtype=bool)
    remote = np.asarray(population["remoteLTBIAtBaseline"], dtype=bool)
    baseline_active = np.asarray(
        population.get("baselineActiveTBAtBaseline", np.zeros(len(infected), dtype=bool)),
        dtype=bool,
    )
    explicit_state = np.asarray(
        population.get("explicitTBIState", np.full(len(infected), STATE_UNINFECTED)),
        dtype=object,
    )
    prior_remote_plus_recent = np.asarray(
        population.get("priorRemotePlusRecent", np.zeros(len(infected), dtype=bool)),
        dtype=bool,
    )
    t_active = np.asarray(arrays["t_active"], dtype=float)
    screened = np.asarray(arrays["screened"], dtype=bool)
    screen_time = np.asarray(arrays["screen_time"], dtype=float)
    test_positive = np.asarray(arrays["test_positive"], dtype=bool)
    true_positive_latent = np.asarray(arrays["true_positive_latent"], dtype=bool)
    false_positive_test = np.asarray(arrays["false_positive_test"], dtype=bool)
    test_positive_active = np.asarray(arrays["test_positive_active"], dtype=bool)
    active_at_screen = np.asarray(arrays["active_at_screen"], dtype=bool)
    started_tpt = np.asarray(arrays["started_tpt"], dtype=bool)
    completed_tpt = np.asarray(arrays["completed_tpt"], dtype=bool)
    protected_full = np.asarray(arrays["protected_full"], dtype=bool)
    protected_partial = np.asarray(arrays["protected_partial"], dtype=bool)
    prevented_active_tb = np.asarray(arrays["prevented_active_tb"], dtype=bool)
    course_stop_time = np.asarray(arrays["course_stop_time"], dtype=float)
    active_by_horizon = infected & (t_active <= horizon)
    protected_any = protected_full | protected_partial
    intervention_active = active_by_horizon & ~prevented_active_tb
    before_screening = active_at_screen & active_by_horizon
    latent_at_screen = screened & infected & ~active_at_screen
    infected_screened = screened & infected
    uninfected_screened = screened & ~infected & ~baseline_active
    false_negative_latent = latent_at_screen & ~test_positive
    true_negative = uninfected_screened & ~test_positive
    tpt_eligible_after_assessment = true_positive_latent | false_positive_test
    true_positive_tpt_starts = started_tpt & true_positive_latent
    false_positive_tpt_starts = started_tpt & false_positive_test
    after_screen_before_effective_tpt = (
        started_tpt
        & infected
        & (t_active > screen_time)
        & (t_active <= course_stop_time)
        & active_by_horizon
    )
    subsequent_active = intervention_active & ~before_screening & ~after_screen_before_effective_tpt
    n = len(infected)
    uninfected = explicit_state == STATE_UNINFECTED
    base_counts = {
        "population": int(n),
        "baselineActiveTB": int(baseline_active.sum()),
        "eligiblePopulation": int((~baseline_active).sum()),
        "recentTBI": int(recent.sum()),
        "remoteOnlyTBI": int(remote.sum()),
        "uninfected": int(uninfected.sum()),
        "priorRemotePlusRecent": int(prior_remote_plus_recent.sum()),
    }
    comparator = {
        **base_counts,
        "totalActiveTBOverFollowUp": int(active_by_horizon.sum()),
        "directActiveTBCasesAvertedAmongModelledPopulation": 0,
        "annualActiveTB": _annual_counts_from_times(t_active[active_by_horizon], horizon),
    }
    intervention = {
        **base_counts,
        "screened": int(screened.sum()),
        "notScreened": int((~baseline_active).sum() - screened.sum()),
        "infectedScreened": int(infected_screened.sum()),
        "latentInfectedScreened": int(latent_at_screen.sum()),
        "falseNegativeLatent": int(false_negative_latent.sum()),
        "uninfectedScreened": int(uninfected_screened.sum()),
        "trueNegative": int(true_negative.sum()),
        "testPositive": int((screened & test_positive).sum()),
        "truePositiveLatent": int(true_positive_latent.sum()),
        "falsePositive": int(false_positive_test.sum()),
        "activeTBDetectedDuringWorkup": int(test_positive_active.sum()),
        "tptEligibleAfterTestingAndAssessment": int(
            tpt_eligible_after_assessment.sum()
        ),
        "tptStarts": int(started_tpt.sum()),
        "truePositiveTPTStarts": int(true_positive_tpt_starts.sum()),
        "falsePositiveTPTStarts": int(false_positive_tpt_starts.sum()),
        "tptCompletions": int(completed_tpt.sum()),
        "effectiveTPT": int(protected_any.sum()),
        "effectiveTPTFull": int(protected_full.sum()),
        "effectiveTPTPartial": int(protected_partial.sum()),
        "incidentActiveTBBeforeScreening": int(before_screening.sum()),
        "incidentActiveTBAfterScreenBeforeEffectiveTPT": int(
            after_screen_before_effective_tpt.sum()
        ),
        "subsequentIncidentActiveTB": int(subsequent_active.sum()),
        "totalActiveTBOverFollowUp": int(intervention_active.sum()),
        "directActiveTBCasesAvertedAmongModelledPopulation": int(
            prevented_active_tb.sum()
        ),
        "uninfectedFalsePositiveEpidemiologicalBenefit": 0,
        "annualActiveTB": _annual_counts_from_times(
            t_active[intervention_active],
            horizon,
        ),
        "annualDirectActiveTBCasesAverted": _annual_counts_from_times(
            t_active[prevented_active_tb],
            horizon,
        ),
    }
    summary = {
        "contractVersion": EPIDEMIOLOGICAL_INTERVENTION_OUTPUT_CONTRACT_VERSION,
        "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
        ADJUSTMENT_FIELD: adjustment_for_runner(validated),
        "replicateId": replicate_id,
        "replicateSeed": replicate_seed,
        "tptEffectPolicyIdentifier": validated["tptEffectPolicyIdentifier"],
        "tptEffectPolicyHash": validated["tptEffectPolicyHash"],
        "commonRandomNumberContractIdentifier": validated[
            "commonRandomNumberContractIdentifier"
        ],
        "commonRandomNumberContractHash": validated[
            "commonRandomNumberContractHash"
        ],
        "caseAvertedDefinition": "direct active TB cases averted among the modelled population",
        "comparator": comparator,
        "intervention": intervention,
        "byAgeGroup": _explicit_counts_by_group(
            population.get("ageGroup"),
            active_by_horizon,
            intervention_active,
            prevented_active_tb,
            recent,
            remote,
            baseline_active,
        ),
        "byEffectiveState": _explicit_counts_by_group(
            explicit_state,
            active_by_horizon,
            intervention_active,
            prevented_active_tb,
            recent,
            remote,
            baseline_active,
        ),
        "diagnosticMessages": [
            "Comparator and intervention share baseline states, infection clocks and untreated active-TB progression potential.",
            "Intervention-specific screening, treatment and efficacy draws are not shared with the comparator.",
            "No event-ledger costing, economics or DALYs are generated from this summary.",
        ],
    }
    summary["accountingValidation"] = validate_stochastic_intervention_accounting(
        summary
    )
    return summary


def summarise_explicit_epidemiology_replicates(
    replicate_summaries: list[Mapping[str, Any]],
    runner_config: Mapping[str, Any],
) -> dict[str, Any]:
    validated = validate_explicit_recent_remote_runner_config(runner_config)
    rows = []
    for item in replicate_summaries:
        base = {
            "replicateId": item.get("replicateId"),
            "replicateSeed": item.get("replicateSeed"),
        }
        for arm in ("comparator", "intervention"):
            for key, value in (item.get(arm) or {}).items():
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    rows.append({**base, "arm": arm, "metric": key, "value": float(value)})
    metric_summary = []
    for arm in ("comparator", "intervention"):
        arm_rows = [row for row in rows if row["arm"] == arm]
        for metric in sorted({row["metric"] for row in arm_rows}):
            values = np.asarray(
                [row["value"] for row in arm_rows if row["metric"] == metric],
                dtype=float,
            )
            if len(values) == 0:
                continue
            metric_summary.append(
                {
                    "arm": arm,
                    "metric": metric,
                    "mean": float(np.mean(values)),
                    "low95": float(np.quantile(values, 0.025)),
                    "high95": float(np.quantile(values, 0.975)),
                }
            )
    return {
        "contractVersion": EPIDEMIOLOGICAL_INTERVENTION_OUTPUT_CONTRACT_VERSION,
        "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
        "available": True,
        "replicateCount": len(replicate_summaries),
        ADJUSTMENT_FIELD: adjustment_for_runner(validated),
        "tptEffectPolicyIdentifier": validated["tptEffectPolicyIdentifier"],
        "tptEffectPolicyHash": validated["tptEffectPolicyHash"],
        "commonRandomNumberContractIdentifier": validated[
            "commonRandomNumberContractIdentifier"
        ],
        "commonRandomNumberContractHash": validated[
            "commonRandomNumberContractHash"
        ],
        "caseAvertedDefinition": "direct active TB cases averted among the modelled population",
        "metricSummary": metric_summary,
        "replicates": [dict(item) for item in replicate_summaries],
        "diagnosticMessages": [
            "This is a non-economic epidemiological summary; it is not an event ledger and is not costable.",
        ],
    }


def validate_stochastic_intervention_accounting(
    summary: Mapping[str, Any],
    *,
    tolerance: float = 1e-9,
) -> dict[str, Any]:
    comparator = summary.get("comparator") or {}
    intervention = summary.get("intervention") or {}
    by_effective_state = summary.get("byEffectiveState") or []
    counts = {
        "population": comparator.get("population", intervention.get("population", 0)),
        "baselineActiveTB": comparator.get(
            "baselineActiveTB",
            intervention.get("baselineActiveTB", 0),
        ),
        "eligiblePopulation": comparator.get(
            "eligiblePopulation",
            intervention.get("eligiblePopulation", 0),
        ),
        "recentTBI": comparator.get("recentTBI", intervention.get("recentTBI", 0)),
        "remoteOnlyTBI": comparator.get(
            "remoteOnlyTBI",
            intervention.get("remoteOnlyTBI", 0),
        ),
        "uninfected": comparator.get("uninfected", intervention.get("uninfected", 0)),
        "priorRemotePlusRecent": comparator.get(
            "priorRemotePlusRecent",
            intervention.get("priorRemotePlusRecent", 0),
        ),
        "screened": intervention.get("screened", 0),
        "notScreened": intervention.get("notScreened", 0),
        "infectedScreened": intervention.get("infectedScreened", 0),
        "latentInfectedScreened": intervention.get("latentInfectedScreened", 0),
        "falseNegativeLatent": intervention.get("falseNegativeLatent", 0),
        "uninfectedScreened": intervention.get("uninfectedScreened", 0),
        "trueNegative": intervention.get("trueNegative", 0),
        "truePositiveLatent": intervention.get("truePositiveLatent", 0),
        "falsePositive": intervention.get("falsePositive", 0),
        "tptEligibleAfterTestingAndAssessment": intervention.get(
            "tptEligibleAfterTestingAndAssessment",
            0,
        ),
        "tptStarts": intervention.get("tptStarts", 0),
        "truePositiveTPTStarts": intervention.get("truePositiveTPTStarts", 0),
        "falsePositiveTPTStarts": intervention.get("falsePositiveTPTStarts", 0),
        "tptCompletions": intervention.get("tptCompletions", 0),
        "effectiveTPT": intervention.get("effectiveTPT", 0),
        "comparatorActiveTB": comparator.get("totalActiveTBOverFollowUp", 0),
        "interventionActiveTB": intervention.get("totalActiveTBOverFollowUp", 0),
        "directActiveTBCasesAverted": intervention.get(
            "directActiveTBCasesAvertedAmongModelledPopulation",
            0,
        ),
        "recentDirectActiveTBCasesAverted": _group_metric(
            by_effective_state,
            STATE_RECENT,
        ),
        "remoteOnlyDirectActiveTBCasesAverted": _group_metric(
            by_effective_state,
            STATE_REMOTE_ONLY,
        ),
        "incidentActiveTBBeforeScreening": intervention.get(
            "incidentActiveTBBeforeScreening",
            0,
        ),
        "incidentActiveTBAfterScreenBeforeEffectiveTPT": intervention.get(
            "incidentActiveTBAfterScreenBeforeEffectiveTPT",
            0,
        ),
        "subsequentIncidentActiveTB": intervention.get("subsequentIncidentActiveTB", 0),
    }
    errors: list[str] = []
    _append_counting_identity_errors(
        errors,
        counts=counts,
        by_effective_state=by_effective_state,
        tolerance=tolerance,
    )
    return {
        "isValid": not errors,
        "errors": errors,
        "tolerance": float(tolerance),
    }


def validate_deterministic_intervention_accounting(
    deterministic: Mapping[str, Any],
    *,
    tolerance: float = 1e-9,
) -> dict[str, Any]:
    state = deterministic.get("stateProportions") or {}
    by_horizon = deterministic.get("byHorizon") or {}
    errors: list[str] = []
    state_total = (
        float(state.get(STATE_BASELINE_ACTIVE_TB, 0.0))
        + float(state.get(STATE_RECENT, 0.0))
        + float(state.get(STATE_REMOTE_ONLY, 0.0))
        + float(state.get(STATE_UNINFECTED, 0.0))
    )
    _check_close(errors, state_total, 1.0, tolerance, "baseline state proportions")
    for horizon_key, horizon in by_horizon.items():
        prefix = f"horizon {horizon_key}: "
        counts = {
            "population": 1.0,
            "baselineActiveTB": float(state.get(STATE_BASELINE_ACTIVE_TB, 0.0)),
            "eligiblePopulation": 1.0
            - float(state.get(STATE_BASELINE_ACTIVE_TB, 0.0)),
            "recentTBI": float(state.get(STATE_RECENT, 0.0)),
            "remoteOnlyTBI": float(state.get(STATE_REMOTE_ONLY, 0.0)),
            "uninfected": float(state.get(STATE_UNINFECTED, 0.0)),
            "priorRemotePlusRecent": float(state.get("priorRemotePlusRecent", 0.0)),
            "screened": float(horizon.get("screenedProportion", 0.0)),
            "notScreened": float(horizon.get("notScreenedProportion", 0.0)),
            "infectedScreened": float(horizon.get("infectedScreenedProportion", 0.0)),
            "latentInfectedScreened": float(
                horizon.get("latentInfectedScreenedProportion", 0.0)
            ),
            "falseNegativeLatent": float(
                horizon.get("falseNegativeLatentProportion", 0.0)
            ),
            "uninfectedScreened": float(
                horizon.get("uninfectedScreenedProportion", 0.0)
            ),
            "trueNegative": float(horizon.get("trueNegativeProportion", 0.0)),
            "truePositiveLatent": float(
                horizon.get("truePositiveLatentProportion", 0.0)
            ),
            "falsePositive": float(horizon.get("falsePositiveProportion", 0.0)),
            "tptEligibleAfterTestingAndAssessment": float(
                horizon.get("tptEligibleAfterTestingAndAssessmentProportion", 0.0)
            ),
            "tptStarts": float(horizon.get("tptStartsProportion", 0.0)),
            "truePositiveTPTStarts": float(
                horizon.get("truePositiveTPTStartsProportion", 0.0)
            ),
            "falsePositiveTPTStarts": float(
                horizon.get("falsePositiveTPTStartsProportion", 0.0)
            ),
            "tptCompletions": float(horizon.get("tptCompletionsProportion", 0.0)),
            "effectiveTPT": float(horizon.get("effectiveTPTProportion", 0.0)),
            "comparatorActiveTB": float(
                horizon.get("comparatorActiveTBProportion", 0.0)
            ),
            "interventionActiveTB": float(
                horizon.get("interventionActiveTBProportion", 0.0)
            ),
            "directActiveTBCasesAverted": float(
                horizon.get("directActiveTBCasesAvertedProportion", 0.0)
            ),
            "recentDirectActiveTBCasesAverted": float(
                horizon.get("recentDirectActiveTBCasesAvertedProportion", 0.0)
            ),
            "remoteOnlyDirectActiveTBCasesAverted": float(
                horizon.get("remoteOnlyDirectActiveTBCasesAvertedProportion", 0.0)
            ),
            "incidentActiveTBBeforeScreening": float(
                horizon.get("incidentActiveTBBeforeScreeningProportion", 0.0)
            ),
            "incidentActiveTBAfterScreenBeforeEffectiveTPT": float(
                horizon.get(
                    "incidentActiveTBAfterScreenBeforeEffectiveTPTProportion",
                    0.0,
                )
            ),
            "subsequentIncidentActiveTB": float(
                horizon.get("subsequentIncidentActiveTBProportion", 0.0)
            ),
        }
        horizon_errors: list[str] = []
        _append_counting_identity_errors(
            horizon_errors,
            counts=counts,
            by_effective_state=None,
            tolerance=tolerance,
        )
        errors.extend(prefix + error for error in horizon_errors)
    return {
        "isValid": not errors,
        "errors": errors,
        "tolerance": float(tolerance),
    }


def _append_counting_identity_errors(
    errors: list[str],
    *,
    counts: Mapping[str, Any],
    by_effective_state: list[Mapping[str, Any]] | None,
    tolerance: float,
) -> None:
    population = float(counts.get("population", 0.0))
    baseline = float(counts.get("baselineActiveTB", 0.0))
    recent = float(counts.get("recentTBI", 0.0))
    remote = float(counts.get("remoteOnlyTBI", 0.0))
    uninfected = float(counts.get("uninfected", 0.0))
    eligible = float(counts.get("eligiblePopulation", population - baseline))
    screened = float(counts.get("screened", 0.0))
    not_screened = float(counts.get("notScreened", eligible - screened))
    true_positive = float(counts.get("truePositiveLatent", 0.0))
    false_negative = float(counts.get("falseNegativeLatent", 0.0))
    latent_screened = float(counts.get("latentInfectedScreened", 0.0))
    false_positive = float(counts.get("falsePositive", 0.0))
    true_negative = float(counts.get("trueNegative", 0.0))
    uninfected_screened = float(counts.get("uninfectedScreened", 0.0))
    tpt_eligible = float(
        counts.get(
            "tptEligibleAfterTestingAndAssessment",
            true_positive + false_positive,
        )
    )
    tpt_starts = float(counts.get("tptStarts", 0.0))
    tpt_completions = float(counts.get("tptCompletions", 0.0))
    effective_tpt = float(counts.get("effectiveTPT", 0.0))
    true_positive_tpt_starts = float(
        counts.get("truePositiveTPTStarts", tpt_starts)
    )
    comparator_tb = float(counts.get("comparatorActiveTB", 0.0))
    intervention_tb = float(counts.get("interventionActiveTB", 0.0))
    averted = float(counts.get("directActiveTBCasesAverted", 0.0))
    recent_averted = float(counts.get("recentDirectActiveTBCasesAverted", 0.0))
    remote_averted = float(counts.get("remoteOnlyDirectActiveTBCasesAverted", 0.0))
    before_screen = float(counts.get("incidentActiveTBBeforeScreening", 0.0))
    before_effective = float(
        counts.get("incidentActiveTBAfterScreenBeforeEffectiveTPT", 0.0)
    )
    subsequent = float(counts.get("subsequentIncidentActiveTB", 0.0))

    _check_close(
        errors,
        baseline + recent + remote + uninfected,
        population,
        tolerance,
        "baseline active TB + recent + remote-only + uninfected",
    )
    _check_close(
        errors,
        population - baseline,
        eligible,
        tolerance,
        "eligible population after baseline active TB exclusion",
    )
    _check_close(
        errors,
        screened + not_screened,
        eligible,
        tolerance,
        "screened + not screened",
    )
    _check_close(
        errors,
        true_positive + false_negative,
        latent_screened,
        tolerance,
        "true-positive latent + false-negative latent",
    )
    _check_close(
        errors,
        false_positive + true_negative,
        uninfected_screened,
        tolerance,
        "false-positive + true-negative",
    )
    _check_close(
        errors,
        true_positive + false_positive,
        tpt_eligible,
        tolerance,
        "TPT-eligible after testing and clinical assessment",
    )
    _check_close(
        errors,
        comparator_tb - intervention_tb,
        averted,
        tolerance,
        "cases averted",
    )
    _check_close(
        errors,
        before_screen + before_effective + subsequent,
        intervention_tb,
        tolerance,
        "timing-specific active TB components",
    )
    _check_close(
        errors,
        recent_averted + remote_averted,
        averted,
        tolerance,
        "recent plus remote-only cases averted",
    )
    if by_effective_state is not None:
        grouped_averted = sum(
            float(row.get("directActiveTBCasesAverted", 0.0))
            for row in by_effective_state
        )
        _check_close(
            errors,
            grouped_averted,
            averted,
            tolerance,
            "effective-state cases averted",
        )
    _check_less_equal(errors, screened, eligible, tolerance, "screened")
    _check_less_equal(
        errors,
        true_positive,
        latent_screened,
        tolerance,
        "true-positive latent",
    )
    _check_less_equal(
        errors,
        tpt_starts,
        tpt_eligible,
        tolerance,
        "TPT starts",
    )
    _check_less_equal(
        errors,
        tpt_completions,
        tpt_starts,
        tolerance,
        "TPT completions",
    )
    _check_less_equal(
        errors,
        effective_tpt,
        true_positive_tpt_starts,
        tolerance,
        "effective TPT",
    )
    _check_less_equal(
        errors,
        averted,
        comparator_tb,
        tolerance,
        "cases averted",
    )
    if intervention_tb < -tolerance:
        errors.append("intervention active TB cannot be negative.")


def _check_close(
    errors: list[str],
    observed: float,
    expected: float,
    tolerance: float,
    label: str,
) -> None:
    if abs(float(observed) - float(expected)) > float(tolerance):
        errors.append(
            f"{label} identity failed: observed {observed:.12g}, expected {expected:.12g}."
        )


def _check_less_equal(
    errors: list[str],
    observed: float,
    upper: float,
    tolerance: float,
    label: str,
) -> None:
    if float(observed) - float(upper) > float(tolerance):
        errors.append(
            f"{label} exceeded its bound: observed {observed:.12g}, upper {upper:.12g}."
        )


def deterministic_explicit_recent_remote_intervention_expectations(
    calibration: Mapping[str, Any],
    *,
    config: Mapping[str, Any],
    regimen: Mapping[str, Any],
    horizons_years: tuple[float, ...] = (20.0,),
    infection_quadrature_points: int = 32,
    screen_quadrature_points: int = 24,
    competing_survival_curve: Any = None,
    competing_mortality_hazard: Any = None,
    integration_steps: int = 512,
) -> dict[str, Any]:
    if not is_explicit_recent_remote_calibration(calibration):
        raise ValueError(
            "deterministic intervention expectations require explicit recent/remote calibration."
        )
    if competing_survival_curve is not None and competing_mortality_hazard is not None:
        raise ValueError(
            "Provide either competing_survival_curve or competing_mortality_hazard, not both."
        )
    runner_config = validate_explicit_recent_remote_runner_config(
        calibration["explicitRecentRemoteRunnerConfig"]
    )
    tpt_policy = validate_tpt_effect_policy_contract(
        runner_config["tptEffectPolicy"]
    )
    curve = validate_progression_curve_contract(calibration["progressionCurve"])
    risk_adjustment = adjustment_for_runner(runner_config)
    horizons = tuple(
        _finite_nonnegative_float(value, "horizons_years")
        for value in horizons_years
    )
    q_infection = _positive_int(
        infection_quadrature_points,
        "infection_quadrature_points",
    )
    q_screen = _positive_int(screen_quadrature_points, "screen_quadrature_points")
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
    screen_coverage = _probability(config.get("screenCoverage", 0.0), "screenCoverage")
    screen_window = _positive_float(
        config.get("screeningWindowYears", config.get("screenWindow", 0.0)),
        "screeningWindowYears",
    )
    p_start = _probability(_default_value(config.get("pStartTPT"), 0.85), "pStartTPT")
    test_sensitivity = _deterministic_test_sensitivity(config)
    false_positive_rate = 1.0 - _deterministic_test_specificity(config)
    p_adr = _probability(regimen.get("pADRstop", 0.0), "regimen.pADRstop")
    p_complete = _probability(regimen.get("pComplete", 0.0), "regimen.pComplete")
    if p_adr + p_complete > 1.0 + 1e-12:
        raise ValueError("regimen pADRstop plus pComplete cannot exceed one.")
    p_stop_other = max(0.0, 1.0 - p_adr - p_complete)
    months = _finite_nonnegative_float(regimen.get("months", 0.0), "regimen.months")
    eff_full = _probability(regimen.get("effFull", 0.0), "regimen.effFull")
    adr_fraction = _probability(
        _default_value(config.get("partialDoseFractionADR"), 0.30),
        "partialDoseFractionADR",
    )
    other_fraction = _probability(
        _default_value(config.get("partialDoseFractionOther"), 0.60),
        "partialDoseFractionOther",
    )
    treatment_paths = (
        {
            "path": "completed",
            "probability": p_complete,
            "doseFraction": 1.0,
            "efficacy": eff_full,
        },
        {
            "path": "adr_stop",
            "probability": p_adr,
            "doseFraction": adr_fraction,
            "efficacy": float(regimen_partial_efficacy(regimen, adr_fraction)),
        },
        {
            "path": "stopped_other",
            "probability": p_stop_other,
            "doseFraction": other_fraction,
            "efficacy": float(regimen_partial_efficacy(regimen, other_fraction)),
        },
    )
    p_recent = np.asarray(probabilities["recent"], dtype=float)
    p_remote = np.asarray(probabilities["remote_only"], dtype=float)
    p_uninfected = np.asarray(probabilities["uninfected"], dtype=float)
    p_recent_with_prior = np.asarray(
        probabilities["recent_with_prior_remote"],
        dtype=float,
    )
    totals_by_horizon: dict[str, dict[str, Any]] = {}
    for horizon in horizons:
        comparator_active = 0.0
        comparator_recent = 0.0
        comparator_remote = 0.0
        prevented = 0.0
        prevented_recent = 0.0
        prevented_remote = 0.0
        active_before_screen = 0.0
        active_after_screen_before_effective = 0.0
        subsequent_active = 0.0
        true_positive_latent = 0.0
        latent_infected_screened = 0.0
        infected_screened = 0.0
        active_detected_workup = 0.0
        for idx, age in enumerate(ages):
            state_components = (
                (
                    STATE_RECENT,
                    float(p_recent[idx]),
                    0.0,
                    float(durations.recent_years[idx]),
                    recent_hazard,
                ),
                (
                    STATE_REMOTE_ONLY,
                    float(p_remote[idx]),
                    recent_window,
                    float(durations.remote_years[idx]),
                    remote_hazard,
                ),
            )
            if runner_config["reinfectionPolicy"] == REINFECTION_POLICY_NO_RESET_CLOCK:
                # Recently reinfected people with remote history retain that
                # remote clock; they remain in the recent prevalence group.
                prior = float(p_recent_with_prior[idx])
                state_components = (
                    (STATE_RECENT, float(p_recent[idx]) - prior, 0.0, float(durations.recent_years[idx]), recent_hazard),
                    (STATE_RECENT, prior, recent_window, float(durations.remote_years[idx]), remote_hazard),
                    state_components[1],
                )
            age_weight = at_risk_fraction * float(weights[idx])
            for state_name, state_probability, window_start, duration, hazard in state_components:
                if state_probability <= 0.0:
                    continue
                moments = _infected_intervention_moments(
                    curve,
                    age=float(age),
                    horizon=horizon,
                    window_start_years_before_baseline=window_start,
                    window_duration_years=duration,
                    infection_hazard_per_year=hazard,
                    screen_window_years=screen_window,
                    treatment_months=months,
                    treatment_paths=treatment_paths,
                    infection_quadrature_points=q_infection,
                    screen_quadrature_points=q_screen,
                    competing_survival_curve=competing_survival_curve,
                    competing_mortality_hazard=competing_mortality_hazard,
                    integration_steps=integration_steps,
                    progression_hazard_multiplier=risk_adjustment["multiplier"],
                )
                scale = age_weight * state_probability
                comparator_component = scale * moments["activeByHorizon"]
                prevented_component = (
                    scale
                    * screen_coverage
                    * test_sensitivity
                    * p_start
                    * moments["preventedByEffectiveTPT"]
                )
                comparator_active += comparator_component
                prevented += prevented_component
                if state_name == STATE_RECENT:
                    comparator_recent += comparator_component
                    prevented_recent += prevented_component
                else:
                    comparator_remote += comparator_component
                    prevented_remote += prevented_component
                active_before_screen += (
                    scale * screen_coverage * moments["activeBeforeScreening"]
                )
                latent_infected_screened += (
                    scale * screen_coverage * moments["latentAtScreen"]
                )
                infected_screened += (
                    scale
                    * screen_coverage
                    * (
                        moments["latentAtScreen"]
                        + moments["activeBeforeScreening"]
                    )
                )
                active_after_screen_before_effective += (
                    scale
                    * screen_coverage
                    * test_sensitivity
                    * p_start
                    * moments["activeAfterScreenBeforeEffectiveTPT"]
                )
                true_positive_latent += (
                    scale
                    * screen_coverage
                    * test_sensitivity
                    * moments["latentAtScreen"]
                )
                active_detected_workup += (
                    scale
                    * screen_coverage
                    * test_sensitivity
                    * moments["activeBeforeScreening"]
                )
        intervention_active = max(0.0, comparator_active - prevented)
        subsequent_active = max(
            0.0,
            intervention_active
            - active_before_screen
            - active_after_screen_before_effective,
        )
        uninfected_prop = at_risk_fraction * float(np.sum(weights * p_uninfected))
        screened_prop = at_risk_fraction * screen_coverage
        uninfected_screened = uninfected_prop * screen_coverage
        false_negative_latent = latent_infected_screened * (1.0 - test_sensitivity)
        false_positive = uninfected_prop * screen_coverage * false_positive_rate
        true_negative = uninfected_screened * (1.0 - false_positive_rate)
        tpt_eligible = true_positive_latent + false_positive
        tpt_starts = p_start * (true_positive_latent + false_positive)
        true_positive_tpt_starts = p_start * true_positive_latent
        false_positive_tpt_starts = p_start * false_positive
        tpt_completions = p_complete * tpt_starts
        effective_tpt = p_start * true_positive_latent * sum(
            path["probability"] * path["efficacy"] for path in treatment_paths
        )
        totals_by_horizon[f"{horizon:g}"] = {
            "horizonYears": horizon,
            "comparatorActiveTBProportion": comparator_active,
            "interventionActiveTBProportion": intervention_active,
            "directActiveTBCasesAvertedProportion": prevented,
            "recentComparatorActiveTBProportion": comparator_recent,
            "remoteOnlyComparatorActiveTBProportion": comparator_remote,
            "recentDirectActiveTBCasesAvertedProportion": prevented_recent,
            "remoteOnlyDirectActiveTBCasesAvertedProportion": prevented_remote,
            "screenedProportion": screened_prop,
            "notScreenedProportion": at_risk_fraction - screened_prop,
            "infectedScreenedProportion": infected_screened,
            "latentInfectedScreenedProportion": latent_infected_screened,
            "falseNegativeLatentProportion": false_negative_latent,
            "uninfectedScreenedProportion": uninfected_screened,
            "trueNegativeProportion": true_negative,
            "testPositiveProportion": true_positive_latent
            + active_detected_workup
            + false_positive,
            "truePositiveLatentProportion": true_positive_latent,
            "falsePositiveProportion": false_positive,
            "activeTBDetectedDuringWorkupProportion": active_detected_workup,
            "tptEligibleAfterTestingAndAssessmentProportion": tpt_eligible,
            "tptStartsProportion": tpt_starts,
            "truePositiveTPTStartsProportion": true_positive_tpt_starts,
            "falsePositiveTPTStartsProportion": false_positive_tpt_starts,
            "tptCompletionsProportion": tpt_completions,
            "effectiveTPTProportion": effective_tpt,
            "incidentActiveTBBeforeScreeningProportion": active_before_screen,
            "incidentActiveTBAfterScreenBeforeEffectiveTPTProportion": (
                active_after_screen_before_effective
            ),
            "subsequentIncidentActiveTBProportion": subsequent_active,
        }
    result = {
        "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
        "method": "deterministic_epidemiological_intervention_expected_value",
        ADJUSTMENT_FIELD: risk_adjustment,
        "integrationMethod": DETERMINISTIC_INTERVENTION_INTEGRATION_METHOD,
        "infectionQuadraturePoints": q_infection,
        "screenQuadraturePoints": q_screen,
        "usesMeanTimeSinceInfection": False,
        "tptEffectPolicyIdentifier": tpt_policy["policyIdentifier"],
        "caseAvertedDefinition": "direct active TB cases averted among the modelled population",
        "riskFactorProgressionPolicy": runner_config["riskFactorProgressionPolicy"],
        "mortalityPolicy": (
            COMPETING_MORTALITY_CONSTANT_HAZARD
            if competing_mortality_hazard is not None
            else (
                "external_survival_curve"
                if competing_survival_curve is not None
                else COMPETING_MORTALITY_NOT_MODELLED
            )
        ),
        "activeTBObservationPolicy": runner_config["activeTBObservationPolicy"],
        "baselineStateSequencing": [
            STATE_BASELINE_ACTIVE_TB,
            STATE_RECENT,
            STATE_REMOTE_ONLY,
            STATE_UNINFECTED,
        ],
        "stateProportions": {
            STATE_BASELINE_ACTIVE_TB: baseline_active,
            STATE_RECENT: at_risk_fraction * float(np.sum(weights * p_recent)),
            STATE_REMOTE_ONLY: at_risk_fraction * float(np.sum(weights * p_remote)),
            STATE_UNINFECTED: at_risk_fraction * float(np.sum(weights * p_uninfected)),
            "priorRemotePlusRecent": at_risk_fraction
            * float(np.sum(weights * p_recent_with_prior)),
        },
        "byHorizon": totals_by_horizon,
        "diagnosticMessages": [
            "Deterministic intervention expectations integrate over infection time and screening time.",
            "Effective TPT removes only future progression after the treatment stop time.",
            "Uninfected false positives receive no epidemiological TB-prevention benefit.",
            "No event-ledger costing, economics, DALYs or report-facing outputs are produced.",
        ],
    }
    result["accountingValidation"] = validate_deterministic_intervention_accounting(
        result
    )
    return result


def stochastic_explicit_recent_remote_intervention_summary(
    calibration: Mapping[str, Any],
    *,
    config: Mapping[str, Any],
    regimen: Mapping[str, Any],
    population_size: int | None = None,
    seed: int | None = None,
    rng: Any = None,
    follow_horizon_years: float | None = None,
    competing_mortality_hazard: Any = None,
    competing_survival_curve: Any = None,
) -> dict[str, Any]:
    """Isolated stochastic intervention summary used for 3D reconciliation tests."""
    if rng is None and seed is None:
        raise ValueError("A seed or explicit rng is required.")
    if rng is not None and seed is not None:
        raise ValueError("Provide either seed or rng, not both.")
    if competing_mortality_hazard is not None and competing_survival_curve is not None:
        raise ValueError(
            "Provide either competing_mortality_hazard or competing_survival_curve, not both."
        )
    if not is_explicit_recent_remote_calibration(calibration):
        raise ValueError(
            "stochastic intervention summaries require explicit recent/remote calibration."
        )
    generator = rng if rng is not None else np.random.default_rng(int(seed))
    n = _positive_int(
        population_size if population_size is not None else config.get("N", 0),
        "population_size",
    )
    horizon = _finite_nonnegative_float(
        follow_horizon_years
        if follow_horizon_years is not None
        else config.get("followHorizon", config.get("followUpHorizonYears", 20.0)),
        "follow_horizon_years",
    )
    pars = calibration["parameters"]
    age_values = np.asarray(pars["exactAgeValues"], dtype=float)
    age_prob = np.asarray(pars["exactAgeProb"], dtype=float)
    age_prob = age_prob / age_prob.sum()
    ages = generator.choice(age_values, size=n, p=age_prob)
    risk_prev = config.get("riskPrev") if isinstance(config.get("riskPrev"), Mapping) else {}
    bcg_probability = _probability(
        _default_value(risk_prev.get("BCG", 0.0), 0.0),
        "riskPrev.BCG",
    )
    population = {
        "ageYears": ages,
        "ageGroup": np.asarray([_age_group_label(float(age)) for age in ages], dtype=object),
        "BCG": generator.random(n) < bcg_probability if bcg_probability > 0.0 else np.zeros(n, dtype=bool),
    }
    population, recent, remote, _transition_time, active_time = apply_explicit_recent_remote_population_assignment(
        population,
        calibration,
        generator,
    )
    population["recentLTBIAtBaseline"] = recent
    population["remoteLTBIAtBaseline"] = remote
    baseline_active = np.asarray(population["baselineActiveTBAtBaseline"], dtype=bool)
    infected = np.asarray(population["infected"], dtype=bool)
    death_time = np.full(n, np.inf, dtype=float)
    if competing_mortality_hazard is not None:
        death_time = _draw_death_times_from_hazard(
            competing_mortality_hazard,
            ages,
            generator,
            horizon,
        )
    elif competing_survival_curve is not None:
        death_time = _draw_death_times_from_survival_curve(
            competing_survival_curve,
            ages,
            generator,
            horizon,
        )
    observed_active_time = np.where(active_time < death_time, active_time, np.inf)
    screen_coverage = _probability(config.get("screenCoverage", 0.0), "screenCoverage")
    screen_window = _positive_float(
        config.get("screeningWindowYears", config.get("screenWindow", 0.0)),
        "screeningWindowYears",
    )
    screened = (
        (~baseline_active)
        & (generator.random(n) < screen_coverage)
    )
    screen_time = np.full(n, np.inf, dtype=float)
    if screened.any():
        screen_time[screened] = generator.random(int(screened.sum())) * screen_window
    screened = screened & (screen_time < death_time)
    screen_time[~screened] = np.inf
    active_at_screen = screened & infected & (observed_active_time <= screen_time)
    latent_at_screen = screened & infected & (observed_active_time > screen_time)
    sensitivity = _deterministic_test_sensitivity(config)
    specificity = _deterministic_test_specificity(config)
    test_positive = np.zeros(n, dtype=bool)
    inf_screened = screened & infected
    uninf_screened = screened & ~infected
    if inf_screened.any():
        test_positive[inf_screened] = (
            generator.random(int(inf_screened.sum())) < sensitivity
        )
    if uninf_screened.any():
        test_positive[uninf_screened] = (
            generator.random(int(uninf_screened.sum())) < (1.0 - specificity)
        )
    true_positive_latent = latent_at_screen & test_positive
    false_positive_test = uninf_screened & test_positive
    test_positive_active = active_at_screen & test_positive
    p_start = _probability(_default_value(config.get("pStartTPT"), 0.85), "pStartTPT")
    eligible_tpt = screened & test_positive & ~active_at_screen & ~baseline_active
    started_tpt = eligible_tpt & (generator.random(n) < p_start)
    completed_tpt = np.zeros(n, dtype=bool)
    stopped_other = np.zeros(n, dtype=bool)
    adr_stop = np.zeros(n, dtype=bool)
    course_stop_time = np.full(n, np.inf, dtype=float)
    full_eff_assigned = np.zeros(n, dtype=float)
    partial_eff_assigned = np.zeros(n, dtype=float)
    treat_idx = np.flatnonzero(started_tpt)
    if len(treat_idx) > 0:
        p_adr = _probability(regimen.get("pADRstop", 0.0), "regimen.pADRstop")
        p_complete = _probability(regimen.get("pComplete", 0.0), "regimen.pComplete")
        adr_draw = generator.random(len(treat_idx)) < p_adr
        adr_stop[treat_idx[adr_draw]] = True
        remain_idx = treat_idx[~adr_draw]
        p_complete_given_no_adr = 0.0
        if (1.0 - p_adr) > 0.0:
            p_complete_given_no_adr = min(max(p_complete / (1.0 - p_adr), 0.0), 1.0)
        if len(remain_idx) > 0:
            comp_draw = generator.random(len(remain_idx)) < p_complete_given_no_adr
            completed_tpt[remain_idx[comp_draw]] = True
            stopped_other[remain_idx[~comp_draw]] = True
        dose_fraction = np.zeros(n, dtype=float)
        dose_fraction[completed_tpt] = 1.0
        dose_fraction[adr_stop] = _probability(
            _default_value(config.get("partialDoseFractionADR"), 0.30),
            "partialDoseFractionADR",
        )
        dose_fraction[stopped_other] = _probability(
            _default_value(config.get("partialDoseFractionOther"), 0.60),
            "partialDoseFractionOther",
        )
        course_stop_time[started_tpt] = (
            screen_time[started_tpt]
            + float(regimen.get("months", 0.0)) * dose_fraction[started_tpt] / 12.0
        )
        full_eff_assigned[completed_tpt] = _probability(
            regimen.get("effFull", 0.0),
            "regimen.effFull",
        )
        incomplete = started_tpt & ~completed_tpt
        if incomplete.any():
            partial_eff_assigned[incomplete] = regimen_partial_efficacy(
                regimen,
                dose_fraction[incomplete].tolist(),
            )
    comparator_case = infected & (observed_active_time <= horizon)
    eligible_for_future_tpt_protection = (
        infected
        & (course_stop_time < death_time)
        & (observed_active_time > course_stop_time)
    )
    protected_full = (
        completed_tpt
        & eligible_for_future_tpt_protection
        & (generator.random(n) < full_eff_assigned)
    )
    protected_partial = (
        (started_tpt & ~completed_tpt)
        & eligible_for_future_tpt_protection
        & (generator.random(n) < partial_eff_assigned)
    )
    prevented_active_tb = (protected_full | protected_partial) & comparator_case
    return build_explicit_epidemiology_replicate_summary(
        population=population,
        arrays={
            "t_active": observed_active_time,
            "screened": screened,
            "screen_time": screen_time,
            "test_positive": test_positive,
            "true_positive_latent": true_positive_latent,
            "false_positive_test": false_positive_test,
            "test_positive_active": test_positive_active,
            "active_at_screen": active_at_screen,
            "started_tpt": started_tpt,
            "completed_tpt": completed_tpt,
            "protected_full": protected_full,
            "protected_partial": protected_partial,
            "prevented_active_tb": prevented_active_tb,
            "course_stop_time": course_stop_time,
        },
        follow_horizon_years=horizon,
        runner_config=calibration["explicitRecentRemoteRunnerConfig"],
        replicate_seed=None if seed is None else int(seed),
    )


def reconcile_deterministic_stochastic_explicit_intervention(
    calibration: Mapping[str, Any],
    *,
    config: Mapping[str, Any],
    regimen: Mapping[str, Any],
    population_size: int = 50000,
    seeds: tuple[int, ...] = (20261013, 20261014, 20261015),
    horizon_years: float | None = None,
    z_score: float = 4.0,
    absolute_floor_count: float = 6.0,
    competing_mortality_hazard: Any = None,
    competing_survival_curve: Any = None,
    infection_quadrature_points: int = 32,
    screen_quadrature_points: int = 24,
    integration_steps: int = 512,
) -> dict[str, Any]:
    seed_values = tuple(int(value) for value in seeds)
    if not seed_values:
        raise ValueError("At least one stochastic seed is required.")
    horizon = _finite_nonnegative_float(
        horizon_years
        if horizon_years is not None
        else config.get("followHorizon", config.get("followUpHorizonYears", 20.0)),
        "horizon_years",
    )
    deterministic = deterministic_explicit_recent_remote_intervention_expectations(
        calibration,
        config=config,
        regimen=regimen,
        horizons_years=(horizon,),
        competing_mortality_hazard=competing_mortality_hazard,
        competing_survival_curve=competing_survival_curve,
        infection_quadrature_points=infection_quadrature_points,
        screen_quadrature_points=screen_quadrature_points,
        integration_steps=integration_steps,
    )
    stochastic_runs = [
        stochastic_explicit_recent_remote_intervention_summary(
            calibration,
            config=config,
            regimen=regimen,
            population_size=population_size,
            seed=seed_value,
            follow_horizon_years=horizon,
            competing_mortality_hazard=competing_mortality_hazard,
            competing_survival_curve=competing_survival_curve,
        )
        for seed_value in seed_values
    ]
    horizon_key = f"{horizon:g}"
    det_horizon = deterministic["byHorizon"][horizon_key]
    n = float(population_size)
    deterministic_counts = {
        "baselineActiveTB": n
        * float(deterministic["stateProportions"][STATE_BASELINE_ACTIVE_TB]),
        "recentTBI": n * float(deterministic["stateProportions"][STATE_RECENT]),
        "remoteOnlyTBI": n
        * float(deterministic["stateProportions"][STATE_REMOTE_ONLY]),
        "uninfected": n * float(deterministic["stateProportions"][STATE_UNINFECTED]),
        "priorRemotePlusRecent": n
        * float(deterministic["stateProportions"]["priorRemotePlusRecent"]),
        "eligiblePopulation": n * (
            1.0 - float(deterministic["stateProportions"][STATE_BASELINE_ACTIVE_TB])
        ),
        "screened": n * det_horizon["screenedProportion"],
        "notScreened": n * det_horizon["notScreenedProportion"],
        "infectedScreened": n * det_horizon["infectedScreenedProportion"],
        "latentInfectedScreened": n
        * det_horizon["latentInfectedScreenedProportion"],
        "truePositiveLatent": n * det_horizon["truePositiveLatentProportion"],
        "falseNegativeLatent": n * det_horizon["falseNegativeLatentProportion"],
        "uninfectedScreened": n * det_horizon["uninfectedScreenedProportion"],
        "falsePositive": n * det_horizon["falsePositiveProportion"],
        "trueNegative": n * det_horizon["trueNegativeProportion"],
        "activeTBDetectedDuringWorkup": n
        * det_horizon["activeTBDetectedDuringWorkupProportion"],
        "tptEligibleAfterTestingAndAssessment": n
        * det_horizon["tptEligibleAfterTestingAndAssessmentProportion"],
        "tptStarts": n * det_horizon["tptStartsProportion"],
        "truePositiveTPTStarts": n
        * det_horizon["truePositiveTPTStartsProportion"],
        "falsePositiveTPTStarts": n
        * det_horizon["falsePositiveTPTStartsProportion"],
        "tptCompletions": n * det_horizon["tptCompletionsProportion"],
        "effectiveTPT": n * det_horizon["effectiveTPTProportion"],
        "comparatorActiveTB": n * det_horizon["comparatorActiveTBProportion"],
        "interventionActiveTB": n * det_horizon["interventionActiveTBProportion"],
        "directActiveTBCasesAverted": n
        * det_horizon["directActiveTBCasesAvertedProportion"],
        "incidentActiveTBBeforeScreening": n
        * det_horizon["incidentActiveTBBeforeScreeningProportion"],
        "incidentActiveTBAfterScreenBeforeEffectiveTPT": n
        * det_horizon["incidentActiveTBAfterScreenBeforeEffectiveTPTProportion"],
        "subsequentIncidentActiveTB": n
        * det_horizon["subsequentIncidentActiveTBProportion"],
        "recentDirectActiveTBCasesAverted": n
        * det_horizon["recentDirectActiveTBCasesAvertedProportion"],
        "remoteOnlyDirectActiveTBCasesAverted": n
        * det_horizon["remoteOnlyDirectActiveTBCasesAvertedProportion"],
    }
    rows = []
    for metric, expected in deterministic_counts.items():
        values = np.asarray(
            [_stochastic_intervention_metric(run, metric) for run in stochastic_runs],
            dtype=float,
        )
        observed = float(np.mean(values))
        stochastic_sd = (
            float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
        )
        mcse = _monte_carlo_standard_error_counts(values, expected, n)
        tolerance = max(float(absolute_floor_count), float(z_score) * mcse)
        standardized = 0.0 if mcse <= 0.0 else (observed - expected) / mcse
        rows.append(
            {
                "metric": metric,
                "deterministicExpectedCount": float(expected),
                "stochasticMeanCount": observed,
                "stochasticStandardDeviationCount": stochastic_sd,
                "monteCarloStandardErrorCount": mcse,
                "standardizedDifference": standardized,
                "absoluteDifference": abs(observed - expected),
                "minimumStochasticCount": float(np.min(values)),
                "medianStochasticCount": float(np.median(values)),
                "maximumStochasticCount": float(np.max(values)),
                "zeroEventSimulationProportion": float(np.mean(values == 0.0)),
                "tolerance": tolerance,
                "acceptanceZScore": float(z_score),
                "withinTolerance": abs(observed - expected) <= tolerance,
            }
        )
    return {
        "method": INTERVENTION_RECONCILIATION_METHOD,
        "populationSize": int(population_size),
        "horizonYears": horizon,
        "seeds": list(seed_values),
        "stochasticRunCount": len(stochastic_runs),
        "acceptanceRule": (
            "pass when absolute count difference is within "
            "max(absolute_floor_count, z_score * Monte Carlo standard error)"
        ),
        "competingMortality": (
            COMPETING_MORTALITY_CONSTANT_HAZARD
            if competing_mortality_hazard is not None
            else (
                "external_survival_curve"
                if competing_survival_curve is not None
                else COMPETING_MORTALITY_NOT_MODELLED
            )
        ),
        "deterministic": deterministic,
        "stochasticRuns": stochastic_runs,
        "rows": rows,
        "allWithinTolerance": all(row["withinTolerance"] for row in rows),
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
    risk_adjustment = adjustment_for_runner(runner_config)
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
    age_incident_rows = []
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
                    progression_hazard_multiplier=risk_adjustment["multiplier"],
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
                    progression_hazard_multiplier=risk_adjustment["multiplier"],
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
            age_incident_rows.append(
                {
                    "ageYears": float(age),
                    "ageWeight": float(weights[idx]),
                    "horizonYears": horizon,
                    "incidentActiveTBProbabilityTotalPopulation": (
                        at_risk_fraction
                        * float(weights[idx])
                        * (recent_component + remote_component)
                    ),
                    "incidentActiveTBProbabilityWithinAge": (
                        at_risk_fraction * (recent_component + remote_component)
                    ),
                    "recentContributionWithinAge": (
                        at_risk_fraction * recent_component
                    ),
                    "remoteOnlyContributionWithinAge": (
                        at_risk_fraction * remote_component
                    ),
                }
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
        ADJUSTMENT_FIELD: risk_adjustment,
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
        "ageSpecificIncidentActiveTBByHorizon": age_incident_rows,
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
    age_specific_states = []
    age_specific_incident = []
    assigned_state = np.asarray(assigned["explicitTBIState"], dtype=object)
    for age in sorted({float(value) for value in ages}):
        mask = ages == age
        denom = int(mask.sum())
        if denom <= 0:
            continue
        age_specific_states.append(
            {
                "ageYears": age,
                "population": denom,
                "baselineActiveTBProportionWithinAge": float(
                    baseline_active[mask].sum() / denom
                ),
                "recentProportionWithinAge": float(recent[mask].sum() / denom),
                "remoteOnlyProportionWithinAge": float(remote[mask].sum() / denom),
                "uninfectedProportionWithinAge": float(
                    (assigned_state[mask] == STATE_UNINFECTED).sum() / denom
                ),
                "priorRemotePlusRecentProportionWithinAge": float(
                    np.asarray(assigned["priorRemotePlusRecent"], dtype=bool)[
                        mask
                    ].sum()
                    / denom
                ),
            }
        )
        for horizon in horizons:
            incident = (
                at_risk_for_incident
                & (active_time <= horizon)
                & (active_time < death_time)
                & mask
            )
            age_specific_incident.append(
                {
                    "ageYears": age,
                    "population": denom,
                    "horizonYears": horizon,
                    "incidentActiveTBProportionWithinAge": float(
                        incident.sum() / denom
                    ),
                    "incidentActiveTBProportionTotalPopulation": float(
                        incident.sum() / n
                    ),
                }
            )
    return {
        "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
        "method": STOCHASTIC_NATURAL_HISTORY_METHOD,
        ADJUSTMENT_FIELD: adjustment_for_runner(calibration["explicitRecentRemoteRunnerConfig"]),
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
        "ageSpecificStateDistributions": age_specific_states,
        "ageSpecificIncidentActiveTBByHorizon": age_specific_incident,
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
    seeds: tuple[int, ...] | None = None,
    horizons_years: tuple[float, ...] = (1.0, 5.0, 10.0, 20.0),
    z_score: float = 4.0,
    absolute_floor: float = 0.003,
    competing_mortality_hazard: Any = None,
    competing_survival_curve: Any = None,
) -> dict[str, Any]:
    seed_values = tuple(int(value) for value in (seeds or (seed,)))
    if not seed_values:
        raise ValueError("At least one stochastic seed is required.")
    deterministic = deterministic_explicit_recent_remote_expectations(
        calibration,
        horizons_years=horizons_years,
        population_size=population_size,
        competing_mortality_hazard=competing_mortality_hazard,
        competing_survival_curve=competing_survival_curve,
    )
    stochastic_runs = [
        stochastic_explicit_recent_remote_natural_history(
            calibration,
            population_size=population_size,
            seed=seed_value,
            horizons_years=horizons_years,
            competing_mortality_hazard=competing_mortality_hazard,
            competing_survival_curve=competing_survival_curve,
        )
        for seed_value in seed_values
    ]
    rows = []
    n = float(population_size)
    for key, expected in deterministic["stateProportions"].items():
        observed_values = np.asarray(
            [run["stateProportions"][key] for run in stochastic_runs],
            dtype=float,
        )
        observed = float(np.mean(observed_values))
        mcse = _monte_carlo_standard_error(observed_values, expected, n)
        tolerance = max(absolute_floor, z_score * mcse)
        standardized = 0.0 if mcse <= 0.0 else (observed - float(expected)) / mcse
        rows.append(
            {
                "quantity": f"state:{key}",
                "deterministic": expected,
                "stochastic": observed,
                "stochasticMean": observed,
                "monteCarloStandardError": mcse,
                "standardizedDifference": standardized,
                "absoluteDifference": abs(observed - expected),
                "tolerance": tolerance,
                "acceptanceZScore": z_score,
                "withinTolerance": abs(observed - expected) <= tolerance,
            }
        )
    for key, expected_row in deterministic["incidentActiveTBByHorizon"].items():
        expected = float(expected_row["incidentActiveTBProportion"])
        observed_values = np.asarray(
            [
                run["incidentActiveTBByHorizon"][key][
                    "incidentActiveTBProportion"
                ]
                for run in stochastic_runs
            ],
            dtype=float,
        )
        observed = float(np.mean(observed_values))
        mcse = _monte_carlo_standard_error(observed_values, expected, n)
        tolerance = max(absolute_floor, z_score * mcse)
        standardized = 0.0 if mcse <= 0.0 else (observed - expected) / mcse
        rows.append(
            {
                "quantity": f"incidentActiveTB:{key}y",
                "deterministic": expected,
                "stochastic": observed,
                "stochasticMean": observed,
                "monteCarloStandardError": mcse,
                "standardizedDifference": standardized,
                "absoluteDifference": abs(observed - expected),
                "tolerance": tolerance,
                "acceptanceZScore": z_score,
                "withinTolerance": abs(observed - expected) <= tolerance,
            }
        )
    return {
        "method": "deterministic_stochastic_reconciliation_with_mcse_v2",
        "populationSize": int(population_size),
        "seed": int(seed_values[0]),
        "seeds": list(seed_values),
        "stochasticRunCount": len(stochastic_runs),
        "acceptanceRule": (
            "pass when absolute difference is within max(absolute_floor, "
            "z_score * Monte Carlo standard error)"
        ),
        "rows": rows,
        "allWithinTolerance": all(row["withinTolerance"] for row in rows),
        "deterministic": deterministic,
        "stochastic": stochastic_runs[0],
        "stochasticRuns": stochastic_runs,
    }


def _annual_counts_from_times(times: np.ndarray, horizon: float) -> list[dict[str, Any]]:
    finite = np.asarray(times, dtype=float)
    finite = finite[np.isfinite(finite) & (finite <= horizon)]
    rows = []
    for year in range(1, int(math.ceil(horizon)) + 1):
        start = float(year - 1)
        end = min(float(year), horizon)
        rows.append(
            {
                "year": year,
                "intervalStart": start,
                "intervalEnd": end,
                "count": int(((finite > start) & (finite <= end)).sum()),
            }
        )
    return rows


def _age_group_label(age: float) -> str:
    edges = (
        (0.0, 5.0, "0-4"),
        (5.0, 15.0, "5-14"),
        (15.0, 25.0, "15-24"),
        (25.0, 35.0, "25-34"),
        (35.0, 45.0, "35-44"),
        (45.0, 55.0, "45-54"),
        (55.0, 65.0, "55-64"),
        (65.0, 75.0, "65-74"),
        (75.0, 85.0, "75-84"),
    )
    for lower, upper, label in edges:
        if lower <= age < upper:
            return label
    return "85+"


def _stochastic_intervention_metric(summary: Mapping[str, Any], metric: str) -> float:
    if metric == "comparatorActiveTB":
        return float(summary["comparator"]["totalActiveTBOverFollowUp"])
    if metric == "interventionActiveTB":
        return float(summary["intervention"]["totalActiveTBOverFollowUp"])
    if metric == "directActiveTBCasesAverted":
        return float(
            summary["intervention"][
                "directActiveTBCasesAvertedAmongModelledPopulation"
            ]
        )
    if metric == "recentDirectActiveTBCasesAverted":
        return _group_metric(summary.get("byEffectiveState") or [], STATE_RECENT)
    if metric == "remoteOnlyDirectActiveTBCasesAverted":
        return _group_metric(summary.get("byEffectiveState") or [], STATE_REMOTE_ONLY)
    if metric in summary.get("intervention", {}):
        return float(summary["intervention"][metric])
    if metric in summary.get("comparator", {}):
        return float(summary["comparator"][metric])
    raise KeyError(f"Unknown stochastic intervention metric: {metric}")


def _group_metric(rows: list[Mapping[str, Any]], group: str) -> float:
    for row in rows:
        if str(row.get("group")) == group:
            return float(row.get("directActiveTBCasesAverted", 0.0))
    return 0.0


def _monte_carlo_standard_error_counts(
    observed_counts: np.ndarray,
    expected_count: float,
    population_size: float,
) -> float:
    if population_size <= 0.0:
        return 0.0
    observed_proportions = np.asarray(observed_counts, dtype=float) / population_size
    expected_proportion = float(expected_count) / population_size
    return population_size * _monte_carlo_standard_error(
        observed_proportions,
        expected_proportion,
        population_size,
    )


def _explicit_counts_by_group(
    group_values: Any,
    comparator_active: np.ndarray,
    intervention_active: np.ndarray,
    prevented_active_tb: np.ndarray,
    recent: np.ndarray,
    remote: np.ndarray,
    baseline_active: np.ndarray,
) -> list[dict[str, Any]]:
    if group_values is None:
        return []
    groups = np.asarray(group_values)
    if len(groups) != len(comparator_active):
        return []
    rows = []
    for group in sorted({str(value) for value in groups}):
        mask = np.asarray([str(value) == group for value in groups], dtype=bool)
        rows.append(
            {
                "group": group,
                "population": int(mask.sum()),
                "baselineActiveTB": int((baseline_active & mask).sum()),
                "recentTBI": int((recent & mask).sum()),
                "remoteOnlyTBI": int((remote & mask).sum()),
                "comparatorActiveTB": int((comparator_active & mask).sum()),
                "interventionActiveTB": int((intervention_active & mask).sum()),
                "directActiveTBCasesAverted": int(
                    (prevented_active_tb & mask).sum()
                ),
            }
        )
    return rows


def _infected_intervention_moments(
    curve: Mapping[str, Any],
    *,
    age: float,
    horizon: float,
    window_start_years_before_baseline: float,
    window_duration_years: float,
    infection_hazard_per_year: float,
    screen_window_years: float,
    treatment_months: float,
    treatment_paths: tuple[Mapping[str, float], ...],
    infection_quadrature_points: int,
    screen_quadrature_points: int,
    competing_survival_curve: Any,
    competing_mortality_hazard: Any,
    integration_steps: int,
    progression_hazard_multiplier: float = 1.0,
) -> dict[str, float]:
    def by_clock(clock: float) -> dict[str, float]:
        active_by_horizon = _progression_probability_for_clock(
            curve,
            clock,
            horizon,
            age,
            competing_survival_curve=competing_survival_curve,
            competing_mortality_hazard=competing_mortality_hazard,
            integration_steps=integration_steps,
            progression_hazard_multiplier=progression_hazard_multiplier,
        )

        def by_screen(screen_time: float) -> dict[str, float]:
            active_before_screen = _progression_probability_for_clock(
                curve,
                clock,
                min(screen_time, horizon),
                age,
                competing_survival_curve=competing_survival_curve,
                competing_mortality_hazard=competing_mortality_hazard,
                integration_steps=integration_steps,
                progression_hazard_multiplier=progression_hazard_multiplier,
            )
            latent_at_screen = max(0.0, 1.0 - active_before_screen)
            prevented = 0.0
            after_screen_before_effective = 0.0
            for path in treatment_paths:
                path_probability = float(path["probability"])
                if path_probability <= 0.0:
                    continue
                effect_time = screen_time + treatment_months * float(
                    path["doseFraction"]
                ) / 12.0
                active_by_effect = _progression_probability_for_clock(
                    curve,
                    clock,
                    min(effect_time, horizon),
                    age,
                    competing_survival_curve=competing_survival_curve,
                    competing_mortality_hazard=competing_mortality_hazard,
                    integration_steps=integration_steps,
                    progression_hazard_multiplier=progression_hazard_multiplier,
                )
                after_screen_before_effective += path_probability * max(
                    0.0,
                    active_by_effect - active_before_screen,
                )
                prevented += path_probability * float(path["efficacy"]) * max(
                    0.0,
                    active_by_horizon - active_by_effect,
                )
            return {
                "activeBeforeScreening": active_before_screen,
                "latentAtScreen": latent_at_screen,
                "preventedByEffectiveTPT": prevented,
                "activeAfterScreenBeforeEffectiveTPT": (
                    after_screen_before_effective
                ),
            }

        screen_moments = _screen_time_expectation(
            screen_window_years,
            screen_quadrature_points,
            by_screen,
        )
        screen_moments["activeByHorizon"] = active_by_horizon
        return screen_moments

    return _conditional_time_moment_dict(
        window_start_years_before_baseline=window_start_years_before_baseline,
        window_duration_years=window_duration_years,
        hazard_per_year=infection_hazard_per_year,
        quadrature_points=infection_quadrature_points,
        fn=by_clock,
    )


def _screen_time_expectation(
    screen_window_years: float,
    quadrature_points: int,
    fn,
) -> dict[str, float]:
    window = _finite_nonnegative_float(screen_window_years, "screen_window_years")
    if window <= 0.0:
        return dict(fn(0.0))
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_points)
    screen_times = 0.5 * window * (nodes + 1.0)
    totals: dict[str, float] = {}
    for weight, screen_time in zip(weights, screen_times):
        values = fn(float(screen_time))
        for key, value in values.items():
            totals[key] = totals.get(key, 0.0) + 0.5 * float(weight) * float(value)
    return totals


def _conditional_time_moment_dict(
    *,
    window_start_years_before_baseline: float,
    window_duration_years: float,
    hazard_per_year: float,
    quadrature_points: int,
    fn,
) -> dict[str, float]:
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
        return {
            "activeByHorizon": 0.0,
            "activeBeforeScreening": 0.0,
            "latentAtScreen": 0.0,
            "preventedByEffectiveTPT": 0.0,
            "activeAfterScreenBeforeEffectiveTPT": 0.0,
        }
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_points)
    years_after_start = 0.5 * duration * (nodes + 1.0)
    clocks = start + years_after_start
    if hazard <= 0.0:
        densities = np.full_like(years_after_start, 1.0 / duration, dtype=float)
    else:
        normaliser = -math.expm1(-hazard * duration)
        densities = hazard * np.exp(-hazard * years_after_start) / normaliser
    totals: dict[str, float] = {}
    for weight, density, clock in zip(weights, densities, clocks):
        values = fn(float(clock))
        scale = 0.5 * duration * float(weight) * float(density)
        for key, value in values.items():
            totals[key] = totals.get(key, 0.0) + scale * float(value)
    return totals


def _deterministic_test_sensitivity(config: Mapping[str, Any]) -> float:
    test_type = str(config.get("testType", "IGRA")).upper()
    if test_type == "TST":
        return _probability(
            _default_value(config.get("tstSensitivity"), 0.70),
            "tstSensitivity",
        )
    return _probability(
        _default_value(config.get("testSensitivity"), 0.84),
        "testSensitivity",
    )


def _deterministic_test_specificity(config: Mapping[str, Any]) -> float:
    test_type = str(config.get("testType", "IGRA")).upper()
    if test_type != "TST":
        return _probability(
            _default_value(config.get("testSpecificity"), 0.99),
            "testSpecificity",
        )
    risk_prev = config.get("riskPrev") if isinstance(config.get("riskPrev"), Mapping) else {}
    bcg_probability = _probability(risk_prev.get("BCG", 0.0), "riskPrev.BCG")
    spec_bcg = _probability(
        _default_value(config.get("tstSpecificityBCG"), 0.59),
        "tstSpecificityBCG",
    )
    spec_no_bcg = _probability(
        _default_value(config.get("tstSpecificityNoBCG"), 0.97),
        "tstSpecificityNoBCG",
    )
    return bcg_probability * spec_bcg + (1.0 - bcg_probability) * spec_no_bcg


def _default_value(value: Any, default: float) -> Any:
    if value is None or value == []:
        return default
    return value


def _monte_carlo_standard_error(
    observed_values: np.ndarray,
    expected_probability: float,
    population_size: float,
) -> float:
    values = np.asarray(observed_values, dtype=float)
    if len(values) > 1 and float(np.std(values, ddof=1)) > 0.0:
        return float(np.std(values, ddof=1) / math.sqrt(len(values)))
    p = min(max(float(expected_probability), 0.0), 1.0)
    n = max(float(population_size), 1.0)
    return math.sqrt(max(p * (1.0 - p), 0.0) / n)


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
            progression_hazard_multiplier=adjustment_for_runner(runner_config)["multiplier"],
        )
    population["ltbiRiskScore"] = p_tbi
    population["cureTargetScore"] = p_tbi
    population["preventTargetScore"] = p_tbi * future_risk


def _draw_active_times_from_progression_curve(
    states: np.ndarray,
    progression_clock: np.ndarray,
    curve: Mapping[str, Any],
    rng: Any,
    *,
    progression_hazard_multiplier: float = 1.0,
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
            progression_hazard_multiplier=progression_hazard_multiplier,
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
    progression_hazard_multiplier: float = 1.0,
) -> float:
    if competing_survival_curve is None and competing_mortality_hazard is None:
        return time_since_curve_conditional_future_progression_probability(
            curve,
            time_since_infection_at_baseline=clock,
            horizon_years=horizon,
            progression_hazard_multiplier=progression_hazard_multiplier,
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
            progression_hazard_multiplier=progression_hazard_multiplier,
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
            progression_hazard_multiplier=progression_hazard_multiplier,
        )
        return float(result["cumulativeIncidence"])
    result = competing_risk_progression_probability_for_curve(
        curve,
        time_since_infection_at_baseline=clock,
        horizon_years=horizon,
        competing_mortality_hazard=competing_mortality_hazard,
        integration_steps=integration_steps,
        progression_hazard_multiplier=progression_hazard_multiplier,
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


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


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
