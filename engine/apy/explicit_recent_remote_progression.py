from __future__ import annotations

import hashlib
import json
import math
from typing import Any, Iterable, Mapping

from engine.apy.explicit_recent_remote_tbi import (
    STATE_RECENT,
    STATE_REMOTE_ONLY,
    STATE_UNINFECTED,
    validate_active_tb_observation,
)


PROGRESSION_CONTRACT_VERSION = "explicit_recent_remote_tbi_progression_v1"
PROGRESSION_CALIBRATION_CONTRACT_VERSION = (
    "explicit_recent_remote_tbi_progression_calibration_policy_v1"
)

TARGET_BASELINE_PREVALENCE = "baseline_prevalence_target"
TARGET_SCREEN_DETECTED = "screen_detected_disease_target"
TARGET_PROSPECTIVE_INCIDENT = "prospective_incident_disease_target"
TARGET_RETROSPECTIVE_NOTIFICATION = "retrospective_notification_incidence_target"
TARGET_MIXED_OR_INSUFFICIENT = "mixed_or_insufficiently_defined"

POLICY_EXTERNAL_HAZARDS = "external_progression_hazards_v1"
POLICY_FIXED_RATIO_FIT_SCALE = "fixed_early_remote_ratio_fit_scale_v1"
POLICY_VALIDATION_ONLY = "validation_only_v1"
POLICY_JOINT_EARLY_REMOTE_HAZARDS = "joint_early_remote_hazards_v1"

POLICY_BOTH_HAZARDS_SUPPLIED = "both_hazards_externally_supplied"
POLICY_ONE_HAZARD_SUPPLIED = "one_hazard_externally_supplied_fit_other"
POLICY_FIXED_EARLY_LATE_RATIO = "fixed_early_late_ratio_fit_common_scale"
POLICY_JOINT_LIKELIHOOD = "joint_likelihood_multiple_targets"
POLICY_EXTERNAL_VALIDATION_ONLY = "external_validation_only"

OBSERVATION_MODEL_BINOMIAL = "binomial_person_event_v1"
OBSERVATION_MODEL_POISSON = "poisson_count_v1"

RISK_POLICY_NONE = "none"
RISK_POLICY_REVIEWED_HAZARD_MULTIPLIERS = "reviewed_hazard_multipliers"
RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC = "legacy_or_as_hazard_diagnostic_only"

STATE_BASELINE_ACTIVE_TB = "baseline_prevalent_active_tb"

ROOT_TOLERANCE = 1e-12
MAX_ROOT_ITERATIONS = 200
MAX_HAZARD_PER_YEAR = 1e6

DEFAULT_REVIEW_WARNING_MULTIPLIER = 10.0
DEFAULT_REVIEW_BLOCK_MULTIPLIER = 100.0

LEGACY_PROGRESSION_RISK_FACTORS = (
    {
        "factorName": "MJ",
        "displayLabel": "Other current risk-factor prevalence",
        "effect": 3.0,
        "effectMeasure": "OR",
        "source": "inherited APY default disOR",
        "reviewedAsHazardMultiplier": False,
    },
    {
        "factorName": "contact",
        "displayLabel": "Contact history",
        "effect": 5.0,
        "effectMeasure": "OR",
        "source": "inherited APY default disOR",
        "reviewedAsHazardMultiplier": False,
    },
    {
        "factorName": "renal",
        "displayLabel": "Renal impairment",
        "effect": 3.6,
        "effectMeasure": "OR",
        "source": "inherited APY default disOR",
        "reviewedAsHazardMultiplier": False,
    },
    {
        "factorName": "diabetes",
        "displayLabel": "Diabetes",
        "effect": 3.0,
        "effectMeasure": "OR",
        "source": "inherited APY default disOR",
        "reviewedAsHazardMultiplier": False,
    },
    {
        "factorName": "smoking",
        "displayLabel": "Smoking",
        "effect": 2.0,
        "effectMeasure": "OR",
        "source": "inherited APY default disOR",
        "reviewedAsHazardMultiplier": False,
    },
    {
        "factorName": "cld",
        "displayLabel": "Chronic lung disease",
        "effect": 3.0,
        "effectMeasure": "OR",
        "source": "inherited APY default disOR",
        "reviewedAsHazardMultiplier": False,
    },
    {
        "factorName": "alcohol",
        "displayLabel": "Alcohol/drug exposure",
        "effect": 3.0,
        "effectMeasure": "OR",
        "source": "inherited APY default disOR",
        "reviewedAsHazardMultiplier": False,
    },
)


class ProgressionCalibrationError(ValueError):
    def __init__(self, message: str, diagnostics: Mapping[str, Any] | None = None) -> None:
        super().__init__(message)
        self.diagnostics = dict(diagnostics or {})


def build_progression_calibration_policy(
    *,
    policy_id: str,
    hazard_units: str = "per_person_year",
    source: str = "",
    reference_population: str = "",
    review_status: str = "",
    notes: str = "",
    early_hazard: float | None = None,
    remote_hazard: float | None = None,
    early_to_remote_ratio: float | None = None,
    ratio_source: str = "",
    ratio_review_status: str = "",
    ratio_provenance: str = "",
    observation_model_id: str | None = None,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "contractVersion": PROGRESSION_CALIBRATION_CONTRACT_VERSION,
        "policyId": str(policy_id),
        "hazardUnits": str(hazard_units),
        "source": str(source),
        "referencePopulation": str(reference_population),
        "reviewStatus": str(review_status),
        "notes": str(notes),
        "earlyHazard": early_hazard,
        "remoteHazard": remote_hazard,
        "earlyToRemoteRatio": early_to_remote_ratio,
        "ratioSource": str(ratio_source),
        "ratioReviewStatus": str(ratio_review_status),
        "ratioProvenance": str(ratio_provenance),
        "observationModelId": observation_model_id,
    }
    return validate_progression_calibration_policy(payload)


def validate_progression_calibration_policy(payload: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("Progression calibration policy must be a mapping.")
    contract = str(
        payload.get("contractVersion", PROGRESSION_CALIBRATION_CONTRACT_VERSION)
    )
    if contract != PROGRESSION_CALIBRATION_CONTRACT_VERSION:
        raise ValueError(
            "contractVersion must be "
            f"{PROGRESSION_CALIBRATION_CONTRACT_VERSION!r}."
        )
    policy_id = str(payload.get("policyId") or "")
    if policy_id not in {
        POLICY_EXTERNAL_HAZARDS,
        POLICY_FIXED_RATIO_FIT_SCALE,
        POLICY_VALIDATION_ONLY,
        POLICY_JOINT_EARLY_REMOTE_HAZARDS,
    }:
        raise ValueError("Unsupported progression calibration policyId.")
    hazard_units = str(payload.get("hazardUnits") or "").strip()
    if not hazard_units:
        raise ValueError("hazardUnits must be supplied.")
    source = str(payload.get("source") or "").strip()
    reference_population = str(payload.get("referencePopulation") or "").strip()
    review_status = str(payload.get("reviewStatus") or "").strip()
    notes = str(payload.get("notes") or "")
    observation_model_id = payload.get("observationModelId")
    if observation_model_id not in (None, ""):
        observation_model_id = _observation_model_id(observation_model_id)
    else:
        observation_model_id = None

    early_hazard = payload.get("earlyHazard")
    remote_hazard = payload.get("remoteHazard")
    ratio = payload.get("earlyToRemoteRatio")
    ratio_source = str(payload.get("ratioSource") or "").strip()
    ratio_review_status = str(payload.get("ratioReviewStatus") or "").strip()
    ratio_provenance = str(payload.get("ratioProvenance") or "").strip()

    if policy_id == POLICY_EXTERNAL_HAZARDS:
        if early_hazard is None or remote_hazard is None:
            raise ValueError("External-hazard policy requires early and remote hazards.")
        early_hazard = _finite_nonnegative_float(early_hazard, "earlyHazard")
        remote_hazard = _finite_nonnegative_float(remote_hazard, "remoteHazard")
        if not source or not reference_population or not review_status:
            raise ValueError(
                "External-hazard policy requires source, referencePopulation and reviewStatus."
            )
    else:
        early_hazard = None if early_hazard in (None, "") else _finite_nonnegative_float(
            early_hazard,
            "earlyHazard",
        )
        remote_hazard = None if remote_hazard in (None, "") else _finite_nonnegative_float(
            remote_hazard,
            "remoteHazard",
        )

    if policy_id == POLICY_FIXED_RATIO_FIT_SCALE:
        if ratio in (None, ""):
            raise ValueError("Fixed-ratio policy requires earlyToRemoteRatio.")
        ratio = _positive_float(ratio, "earlyToRemoteRatio")
        if not ratio_source or not ratio_review_status or not ratio_provenance:
            raise ValueError(
                "Fixed-ratio policy requires ratio source, review status and provenance."
            )
    else:
        ratio = None if ratio in (None, "") else _positive_float(
            ratio,
            "earlyToRemoteRatio",
        )

    if policy_id == POLICY_JOINT_EARLY_REMOTE_HAZARDS:
        warnings = [
            "Joint early/remote hazard fitting is unavailable unless an explicit "
            "identifiability rank check passes."
        ]
    elif policy_id == POLICY_VALIDATION_ONLY:
        warnings = ["Validation-only policy performs no progression-hazard fitting."]
    else:
        warnings = []

    return {
        "contractVersion": contract,
        "policyId": policy_id,
        "hazardUnits": hazard_units,
        "source": source,
        "referencePopulation": reference_population,
        "reviewStatus": review_status,
        "notes": notes,
        "earlyHazard": early_hazard,
        "remoteHazard": remote_hazard,
        "earlyToRemoteRatio": ratio,
        "ratioSource": ratio_source,
        "ratioReviewStatus": ratio_review_status,
        "ratioProvenance": ratio_provenance,
        "observationModelId": observation_model_id,
        "warnings": warnings,
    }


def progression_calibration_policy_json(policy: Mapping[str, Any]) -> str:
    return _canonical_json(validate_progression_calibration_policy(policy))


def progression_calibration_policy_hash(policy: Mapping[str, Any]) -> str:
    return hashlib.sha256(
        progression_calibration_policy_json(policy).encode("utf-8")
    ).hexdigest()


def remaining_early_risk_years(
    time_since_most_recent_infection: float,
    *,
    recent_window_years: float,
) -> float:
    time_since = _finite_nonnegative_float(
        time_since_most_recent_infection,
        "time_since_most_recent_infection",
    )
    window = _positive_float(recent_window_years, "recent_window_years")
    return max(0.0, window - time_since)


def progression_cumulative_hazard(
    *,
    state: str,
    horizon_years: float,
    early_hazard: float,
    remote_hazard: float,
    multiplier: float = 1.0,
    remaining_early_risk_years: float = 0.0,
) -> float:
    t = _finite_nonnegative_float(horizon_years, "horizon_years")
    early = _finite_nonnegative_float(early_hazard, "early_hazard")
    late = _finite_nonnegative_float(remote_hazard, "remote_hazard")
    mult = _finite_nonnegative_float(multiplier, "multiplier")
    remaining = _finite_nonnegative_float(
        remaining_early_risk_years,
        "remaining_early_risk_years",
    )
    state_key = _state_key(state)
    if state_key == STATE_UNINFECTED:
        return 0.0
    if state_key == STATE_REMOTE_ONLY:
        return mult * late * t
    early_time = min(t, remaining)
    late_time = max(0.0, t - remaining)
    return mult * (early * early_time + late * late_time)


def progression_piecewise_hazard(
    *,
    state: str,
    time_years: float,
    early_hazard: float,
    remote_hazard: float,
    multiplier: float = 1.0,
    remaining_early_risk_years: float = 0.0,
) -> float:
    t = _finite_nonnegative_float(time_years, "time_years")
    early = _finite_nonnegative_float(early_hazard, "early_hazard")
    late = _finite_nonnegative_float(remote_hazard, "remote_hazard")
    mult = _finite_nonnegative_float(multiplier, "multiplier")
    remaining = _finite_nonnegative_float(
        remaining_early_risk_years,
        "remaining_early_risk_years",
    )
    state_key = _state_key(state)
    if state_key == STATE_UNINFECTED:
        return 0.0
    if state_key == STATE_REMOTE_ONLY:
        return mult * late
    return mult * (early if t < remaining else late)


def progression_survival_probability(
    *,
    state: str,
    horizon_years: float,
    early_hazard: float,
    remote_hazard: float,
    multiplier: float = 1.0,
    remaining_early_risk_years: float = 0.0,
) -> float:
    cumulative = progression_cumulative_hazard(
        state=state,
        horizon_years=horizon_years,
        early_hazard=early_hazard,
        remote_hazard=remote_hazard,
        multiplier=multiplier,
        remaining_early_risk_years=remaining_early_risk_years,
    )
    return math.exp(-cumulative)


def progression_probability(
    *,
    state: str,
    horizon_years: float,
    early_hazard: float,
    remote_hazard: float,
    multiplier: float = 1.0,
    remaining_early_risk_years: float = 0.0,
) -> float:
    return 1.0 - progression_survival_probability(
        state=state,
        horizon_years=horizon_years,
        early_hazard=early_hazard,
        remote_hazard=remote_hazard,
        multiplier=multiplier,
        remaining_early_risk_years=remaining_early_risk_years,
    )


def progression_time_quantile(
    *,
    state: str,
    quantile: float,
    early_hazard: float,
    remote_hazard: float,
    multiplier: float = 1.0,
    remaining_early_risk_years: float = 0.0,
) -> float:
    q = _probability(quantile, "quantile")
    state_key = _state_key(state)
    if state_key == STATE_UNINFECTED:
        return math.inf
    if q >= 1.0:
        return math.inf
    if q == 0.0:
        return 0.0
    target = -math.log1p(-q)
    early = _finite_nonnegative_float(early_hazard, "early_hazard")
    late = _finite_nonnegative_float(remote_hazard, "remote_hazard")
    mult = _finite_nonnegative_float(multiplier, "multiplier")
    remaining = _finite_nonnegative_float(
        remaining_early_risk_years,
        "remaining_early_risk_years",
    )
    if mult == 0.0:
        return math.inf
    if state_key == STATE_REMOTE_ONLY:
        rate = mult * late
        return math.inf if rate == 0.0 else target / rate
    early_rate = mult * early
    late_rate = mult * late
    early_cumulative = early_rate * remaining
    if target <= early_cumulative and early_rate > 0.0:
        return target / early_rate
    if late_rate == 0.0:
        return math.inf
    return remaining + (target - early_cumulative) / late_rate


def expected_progression_events(
    strata: Iterable[Mapping[str, Any]],
    *,
    horizon_years: float,
    early_hazard: float,
    remote_hazard: float,
    default_ascertainment_probability: float = 1.0,
    competing_survival_probability: Any = None,
    exclude_baseline_active_tb: bool = True,
) -> dict[str, Any]:
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    ascertainment_default = _probability(
        default_ascertainment_probability,
        "default_ascertainment_probability",
    )
    rows = []
    total = 0.0
    total_weight = 0.0
    included_weight = 0.0
    excluded_baseline_weight = 0.0
    survival_mode = (
        "not_modelled"
        if competing_survival_probability is None
        else "external_survival_probability"
    )
    diagnostic_messages: list[str] = []
    if competing_survival_probability is None:
        diagnostic_messages.append(
            "Competing mortality is not modelled; long-horizon expected cases may be overstated."
        )
    for idx, row in enumerate(strata):
        weight = _row_weight(row, idx)
        raw_state = str(row.get("state", STATE_UNINFECTED))
        if exclude_baseline_active_tb and (
            raw_state == STATE_BASELINE_ACTIVE_TB
            or bool(row.get("baselineActiveTB", False))
        ):
            excluded_baseline_weight += weight
            total_weight += weight
            continue
        state = _state_key(raw_state)
        multiplier = _finite_nonnegative_float(
            row.get("multiplier", 1.0),
            f"strata[{idx}].multiplier",
        )
        remaining = _finite_nonnegative_float(
            row.get("remainingEarlyRiskYears", row.get("remaining_early_risk_years", 0.0)),
            f"strata[{idx}].remainingEarlyRiskYears",
        )
        ascertainment = _probability(
            row.get("ascertainmentProbability", ascertainment_default),
            f"strata[{idx}].ascertainmentProbability",
        )
        competing_survival = _survival_probability_at(
            competing_survival_probability,
            horizon,
            row,
            f"strata[{idx}].competingSurvivalProbability",
        )
        probability = progression_probability(
            state=state,
            horizon_years=horizon,
            early_hazard=early_hazard,
            remote_hazard=remote_hazard,
            multiplier=multiplier,
            remaining_early_risk_years=remaining,
        )
        adjusted_probability = probability * competing_survival
        expected = weight * adjusted_probability * ascertainment
        total += expected
        total_weight += weight
        included_weight += weight
        rows.append(
            {
                "state": state,
                "weight": weight,
                "multiplier": multiplier,
                "remainingEarlyRiskYears": remaining,
                "ascertainmentProbability": ascertainment,
                "competingSurvivalProbability": competing_survival,
                "progressionProbability": probability,
                "survivalAdjustedProgressionProbability": adjusted_probability,
                "expectedCases": expected,
            }
        )
    return {
        "contractVersion": PROGRESSION_CONTRACT_VERSION,
        "horizonYears": horizon,
        "earlyHazard": _finite_nonnegative_float(early_hazard, "early_hazard"),
        "remoteHazard": _finite_nonnegative_float(remote_hazard, "remote_hazard"),
        "expectedCases": total,
        "totalWeight": total_weight,
        "includedWeight": included_weight,
        "excludedBaselineActiveTBWeight": excluded_baseline_weight,
        "competingMortality": survival_mode,
        "diagnosticMessages": diagnostic_messages,
        "rows": rows,
    }


def classify_active_tb_observation_target(
    row: Mapping[str, Any],
    *,
    model_baseline_year: int | None = None,
) -> dict[str, Any]:
    observed = validate_active_tb_observation(row)
    meaning = observed["observationWindowMeaning"]
    classification = observed["caseClassification"]
    ascertainment = observed["ascertainmentMethod"]
    target_type = TARGET_MIXED_OR_INSUFFICIENT
    reason = "Observation timing is insufficient for prospective progression calibration."
    if meaning == "baseline_prevalent" or classification == "prevalent_baseline":
        target_type = TARGET_BASELINE_PREVALENCE
        reason = "Prevalent active TB at or near baseline is not prospective incident progression."
    elif meaning == "screen_detected_prevalent" or classification == "screen_detected":
        target_type = TARGET_SCREEN_DETECTED
        reason = "Screen-detected active TB is a screening outcome, not a future incident count."
    elif meaning == "follow_up_incident" or classification == "incident_follow_up":
        if model_baseline_year is None:
            reason = "Model baseline year is required to distinguish retrospective from prospective incidence."
        elif observed["endYear"] < int(model_baseline_year):
            target_type = TARGET_RETROSPECTIVE_NOTIFICATION
            reason = "Incident notification window ends before model baseline."
        elif observed["startYear"] == int(model_baseline_year):
            target_type = TARGET_PROSPECTIVE_INCIDENT
            reason = "Incident observation window begins at model baseline."
        elif observed["startYear"] < int(model_baseline_year):
            target_type = TARGET_MIXED_OR_INSUFFICIENT
            reason = "Observation window spans before and after model baseline."
        else:
            target_type = TARGET_MIXED_OR_INSUFFICIENT
            reason = "Observation window does not begin at model baseline."
    if (
        target_type == TARGET_RETROSPECTIVE_NOTIFICATION
        and ascertainment != "passive_notification"
    ):
        reason = "Retrospective incident window before baseline; ascertainment is not enough to reconstruct a baseline cohort."
    return {
        "schemaVersion": observed["schemaVersion"],
        "targetType": target_type,
        "reason": reason,
        "observation": observed,
    }


def assess_observation_prospective_calibration_eligibility(
    row: Mapping[str, Any],
    *,
    model_baseline_year: int,
    strata: Iterable[Mapping[str, Any]] | None = None,
    ascertainment_probability: float | None = None,
) -> dict[str, Any]:
    classified = classify_active_tb_observation_target(
        row,
        model_baseline_year=model_baseline_year,
    )
    observed = classified["observation"]
    reasons: list[str] = []
    warnings: list[str] = []
    if classified["targetType"] != TARGET_PROSPECTIVE_INCIDENT:
        reasons.append(classified["reason"])
    if observed["caseClassification"] != "incident_follow_up":
        reasons.append("Numerator is not classified as incident follow-up disease.")
    if observed["observationWindowMeaning"] != "follow_up_incident":
        reasons.append("Observation window is not prospective follow-up incidence.")
    if observed["startYear"] != int(model_baseline_year):
        reasons.append("Observation window must begin at model baseline.")
    horizon = _observation_duration_years(observed)
    if horizon <= 0.0:
        reasons.append("Observation duration must be positive and explicit.")
    if ascertainment_probability is None:
        reasons.append("An explicit ascertainment probability is required for fitting.")
    else:
        _probability(ascertainment_probability, "ascertainment_probability")
        if ascertainment_probability < 1.0:
            warnings.append("Incomplete ascertainment is modelled explicitly through q.")
    if strata is None:
        reasons.append("Compatible baseline population composition/strata are required.")
        included_weight = 0.0
        excluded_baseline = 0.0
    else:
        prepared = _prepare_progression_strata(strata)
        included_weight = sum(row_out["weight"] for row_out in prepared["includedRows"])
        excluded_baseline = prepared["excludedBaselineActiveTBWeight"]
        if included_weight <= 0.0:
            reasons.append("Eligible non-baseline-active-TB population weight must be positive.")
        if prepared["invalidStates"]:
            reasons.extend(prepared["invalidStates"])
        if excluded_baseline > 0.0:
            warnings.append(
                "Baseline active-TB strata are excluded from prospective progression calibration."
            )
    if observed["observedActiveTBCaseCount"] > observed["populationDenominator"]:
        reasons.append("Observed cases cannot exceed the population denominator.")
    return {
        "contractVersion": PROGRESSION_CALIBRATION_CONTRACT_VERSION,
        "eligible": len(reasons) == 0,
        "targetType": classified["targetType"],
        "reasons": reasons if reasons else ["Eligible prospective incident calibration target."],
        "warnings": warnings,
        "observation": observed,
        "horizonYears": horizon,
        "includedPopulationWeight": included_weight,
        "excludedBaselineActiveTBWeight": excluded_baseline,
    }


def expected_prospective_incident_cases_for_observation(
    row: Mapping[str, Any],
    strata: Iterable[Mapping[str, Any]],
    *,
    model_baseline_year: int,
    early_hazard: float,
    remote_hazard: float,
    ascertainment_probability: float = 1.0,
    competing_survival_function: Any = None,
    competing_survival_probability: Any = None,
) -> dict[str, Any]:
    if competing_survival_function is not None and competing_survival_probability is not None:
        raise ValueError(
            "Provide either competing_survival_function or competing_survival_probability, not both."
        )
    survival_spec = (
        competing_survival_function
        if competing_survival_function is not None
        else competing_survival_probability
    )
    strata_tuple = tuple(strata)
    eligibility = assess_observation_prospective_calibration_eligibility(
        row,
        model_baseline_year=model_baseline_year,
        strata=strata_tuple,
        ascertainment_probability=ascertainment_probability,
    )
    if not eligibility["eligible"]:
        raise ValueError(
            "Only genuinely prospective incident targets can be used for this expected-count function."
        )
    observed = eligibility["observation"]
    horizon = eligibility["horizonYears"]
    expected = expected_progression_events(
        strata_tuple,
        horizon_years=horizon,
        early_hazard=early_hazard,
        remote_hazard=remote_hazard,
        default_ascertainment_probability=ascertainment_probability,
        competing_survival_probability=survival_spec,
    )
    return {
        "contractVersion": PROGRESSION_CONTRACT_VERSION,
        "calibrationContractVersion": PROGRESSION_CALIBRATION_CONTRACT_VERSION,
        "targetType": eligibility["targetType"],
        "observationId": observed["observationId"],
        "observationStartYear": observed["startYear"],
        "observationEndYear": observed["endYear"],
        "observationHorizonYears": horizon,
        "observedActiveTBCaseCount": observed["observedActiveTBCaseCount"],
        "populationDenominator": observed["populationDenominator"],
        "personYears": observed["personYears"],
        "ascertainmentProbability": _probability(
            ascertainment_probability,
            "ascertainment_probability",
        ),
        "expectedActiveTBCaseCount": expected["expectedCases"],
        "expectedProgression": expected,
        "competingMortality": expected["competingMortality"],
        "competingMortalityApplied": expected["competingMortality"] != "not_modelled",
        "eligibility": eligibility,
    }


def evaluate_external_hazard_policy(
    row: Mapping[str, Any],
    strata: Iterable[Mapping[str, Any]],
    policy: Mapping[str, Any],
    *,
    model_baseline_year: int,
    ascertainment_probability: float = 1.0,
    competing_survival_probability: Any = None,
) -> dict[str, Any]:
    validated = validate_progression_calibration_policy(policy)
    if validated["policyId"] != POLICY_EXTERNAL_HAZARDS:
        raise ValueError("External-hazard evaluation requires external-hazard policy.")
    expected = expected_prospective_incident_cases_for_observation(
        row,
        strata,
        model_baseline_year=model_baseline_year,
        early_hazard=float(validated["earlyHazard"]),
        remote_hazard=float(validated["remoteHazard"]),
        ascertainment_probability=ascertainment_probability,
        competing_survival_probability=competing_survival_probability,
    )
    return {
        **expected,
        "policyId": POLICY_EXTERNAL_HAZARDS,
        "fittingPerformed": False,
        "policy": validated,
        "earlyHazard": validated["earlyHazard"],
        "remoteHazard": validated["remoteHazard"],
        "diagnosticMessages": ["External hazards supplied; no fitting performed."],
    }


def evaluate_validation_only_policy(
    row: Mapping[str, Any],
    strata: Iterable[Mapping[str, Any]],
    policy: Mapping[str, Any],
    *,
    model_baseline_year: int,
    early_hazard: float,
    remote_hazard: float,
    ascertainment_probability: float = 1.0,
    competing_survival_probability: Any = None,
) -> dict[str, Any]:
    validated = validate_progression_calibration_policy(policy)
    if validated["policyId"] != POLICY_VALIDATION_ONLY:
        raise ValueError("Validation-only evaluation requires validation-only policy.")
    classified = classify_active_tb_observation_target(
        row,
        model_baseline_year=model_baseline_year,
    )
    expected = None
    messages = ["Validation-only policy performs no fitting and does not modify hazards."]
    if classified["targetType"] == TARGET_PROSPECTIVE_INCIDENT:
        expected = expected_prospective_incident_cases_for_observation(
            row,
            strata,
            model_baseline_year=model_baseline_year,
            early_hazard=early_hazard,
            remote_hazard=remote_hazard,
            ascertainment_probability=ascertainment_probability,
            competing_survival_probability=competing_survival_probability,
        )
    else:
        messages.append(classified["reason"])
    return {
        "contractVersion": PROGRESSION_CALIBRATION_CONTRACT_VERSION,
        "policyId": POLICY_VALIDATION_ONLY,
        "fittingPerformed": False,
        "policy": validated,
        "targetType": classified["targetType"],
        "classification": classified,
        "earlyHazard": _finite_nonnegative_float(early_hazard, "early_hazard"),
        "remoteHazard": _finite_nonnegative_float(remote_hazard, "remote_hazard"),
        "expectedIncidentCases": expected,
        "diagnosticMessages": messages,
    }


def fit_fixed_ratio_progression_scale(
    row: Mapping[str, Any],
    strata: Iterable[Mapping[str, Any]],
    *,
    model_baseline_year: int,
    early_to_remote_ratio: float,
    ratio_source: str,
    ratio_review_status: str,
    ratio_provenance: str,
    ascertainment_probability: float = 1.0,
    competing_survival_probability: Any = None,
    tolerance: float = ROOT_TOLERANCE,
) -> dict[str, Any]:
    ratio = _positive_float(early_to_remote_ratio, "early_to_remote_ratio")
    if not str(ratio_source).strip():
        raise ValueError("ratio_source must be supplied.")
    if not str(ratio_review_status).strip():
        raise ValueError("ratio_review_status must be supplied.")
    if not str(ratio_provenance).strip():
        raise ValueError("ratio_provenance must be supplied.")
    strata_tuple = tuple(strata)
    eligibility = assess_observation_prospective_calibration_eligibility(
        row,
        model_baseline_year=model_baseline_year,
        strata=strata_tuple,
        ascertainment_probability=ascertainment_probability,
    )
    if not eligibility["eligible"]:
        diagnostics = {
            "contractVersion": PROGRESSION_CALIBRATION_CONTRACT_VERSION,
            "policyId": POLICY_FIXED_RATIO_FIT_SCALE,
            "feasibilityStatus": "ineligible_target",
            "eligibility": eligibility,
        }
        raise ProgressionCalibrationError(
            "Active-TB observation is not eligible for fixed-ratio progression calibration.",
            diagnostics,
        )
    observed = eligibility["observation"]
    requested = float(observed["observedActiveTBCaseCount"])
    denominator = float(observed["populationDenominator"])
    horizon = float(eligibility["horizonYears"])

    def expected_for_scale(scale: float) -> float:
        return expected_progression_events(
            strata_tuple,
            horizon_years=horizon,
            early_hazard=ratio * scale,
            remote_hazard=scale,
            default_ascertainment_probability=ascertainment_probability,
            competing_survival_probability=competing_survival_probability,
        )["expectedCases"]

    maximum = expected_for_scale(MAX_HAZARD_PER_YEAR)
    warnings = list(eligibility["warnings"])
    if requested > maximum + max(float(tolerance), 1e-12):
        diagnostics = _fixed_ratio_result(
            requested=requested,
            denominator=denominator,
            horizon=horizon,
            ascertainment_probability=ascertainment_probability,
            ratio=ratio,
            fitted_remote_hazard=None,
            achieved_expected_cases=maximum,
            convergence_status="not_run",
            feasibility_status="infeasible_above_achievable_range",
            warnings=warnings
            + ["Requested cases exceed the model-achievable range under this policy."],
            ratio_source=ratio_source,
            ratio_review_status=ratio_review_status,
            ratio_provenance=ratio_provenance,
            eligibility=eligibility,
            expected_progression=None,
            concentration=None,
        )
        raise ProgressionCalibrationError(
            "Requested active-TB count is above the achievable range.",
            diagnostics,
        )
    if requested <= float(tolerance):
        scale = 0.0
        convergence_status = "converged_zero_target"
    else:
        lo = 0.0
        hi = 1.0
        while expected_for_scale(hi) < requested and hi < MAX_HAZARD_PER_YEAR:
            hi *= 2.0
        if expected_for_scale(hi) < requested:
            diagnostics = _fixed_ratio_result(
                requested=requested,
                denominator=denominator,
                horizon=horizon,
                ascertainment_probability=ascertainment_probability,
                ratio=ratio,
                fitted_remote_hazard=None,
                achieved_expected_cases=expected_for_scale(hi),
                convergence_status="not_converged",
                feasibility_status="unbracketed",
                warnings=warnings
                + ["Could not bracket fitted scale below maximum hazard bound."],
                ratio_source=ratio_source,
                ratio_review_status=ratio_review_status,
                ratio_provenance=ratio_provenance,
                eligibility=eligibility,
                expected_progression=None,
                concentration=None,
            )
            raise ProgressionCalibrationError(
                "Could not bracket fixed-ratio progression scale.",
                diagnostics,
            )
        scale = hi
        convergence_status = "not_converged"
        for _ in range(MAX_ROOT_ITERATIONS):
            mid = (lo + hi) / 2.0
            value = expected_for_scale(mid)
            if abs(value - requested) <= float(tolerance) or (hi - lo) / 2.0 <= float(
                tolerance
            ):
                scale = mid
                convergence_status = "converged"
                break
            if value < requested:
                lo = mid
            else:
                hi = mid
        else:
            scale = (lo + hi) / 2.0

    expected = expected_progression_events(
        strata_tuple,
        horizon_years=horizon,
        early_hazard=ratio * scale,
        remote_hazard=scale,
        default_ascertainment_probability=ascertainment_probability,
        competing_survival_probability=competing_survival_probability,
    )
    concentration = progression_case_concentration_diagnostics(expected["rows"])
    if concentration["shareExpectedCasesTop1PercentByWeight"] > 0.5:
        warnings.append("Highest 1% of eligible weight contributes more than half of expected cases.")
    return _fixed_ratio_result(
        requested=requested,
        denominator=denominator,
        horizon=horizon,
        ascertainment_probability=ascertainment_probability,
        ratio=ratio,
        fitted_remote_hazard=scale,
        achieved_expected_cases=expected["expectedCases"],
        convergence_status=convergence_status,
        feasibility_status="feasible",
        warnings=warnings,
        ratio_source=ratio_source,
        ratio_review_status=ratio_review_status,
        ratio_provenance=ratio_provenance,
        eligibility=eligibility,
        expected_progression=expected,
        concentration=concentration,
    )


def active_tb_observation_log_likelihood(
    *,
    observed_cases: float,
    expected_cases: float,
    denominator: float,
    observation_model_id: str,
) -> dict[str, Any]:
    observed = _finite_nonnegative_float(observed_cases, "observed_cases")
    expected = _finite_nonnegative_float(expected_cases, "expected_cases")
    denom = _positive_float(denominator, "denominator")
    model_id = _observation_model_id(observation_model_id)
    if abs(observed - round(observed)) > 1e-12:
        raise ValueError("observed_cases must be an integer count for likelihoods.")
    observed_int = int(round(observed))
    if model_id == OBSERVATION_MODEL_BINOMIAL:
        if abs(denom - round(denom)) > 1e-12:
            raise ValueError("denominator must be an integer for binomial likelihood.")
        n = int(round(denom))
        if observed_int > n:
            raise ValueError("observed_cases cannot exceed denominator.")
        if expected > denom + 1e-12:
            raise ValueError("expected_cases cannot exceed denominator for binomial likelihood.")
        p = expected / denom
        if p <= 0.0:
            log_likelihood = 0.0 if observed_int == 0 else -math.inf
        elif p >= 1.0:
            log_likelihood = 0.0 if observed_int == n else -math.inf
        else:
            log_likelihood = (
                math.lgamma(n + 1)
                - math.lgamma(observed_int + 1)
                - math.lgamma(n - observed_int + 1)
                + observed_int * math.log(p)
                + (n - observed_int) * math.log1p(-p)
            )
    else:
        if expected == 0.0:
            log_likelihood = 0.0 if observed_int == 0 else -math.inf
        else:
            log_likelihood = (
                observed_int * math.log(expected)
                - expected
                - math.lgamma(observed_int + 1)
            )
    return {
        "contractVersion": PROGRESSION_CALIBRATION_CONTRACT_VERSION,
        "observationModelId": model_id,
        "observedCases": observed,
        "expectedCases": expected,
        "denominator": denom,
        "logLikelihood": log_likelihood,
        "limitations": (
            "Observation-model choice is explicit; recurrent disease, migration "
            "or changing population composition can invalidate this likelihood."
        ),
    }


def progression_identifiability_assessment(
    targets: Iterable[Mapping[str, Any]],
    *,
    free_parameters: tuple[str, ...] = ("lambda_early", "lambda_late"),
    model_baseline_year: int | None = None,
) -> dict[str, Any]:
    classified = [
        classify_active_tb_observation_target(target, model_baseline_year=model_baseline_year)
        for target in targets
    ]
    informative = [
        item
        for item in classified
        if item["targetType"] == TARGET_PROSPECTIVE_INCIDENT
    ]
    free = tuple(str(parameter) for parameter in free_parameters)
    sufficient = len(informative) >= len(free)
    messages = []
    if len(free) >= 2 and len(informative) < 2:
        messages.append(
            "One aggregate prospective active-TB count generally cannot identify both lambda_early and lambda_late."
        )
    if not informative:
        messages.append("No prospective incident target is available for progression calibration.")
    return {
        "contractVersion": PROGRESSION_CONTRACT_VERSION,
        "freeParameters": list(free),
        "prospectiveIncidentTargetCount": len(informative),
        "isIdentifiedByAggregateTargets": sufficient,
        "identifiabilityStatus": "identified_under_policy" if sufficient else "insufficient",
        "classifiedTargets": classified,
        "diagnosticMessages": messages,
    }


def progression_calibration_policy_options() -> tuple[dict[str, str], ...]:
    return (
        {
            "policyId": POLICY_BOTH_HAZARDS_SUPPLIED,
            "description": "Both progression hazards are externally supplied; active-TB observations are validation data.",
        },
        {
            "policyId": POLICY_ONE_HAZARD_SUPPLIED,
            "description": "One progression hazard is externally supplied and the other is fitted.",
        },
        {
            "policyId": POLICY_FIXED_EARLY_LATE_RATIO,
            "description": "An externally reviewed early-to-late hazard ratio is fixed and a common scale is fitted.",
        },
        {
            "policyId": POLICY_JOINT_LIKELIHOOD,
            "description": "Multiple sufficiently informative targets are used in a joint likelihood.",
        },
        {
            "policyId": POLICY_EXTERNAL_VALIDATION_ONLY,
            "description": "Active-TB data are used only for external validation.",
        },
    )


def progression_calibration_policy_contract_options() -> tuple[dict[str, str], ...]:
    return (
        {
            "policyId": POLICY_EXTERNAL_HAZARDS,
            "description": "Policy A: early and remote progression hazards are supplied externally; no fitting occurs.",
        },
        {
            "policyId": POLICY_FIXED_RATIO_FIT_SCALE,
            "description": "Policy B: lambda_L=k and lambda_E=Rk with externally supplied R and fitted k.",
        },
        {
            "policyId": POLICY_VALIDATION_ONLY,
            "description": "Policy C: observations are compared with predictions but no hazard fitting occurs.",
        },
        {
            "policyId": POLICY_JOINT_EARLY_REMOTE_HAZARDS,
            "description": "Policy D: joint early/remote fitting is unavailable unless identifiability rank criteria pass.",
        },
    )


def assess_progression_policy_identifiability(
    policy_id: str,
    targets: Iterable[Mapping[str, Any]],
    *,
    model_baseline_year: int,
    strata_by_target: Iterable[Iterable[Mapping[str, Any]]] | None = None,
    ascertainment_probability: float = 1.0,
    early_to_remote_ratio: float | None = None,
) -> dict[str, Any]:
    policy = str(policy_id)
    target_tuple = tuple(targets)
    strata_groups = None if strata_by_target is None else tuple(tuple(group) for group in strata_by_target)
    classified = [
        classify_active_tb_observation_target(target, model_baseline_year=model_baseline_year)
        for target in target_tuple
    ]
    eligible = [
        item for item in classified if item["targetType"] == TARGET_PROSPECTIVE_INCIDENT
    ]
    messages: list[str] = []
    sensitivity_vectors: list[tuple[float, float]] = []
    rank = 0
    if strata_groups is not None:
        for target, strata in zip(target_tuple, strata_groups):
            assessment = assess_observation_prospective_calibration_eligibility(
                target,
                model_baseline_year=model_baseline_year,
                strata=strata,
                ascertainment_probability=ascertainment_probability,
            )
            if assessment["eligible"]:
                sensitivity_vectors.append(
                    target_progression_sensitivity_vector(
                        strata,
                        horizon_years=assessment["horizonYears"],
                        ascertainment_probability=ascertainment_probability,
                    )
                )
        rank = _rank_2d(sensitivity_vectors)
    if policy == POLICY_EXTERNAL_HAZARDS:
        status = "no_fitting_required"
        identified = True
        messages.append("Externally supplied hazards require no statistical identification.")
    elif policy == POLICY_FIXED_RATIO_FIT_SCALE:
        _positive_float(early_to_remote_ratio, "early_to_remote_ratio")
        identified = len(eligible) >= 1
        status = "identified_one_dimensional" if identified else "insufficient_targets"
        messages.append(
            "Fixed early-to-remote ratio leaves one fitted scale parameter."
        )
        if not identified:
            messages.append("At least one eligible prospective incident target is required.")
    elif policy == POLICY_VALIDATION_ONLY:
        status = "no_fitting_requested"
        identified = True
        messages.append("Validation-only policy does not estimate progression hazards.")
    elif policy == POLICY_JOINT_EARLY_REMOTE_HAZARDS:
        identified = len(eligible) >= 2 and rank >= 2
        status = "available" if identified else "unavailable_nonidentifiable"
        if len(eligible) < 2:
            messages.append(
                "One aggregate target cannot identify both lambda_E and lambda_L."
            )
        if len(eligible) >= 2 and rank < 2:
            messages.append(
                "Targets with identical or near-identical recent/remote sensitivity do not add independent information."
            )
    else:
        raise ValueError("Unsupported policy_id for identifiability assessment.")
    return {
        "contractVersion": PROGRESSION_CALIBRATION_CONTRACT_VERSION,
        "policyId": policy,
        "identified": identified,
        "identifiabilityStatus": status,
        "prospectiveIncidentTargetCount": len(eligible),
        "sensitivityVectors": [
            {"early": vector[0], "remote": vector[1]} for vector in sensitivity_vectors
        ],
        "sensitivityRank": rank,
        "classifiedTargets": classified,
        "diagnosticMessages": messages,
    }


def target_progression_sensitivity_vector(
    strata: Iterable[Mapping[str, Any]],
    *,
    horizon_years: float,
    ascertainment_probability: float = 1.0,
    competing_survival_probability: Any = None,
) -> tuple[float, float]:
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    q_default = _probability(ascertainment_probability, "ascertainment_probability")
    early = 0.0
    remote = 0.0
    for idx, row in enumerate(strata):
        raw_state = str(row.get("state", STATE_UNINFECTED))
        if raw_state == STATE_BASELINE_ACTIVE_TB or bool(row.get("baselineActiveTB", False)):
            continue
        state = _state_key(raw_state)
        weight = _row_weight(row, idx)
        multiplier = _finite_nonnegative_float(
            row.get("multiplier", 1.0),
            f"strata[{idx}].multiplier",
        )
        remaining = _finite_nonnegative_float(
            row.get("remainingEarlyRiskYears", row.get("remaining_early_risk_years", 0.0)),
            f"strata[{idx}].remainingEarlyRiskYears",
        )
        ascertainment = _probability(
            row.get("ascertainmentProbability", q_default),
            f"strata[{idx}].ascertainmentProbability",
        )
        survival = _survival_probability_at(
            competing_survival_probability,
            horizon,
            row,
            f"strata[{idx}].competingSurvivalProbability",
        )
        if state == STATE_RECENT:
            early += weight * multiplier * min(horizon, remaining) * ascertainment * survival
            remote += weight * multiplier * max(0.0, horizon - remaining) * ascertainment * survival
        elif state == STATE_REMOTE_ONLY:
            remote += weight * multiplier * horizon * ascertainment * survival
    return (early, remote)


def build_risk_factor_application_policy(
    *,
    policy_id: str,
    reviewed_factor_names: Iterable[str] = (),
    review_warning_threshold: float = DEFAULT_REVIEW_WARNING_MULTIPLIER,
    review_block_threshold: float = DEFAULT_REVIEW_BLOCK_MULTIPLIER,
) -> dict[str, Any]:
    policy = str(policy_id)
    if policy not in {
        RISK_POLICY_NONE,
        RISK_POLICY_REVIEWED_HAZARD_MULTIPLIERS,
        RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC,
    }:
        raise ValueError("Unsupported risk-factor application policy.")
    warning = _positive_float(review_warning_threshold, "review_warning_threshold")
    blocking = _positive_float(review_block_threshold, "review_block_threshold")
    if blocking < warning:
        raise ValueError("review_block_threshold must be at least review_warning_threshold.")
    reviewed = tuple(str(name) for name in reviewed_factor_names)
    return {
        "policyId": policy,
        "reviewedFactorNames": list(reviewed),
        "reviewWarningThreshold": warning,
        "reviewBlockThreshold": blocking,
        "scientificStatus": (
            "scientifically_provisional_diagnostic"
            if policy == RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC
            else "requires_explicit_review"
        ),
    }


def risk_factor_multiplier_for_row(
    row: Mapping[str, Any],
    policy: Mapping[str, Any],
    *,
    factor_effects: Iterable[Mapping[str, Any]] = LEGACY_PROGRESSION_RISK_FACTORS,
) -> dict[str, Any]:
    validated_policy = build_risk_factor_application_policy(
        policy_id=str(policy.get("policyId")),
        reviewed_factor_names=policy.get("reviewedFactorNames", ()),
        review_warning_threshold=policy.get(
            "reviewWarningThreshold",
            DEFAULT_REVIEW_WARNING_MULTIPLIER,
        ),
        review_block_threshold=policy.get(
            "reviewBlockThreshold",
            DEFAULT_REVIEW_BLOCK_MULTIPLIER,
        ),
    )
    reviewed_names = set(validated_policy["reviewedFactorNames"])
    multiplier = 1.0
    contributions = []
    for factor in factor_effects:
        name = str(factor["factorName"])
        effect = _positive_float(factor["effect"], f"{name}.effect")
        present = bool(row.get(name, row.get(f"{name}Flag", False)))
        include = False
        reason = "not_present"
        if present and validated_policy["policyId"] == RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC:
            include = True
            reason = "legacy_or_as_hazard_diagnostic_only"
        elif (
            present
            and validated_policy["policyId"] == RISK_POLICY_REVIEWED_HAZARD_MULTIPLIERS
            and (
                name in reviewed_names
                or bool(factor.get("reviewedAsHazardMultiplier", False))
            )
        ):
            include = True
            reason = "reviewed_hazard_multiplier"
        elif present and validated_policy["policyId"] == RISK_POLICY_REVIEWED_HAZARD_MULTIPLIERS:
            reason = "not_reviewed_as_hazard_multiplier"
        if include:
            multiplier *= effect
        contributions.append(
            {
                "factorName": name,
                "displayLabel": str(factor.get("displayLabel", name)),
                "present": present,
                "included": include,
                "effect": effect,
                "effectMeasure": str(factor.get("effectMeasure", "unclear")),
                "source": str(factor.get("source", "")),
                "logContribution": math.log(effect) if include else 0.0,
                "reason": reason,
            }
        )
    status = "ok"
    warnings = []
    if multiplier >= validated_policy["reviewBlockThreshold"]:
        status = "blocking_review_required"
        warnings.append("Combined progression multiplier exceeds blocking review threshold.")
    elif multiplier >= validated_policy["reviewWarningThreshold"]:
        status = "warning_review_required"
        warnings.append("Combined progression multiplier exceeds warning review threshold.")
    if validated_policy["policyId"] == RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC:
        warnings.append(
            "Legacy OR-as-hazard multiplication is diagnostic only and scientifically provisional."
        )
    return {
        "policy": validated_policy,
        "combinedMultiplier": multiplier,
        "contributions": contributions,
        "reviewStatus": status,
        "warnings": warnings,
    }


def risk_factor_multiplier_diagnostics(
    population_rows: Iterable[Mapping[str, Any]],
    policy: Mapping[str, Any],
    *,
    horizon_years: float,
    early_hazard: float,
    remote_hazard: float,
    factor_effects: Iterable[Mapping[str, Any]] = LEGACY_PROGRESSION_RISK_FACTORS,
    review_levels: Iterable[float] = (0.1, 0.2, 0.5),
) -> dict[str, Any]:
    rows = []
    for idx, row in enumerate(population_rows):
        risk = risk_factor_multiplier_for_row(
            row,
            policy,
            factor_effects=factor_effects,
        )
        state = _state_key(str(row.get("state", STATE_UNINFECTED)))
        weight = _row_weight(row, idx)
        remaining = _finite_nonnegative_float(
            row.get("remainingEarlyRiskYears", row.get("remaining_early_risk_years", 0.0)),
            f"population_rows[{idx}].remainingEarlyRiskYears",
        )
        probability = progression_probability(
            state=state,
            horizon_years=horizon_years,
            early_hazard=early_hazard,
            remote_hazard=remote_hazard,
            multiplier=risk["combinedMultiplier"],
            remaining_early_risk_years=remaining,
        )
        rows.append(
            {
                "index": idx,
                "state": state,
                "weight": weight,
                "combinedMultiplier": risk["combinedMultiplier"],
                "progressionProbability": probability,
                "expectedCases": weight * probability,
                "riskDiagnostic": risk,
            }
        )
    multipliers = [row["combinedMultiplier"] for row in rows]
    total_weight = sum(row["weight"] for row in rows)
    total_cases = sum(row["expectedCases"] for row in rows)
    policy_validated = (
        rows[0]["riskDiagnostic"]["policy"]
        if rows
        else build_risk_factor_application_policy(policy_id=str(policy.get("policyId")))
    )
    review_levels_tuple = tuple(_probability(level, "review_level") for level in review_levels)
    return {
        "policy": policy_validated,
        "individualFactorEffects": tuple(dict(item) for item in factor_effects),
        "maximumPossibleMultiplier": _maximum_possible_multiplier(
            factor_effects,
            policy_validated,
        ),
        "rows": rows,
        "multiplierDistribution": {
            "count": len(rows),
            "max": max(multipliers) if multipliers else 0.0,
            "min": min(multipliers) if multipliers else 0.0,
        },
        "proportionAboveWarningThreshold": _weighted_fraction(
            rows,
            lambda item: item["combinedMultiplier"]
            >= policy_validated["reviewWarningThreshold"],
            total_weight,
        ),
        "proportionAboveBlockingThreshold": _weighted_fraction(
            rows,
            lambda item: item["combinedMultiplier"]
            >= policy_validated["reviewBlockThreshold"],
            total_weight,
        ),
        "fractionWithProgressionProbabilityAboveReviewLevels": {
            str(level): _weighted_fraction(
                rows,
                lambda item, cutoff=level: item["progressionProbability"] >= cutoff,
                total_weight,
            )
            for level in review_levels_tuple
        },
        "maximumIndividualProgressionProbability": max(
            (row["progressionProbability"] for row in rows),
            default=0.0,
        ),
        "caseConcentration": progression_case_concentration_diagnostics(rows),
        "expectedCases": total_cases,
    }


def progression_case_concentration_diagnostics(
    rows: Iterable[Mapping[str, Any]],
) -> dict[str, Any]:
    row_tuple = tuple(rows)
    total_cases = sum(float(row.get("expectedCases", 0.0)) for row in row_tuple)
    total_weight = sum(float(row.get("weight", 0.0)) for row in row_tuple)

    def share(top_fraction: float) -> float:
        if total_cases <= 0.0 or total_weight <= 0.0:
            return 0.0
        remaining_weight = total_weight * top_fraction
        acc = 0.0
        for row in sorted(row_tuple, key=lambda item: float(item.get("expectedCases", 0.0)), reverse=True):
            if remaining_weight <= 0.0:
                break
            weight = float(row.get("weight", 0.0))
            take = min(weight, remaining_weight)
            fraction = 0.0 if weight <= 0.0 else take / weight
            acc += float(row.get("expectedCases", 0.0)) * fraction
            remaining_weight -= take
        return acc / total_cases

    return {
        "totalExpectedCases": total_cases,
        "totalWeight": total_weight,
        "shareExpectedCasesTop1PercentByWeight": share(0.01),
        "shareExpectedCasesTop5PercentByWeight": share(0.05),
        "shareExpectedCasesTop10PercentByWeight": share(0.10),
        "strataDominatingExpectedCases": sorted(
            (
                {
                    "index": row.get("index"),
                    "state": row.get("state"),
                    "weight": row.get("weight"),
                    "combinedMultiplier": row.get("combinedMultiplier", row.get("multiplier")),
                    "expectedCases": row.get("expectedCases", 0.0),
                }
                for row in row_tuple
            ),
            key=lambda item: float(item.get("expectedCases", 0.0)),
            reverse=True,
        )[:10],
    }


def fixed_ratio_calibration_multiplier_sensitivity(
    row: Mapping[str, Any],
    scenarios: Iterable[Mapping[str, Any]],
    *,
    model_baseline_year: int,
    early_to_remote_ratio: float,
    ratio_source: str,
    ratio_review_status: str,
    ratio_provenance: str,
    ascertainment_probability: float = 1.0,
) -> tuple[dict[str, Any], ...]:
    results = []
    for scenario in scenarios:
        risk_policy = scenario["riskPolicy"]
        base_rows = tuple(scenario["strata"])
        enriched_rows = []
        for source_row in base_rows:
            risk = risk_factor_multiplier_for_row(source_row, risk_policy)
            enriched_rows.append({**source_row, "multiplier": risk["combinedMultiplier"]})
        fit = fit_fixed_ratio_progression_scale(
            row,
            enriched_rows,
            model_baseline_year=model_baseline_year,
            early_to_remote_ratio=early_to_remote_ratio,
            ratio_source=ratio_source,
            ratio_review_status=ratio_review_status,
            ratio_provenance=ratio_provenance,
            ascertainment_probability=ascertainment_probability,
        )
        diagnostics = risk_factor_multiplier_diagnostics(
            enriched_rows,
            risk_policy,
            horizon_years=fit["horizonYears"],
            early_hazard=fit["derivedEarlyHazard"],
            remote_hazard=fit["fittedRemoteHazard"],
        )
        results.append(
            {
                "scenarioId": str(scenario.get("scenarioId", risk_policy["policyId"])),
                "riskPolicy": risk_policy,
                "fittedRemoteHazard": fit["fittedRemoteHazard"],
                "derivedEarlyHazard": fit["derivedEarlyHazard"],
                "maximumIndividualProgressionProbability": diagnostics[
                    "maximumIndividualProgressionProbability"
                ],
                "caseConcentration": diagnostics["caseConcentration"],
                "riskDiagnostics": diagnostics,
                "calibration": fit,
            }
        )
    return tuple(results)


def worked_progression_calibration_examples() -> tuple[dict[str, Any], ...]:
    target = {
        "observationId": "worked-prospective",
        "startYear": 2026,
        "endYear": 2026,
        "observedActiveTBCaseCount": 2,
        "populationDenominator": 1000,
        "personYears": 1000,
        "denominatorType": "census_population",
        "populationScope": "whole_population",
        "caseClassification": "incident_follow_up",
        "observationWindowMeaning": "follow_up_incident",
        "ascertainmentMethod": "combined",
        "activeTBClassification": "all_active_tb",
        "source": "synthetic worked example",
        "reviewStatus": "unreviewed_synthetic",
        "notes": "not a reviewed default",
        "uncertainty": {},
    }
    strata = (
        {"state": STATE_RECENT, "weight": 100.0, "remainingEarlyRiskYears": 4.5},
        {"state": STATE_REMOTE_ONLY, "weight": 200.0, "remainingEarlyRiskYears": 0.0},
        {"state": STATE_UNINFECTED, "weight": 700.0, "remainingEarlyRiskYears": 0.0},
    )
    external_policy = build_progression_calibration_policy(
        policy_id=POLICY_EXTERNAL_HAZARDS,
        early_hazard=0.02,
        remote_hazard=0.002,
        source="synthetic worked example",
        reference_population="synthetic population",
        review_status="unreviewed_synthetic",
    )
    fixed = fit_fixed_ratio_progression_scale(
        target,
        strata,
        model_baseline_year=2026,
        early_to_remote_ratio=10.0,
        ratio_source="synthetic worked example",
        ratio_review_status="unreviewed_synthetic",
        ratio_provenance="not a reviewed default",
        ascertainment_probability=1.0,
    )
    incomplete = fit_fixed_ratio_progression_scale(
        target,
        strata,
        model_baseline_year=2026,
        early_to_remote_ratio=10.0,
        ratio_source="synthetic worked example",
        ratio_review_status="unreviewed_synthetic",
        ratio_provenance="not a reviewed default",
        ascertainment_probability=0.5,
    )
    survival = fit_fixed_ratio_progression_scale(
        target,
        strata,
        model_baseline_year=2026,
        early_to_remote_ratio=10.0,
        ratio_source="synthetic worked example",
        ratio_review_status="unreviewed_synthetic",
        ratio_provenance="not a reviewed default",
        ascertainment_probability=1.0,
        competing_survival_probability=0.8,
    )
    zero_target = fit_fixed_ratio_progression_scale(
        {**target, "observationId": "worked-zero", "observedActiveTBCaseCount": 0},
        strata,
        model_baseline_year=2026,
        early_to_remote_ratio=10.0,
        ratio_source="synthetic worked example",
        ratio_review_status="unreviewed_synthetic",
        ratio_provenance="not a reviewed default",
        ascertainment_probability=1.0,
    )
    multiplier_sensitivity = fixed_ratio_calibration_multiplier_sensitivity(
        target,
        (
            {
                "scenarioId": "no_multipliers",
                "riskPolicy": build_risk_factor_application_policy(policy_id=RISK_POLICY_NONE),
                "strata": strata,
            },
            {
                "scenarioId": "moderate_reviewed_multipliers",
                "riskPolicy": build_risk_factor_application_policy(
                    policy_id=RISK_POLICY_REVIEWED_HAZARD_MULTIPLIERS,
                    reviewed_factor_names=("smoking",),
                ),
                "strata": (
                    {
                        "state": STATE_RECENT,
                        "weight": 100.0,
                        "remainingEarlyRiskYears": 4.5,
                        "smoking": True,
                    },
                    {
                        "state": STATE_REMOTE_ONLY,
                        "weight": 200.0,
                        "remainingEarlyRiskYears": 0.0,
                    },
                    {
                        "state": STATE_UNINFECTED,
                        "weight": 700.0,
                        "remainingEarlyRiskYears": 0.0,
                    },
                ),
            },
        ),
        model_baseline_year=2026,
        early_to_remote_ratio=10.0,
        ratio_source="synthetic worked example",
        ratio_review_status="unreviewed_synthetic",
        ratio_provenance="not a reviewed default",
    )
    validation_policy = build_progression_calibration_policy(
        policy_id=POLICY_VALIDATION_ONLY,
        source="synthetic worked example",
        reference_population="synthetic population",
        review_status="unreviewed_synthetic",
    )
    legacy_rows = (
        {
            "state": STATE_RECENT,
            "weight": 10.0,
            "remainingEarlyRiskYears": 4.5,
            "contact": True,
            "renal": True,
            "diabetes": True,
        },
        {"state": STATE_REMOTE_ONLY, "weight": 990.0, "remainingEarlyRiskYears": 0.0},
    )
    legacy_fit = fixed_ratio_calibration_multiplier_sensitivity(
        target,
        (
            {
                "scenarioId": "legacy_or_as_hazard_diagnostic",
                "riskPolicy": build_risk_factor_application_policy(
                    policy_id=RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC
                ),
                "strata": legacy_rows,
            },
        ),
        model_baseline_year=2026,
        early_to_remote_ratio=10.0,
        ratio_source="synthetic worked example",
        ratio_review_status="unreviewed_synthetic",
        ratio_provenance="not a reviewed default",
    )[0]
    impossible = None
    try:
        fit_fixed_ratio_progression_scale(
            {**target, "observationId": "worked-impossible", "observedActiveTBCaseCount": 500},
            strata,
            model_baseline_year=2026,
            early_to_remote_ratio=10.0,
            ratio_source="synthetic worked example",
            ratio_review_status="unreviewed_synthetic",
            ratio_provenance="not a reviewed default",
            ascertainment_probability=1.0,
        )
    except ProgressionCalibrationError as exc:
        impossible = exc.diagnostics
    return (
        {
            "exampleId": "eligible_external_hazards",
            "result": evaluate_external_hazard_policy(
                target,
                strata,
                external_policy,
                model_baseline_year=2026,
            ),
        },
        {"exampleId": "eligible_fixed_ratio_fit", "result": fixed},
        {
            "exampleId": "no_multipliers_vs_moderate_reviewed_multipliers",
            "result": multiplier_sensitivity,
        },
        {"exampleId": "incomplete_ascertainment", "result": incomplete},
        {"exampleId": "survival_curve", "result": survival},
        {"exampleId": "zero_case_target", "result": zero_target},
        {"exampleId": "legacy_or_as_hazard_diagnostic", "result": legacy_fit},
        {
            "exampleId": "retrospective_validation_only",
            "result": evaluate_validation_only_policy(
                {**target, "startYear": 2023, "endYear": 2024},
                strata,
                validation_policy,
                model_baseline_year=2026,
                early_hazard=0.02,
                remote_hazard=0.002,
            ),
        },
        {
            "exampleId": "baseline_prevalence_validation_only",
            "result": evaluate_validation_only_policy(
                {
                    **target,
                    "caseClassification": "prevalent_baseline",
                    "observationWindowMeaning": "baseline_prevalent",
                },
                strata,
                validation_policy,
                model_baseline_year=2026,
                early_hazard=0.02,
                remote_hazard=0.002,
            ),
        },
        {"exampleId": "impossible_high_target", "result": impossible},
    )


def baseline_active_tb_state_sequence_specification() -> dict[str, Any]:
    return {
        "contractVersion": PROGRESSION_CONTRACT_VERSION,
        "mutuallyExclusiveBaselineSequence": [
            "baseline_prevalent_active_tb",
            STATE_RECENT,
            STATE_REMOTE_ONLY,
            STATE_UNINFECTED,
        ],
        "rule": (
            "A person assigned baseline/prevalent active TB is removed from the "
            "ordinary latent TBI preventive-treatment cascade and is not counted "
            "simultaneously in TBI prevalence totals."
        ),
        "futureIntegrationStates": {
            "baseline_prevalent_active_tb": "Active TB already present at model baseline.",
            "screen_detected_active_tb": "Active TB detected through screening.",
            "incident_active_tb_follow_up": "Active TB developing after baseline.",
            STATE_RECENT: "Recent TBI without baseline active TB.",
            STATE_REMOTE_ONLY: "Remote-only TBI without baseline active TB.",
            STATE_UNINFECTED: "No TBI in the modelled exposure windows.",
        },
    }


def progression_sampling_specification() -> dict[str, Any]:
    return {
        "contractVersion": PROGRESSION_CONTRACT_VERSION,
        "recent": (
            "Invert A(t)=m[lambda_E min(t,r)+lambda_L max(0,t-r)] with "
            "-log(1-q). If the event falls after remaining early risk, add the "
            "late-hazard time after r."
        ),
        "remoteOnly": "Invert A(t)=m lambda_L t.",
        "uninfected": "No event is generated under the current no-new-infection pathway.",
        "rngRequirement": "Any future stochastic implementation must require an explicit RNG or seed.",
    }


def worked_progression_diagnostic_table(
    *,
    early_hazard: float,
    remote_hazard: float,
    recent_window_years: float = 5.0,
    multipliers: Iterable[float] = (1.0, 3.0),
    horizons: Iterable[float] = (1.0, 2.0, 5.0, 10.0, 20.0),
) -> tuple[dict[str, Any], ...]:
    scenarios = (
        ("recent_s_0_5", STATE_RECENT, 0.5),
        ("recent_s_2_5", STATE_RECENT, 2.5),
        ("recent_s_4_9", STATE_RECENT, 4.9),
        ("remote_only", STATE_REMOTE_ONLY, None),
        ("uninfected", STATE_UNINFECTED, None),
    )
    rows = []
    for label, state, time_since in scenarios:
        remaining = (
            0.0
            if time_since is None
            else remaining_early_risk_years(
                time_since,
                recent_window_years=recent_window_years,
            )
        )
        for multiplier in multipliers:
            for horizon in horizons:
                rows.append(
                    {
                        "scenario": label,
                        "state": state,
                        "timeSinceMostRecentInfection": time_since,
                        "remainingEarlyRiskYears": remaining,
                        "multiplier": float(multiplier),
                        "horizonYears": float(horizon),
                        "progressionProbability": progression_probability(
                            state=state,
                            horizon_years=float(horizon),
                            early_hazard=early_hazard,
                            remote_hazard=remote_hazard,
                            multiplier=float(multiplier),
                            remaining_early_risk_years=remaining,
                        ),
                    }
                )
    return tuple(rows)


def synthetic_natural_history_cost_invariance_example(
    *,
    cost_parameter: float,
    intervention_parameter: float,
    early_hazard: float,
    remote_hazard: float,
) -> dict[str, Any]:
    strata = [
        {
            "state": STATE_RECENT,
            "weight": 0.2,
            "multiplier": 1.0,
            "remainingEarlyRiskYears": 4.5,
        },
        {
            "state": STATE_REMOTE_ONLY,
            "weight": 0.3,
            "multiplier": 1.0,
            "remainingEarlyRiskYears": 0.0,
        },
        {
            "state": STATE_UNINFECTED,
            "weight": 0.5,
            "multiplier": 1.0,
            "remainingEarlyRiskYears": 0.0,
        },
    ]
    expected = expected_progression_events(
        strata,
        horizon_years=20.0,
        early_hazard=early_hazard,
        remote_hazard=remote_hazard,
    )
    return {
        "contractVersion": PROGRESSION_CONTRACT_VERSION,
        "costParameterIgnored": float(cost_parameter),
        "interventionParameterIgnored": float(intervention_parameter),
        "expectedProgressionEvents": expected["expectedCases"],
        "note": "Natural-history progression functions do not read costs or intervention parameters.",
    }


def _fixed_ratio_result(
    *,
    requested: float,
    denominator: float,
    horizon: float,
    ascertainment_probability: float,
    ratio: float,
    fitted_remote_hazard: float | None,
    achieved_expected_cases: float,
    convergence_status: str,
    feasibility_status: str,
    warnings: Iterable[str],
    ratio_source: str,
    ratio_review_status: str,
    ratio_provenance: str,
    eligibility: Mapping[str, Any],
    expected_progression: Mapping[str, Any] | None,
    concentration: Mapping[str, Any] | None,
) -> dict[str, Any]:
    remote_hazard = None if fitted_remote_hazard is None else float(fitted_remote_hazard)
    early_hazard = None if remote_hazard is None else ratio * remote_hazard
    achieved = float(achieved_expected_cases)
    return {
        "contractVersion": PROGRESSION_CALIBRATION_CONTRACT_VERSION,
        "policyId": POLICY_FIXED_RATIO_FIT_SCALE,
        "requestedCases": float(requested),
        "denominator": float(denominator),
        "horizonYears": float(horizon),
        "ascertainmentProbability": _probability(
            ascertainment_probability,
            "ascertainment_probability",
        ),
        "suppliedEarlyToRemoteRatio": float(ratio),
        "fittedRemoteHazard": remote_hazard,
        "derivedEarlyHazard": early_hazard,
        "achievedExpectedCases": achieved,
        "residual": abs(achieved - float(requested)),
        "convergenceStatus": str(convergence_status),
        "feasibilityStatus": str(feasibility_status),
        "warnings": list(warnings),
        "provenance": {
            "ratioSource": str(ratio_source),
            "ratioReviewStatus": str(ratio_review_status),
            "ratioProvenance": str(ratio_provenance),
            "note": (
                "This policy fits one scale parameter only and does not use "
                "the inherited earlyLateRatio, MATLAB-v9 value or 10/770 target."
            ),
        },
        "eligibility": dict(eligibility),
        "expectedProgression": None if expected_progression is None else dict(expected_progression),
        "caseConcentration": None if concentration is None else dict(concentration),
    }


def _prepare_progression_strata(
    strata: Iterable[Mapping[str, Any]],
) -> dict[str, Any]:
    included = []
    invalid_states = []
    excluded_baseline = 0.0
    for idx, row in enumerate(strata):
        weight = _row_weight(row, idx)
        state = str(row.get("state", STATE_UNINFECTED))
        if state == STATE_BASELINE_ACTIVE_TB or bool(row.get("baselineActiveTB", False)):
            excluded_baseline += weight
            continue
        try:
            _state_key(state)
        except ValueError:
            invalid_states.append(f"strata[{idx}] has unsupported state {state!r}.")
            continue
        included.append({"index": idx, "weight": weight, "state": state})
    return {
        "includedRows": included,
        "excludedBaselineActiveTBWeight": excluded_baseline,
        "invalidStates": invalid_states,
    }


def _row_weight(row: Mapping[str, Any], idx: int) -> float:
    if "weight" in row:
        value = row["weight"]
    elif "count" in row:
        value = row["count"]
    elif "populationWeight" in row:
        value = row["populationWeight"]
    else:
        value = 1.0
    return _finite_nonnegative_float(value, f"strata[{idx}].weight")


def _survival_probability_at(
    survival_spec: Any,
    horizon: float,
    row: Mapping[str, Any],
    label: str,
) -> float:
    if survival_spec is None:
        if "competingSurvivalProbability" in row:
            return _probability(row["competingSurvivalProbability"], label)
        return 1.0
    if callable(survival_spec):
        return _probability(survival_spec(horizon, row), label)
    if isinstance(survival_spec, Mapping):
        if "survivalProbability" in survival_spec:
            return _probability(survival_spec["survivalProbability"], label)
        if horizon in survival_spec:
            return _probability(survival_spec[horizon], label)
        horizon_key = str(horizon)
        if horizon_key in survival_spec:
            return _probability(survival_spec[horizon_key], label)
        raise ValueError(f"{label} missing horizon {horizon:g}.")
    return _probability(survival_spec, label)


def _observation_model_id(value: Any) -> str:
    model_id = str(value)
    if model_id not in {OBSERVATION_MODEL_BINOMIAL, OBSERVATION_MODEL_POISSON}:
        raise ValueError(
            "observation_model_id must be binomial_person_event_v1 or poisson_count_v1."
        )
    return model_id


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _rank_2d(vectors: Iterable[tuple[float, float]], *, tolerance: float = 1e-10) -> int:
    vector_tuple = tuple(vectors)
    nonzero = [
        vector for vector in vector_tuple if abs(vector[0]) > tolerance or abs(vector[1]) > tolerance
    ]
    if not nonzero:
        return 0
    reference = nonzero[0]
    for vector in nonzero[1:]:
        determinant = reference[0] * vector[1] - reference[1] * vector[0]
        if abs(determinant) > tolerance:
            return 2
    return 1


def _weighted_fraction(
    rows: Iterable[Mapping[str, Any]],
    predicate,
    total_weight: float,
) -> float:
    if total_weight <= 0.0:
        return 0.0
    return (
        sum(float(row.get("weight", 0.0)) for row in rows if predicate(row))
        / total_weight
    )


def _maximum_possible_multiplier(
    factor_effects: Iterable[Mapping[str, Any]],
    policy: Mapping[str, Any],
) -> float:
    policy_id = str(policy["policyId"])
    reviewed_names = set(policy.get("reviewedFactorNames", ()))
    multiplier = 1.0
    for factor in factor_effects:
        name = str(factor["factorName"])
        include = policy_id == RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC or (
            policy_id == RISK_POLICY_REVIEWED_HAZARD_MULTIPLIERS
            and (name in reviewed_names or bool(factor.get("reviewedAsHazardMultiplier", False)))
        )
        if include:
            multiplier *= _positive_float(factor["effect"], f"{name}.effect")
    return multiplier


def _observation_duration_years(observed: Mapping[str, Any]) -> float:
    start_year = int(observed["startYear"])
    end_year = int(observed["endYear"])
    if observed.get("startDate") and observed.get("endDate"):
        # Inclusive date precision is not needed here; retain row-level period
        # while using calendar years for pure expected-value diagnostics.
        return max(float(end_year - start_year + 1), 0.0)
    return max(float(end_year - start_year + 1), 0.0)


def _state_key(state: str) -> str:
    state_key = str(state)
    if state_key not in {STATE_RECENT, STATE_REMOTE_ONLY, STATE_UNINFECTED}:
        raise ValueError("state must be one of recent, remote_only or uninfected.")
    return state_key


def _finite_nonnegative_float(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{label} must be finite and non-negative.")
    return number


def _positive_float(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{label} must be finite and positive.")
    return number


def _probability(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0 or number > 1.0:
        raise ValueError(f"{label} must be a finite probability in [0,1].")
    return number
