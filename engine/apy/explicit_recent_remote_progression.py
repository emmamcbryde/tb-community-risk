from __future__ import annotations

import math
from typing import Any, Iterable, Mapping

from engine.apy.explicit_recent_remote_tbi import (
    STATE_RECENT,
    STATE_REMOTE_ONLY,
    STATE_UNINFECTED,
    validate_active_tb_observation,
)


PROGRESSION_CONTRACT_VERSION = "explicit_recent_remote_tbi_progression_v1"

TARGET_BASELINE_PREVALENCE = "baseline_prevalence_target"
TARGET_SCREEN_DETECTED = "screen_detected_disease_target"
TARGET_PROSPECTIVE_INCIDENT = "prospective_incident_disease_target"
TARGET_RETROSPECTIVE_NOTIFICATION = "retrospective_notification_incidence_target"
TARGET_MIXED_OR_INSUFFICIENT = "mixed_or_insufficiently_defined"

POLICY_BOTH_HAZARDS_SUPPLIED = "both_hazards_externally_supplied"
POLICY_ONE_HAZARD_SUPPLIED = "one_hazard_externally_supplied_fit_other"
POLICY_FIXED_EARLY_LATE_RATIO = "fixed_early_late_ratio_fit_common_scale"
POLICY_JOINT_LIKELIHOOD = "joint_likelihood_multiple_targets"
POLICY_EXTERNAL_VALIDATION_ONLY = "external_validation_only"


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
) -> dict[str, Any]:
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    ascertainment_default = _probability(
        default_ascertainment_probability,
        "default_ascertainment_probability",
    )
    rows = []
    total = 0.0
    total_weight = 0.0
    for idx, row in enumerate(strata):
        weight = _finite_nonnegative_float(row.get("weight", 0.0), f"strata[{idx}].weight")
        state = _state_key(str(row.get("state", STATE_UNINFECTED)))
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
        probability = progression_probability(
            state=state,
            horizon_years=horizon,
            early_hazard=early_hazard,
            remote_hazard=remote_hazard,
            multiplier=multiplier,
            remaining_early_risk_years=remaining,
        )
        expected = weight * probability * ascertainment
        total += expected
        total_weight += weight
        rows.append(
            {
                "state": state,
                "weight": weight,
                "multiplier": multiplier,
                "remainingEarlyRiskYears": remaining,
                "ascertainmentProbability": ascertainment,
                "progressionProbability": probability,
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


def expected_prospective_incident_cases_for_observation(
    row: Mapping[str, Any],
    strata: Iterable[Mapping[str, Any]],
    *,
    model_baseline_year: int,
    early_hazard: float,
    remote_hazard: float,
    ascertainment_probability: float = 1.0,
    competing_survival_function: Any = None,
) -> dict[str, Any]:
    if competing_survival_function is not None:
        raise NotImplementedError(
            "Competing mortality integration is reserved for a future milestone."
        )
    classified = classify_active_tb_observation_target(
        row,
        model_baseline_year=model_baseline_year,
    )
    if classified["targetType"] != TARGET_PROSPECTIVE_INCIDENT:
        raise ValueError(
            "Only genuinely prospective incident targets can be used for this expected-count function."
        )
    observed = classified["observation"]
    horizon = _observation_duration_years(observed)
    expected = expected_progression_events(
        strata,
        horizon_years=horizon,
        early_hazard=early_hazard,
        remote_hazard=remote_hazard,
        default_ascertainment_probability=ascertainment_probability,
    )
    return {
        "contractVersion": PROGRESSION_CONTRACT_VERSION,
        "targetType": classified["targetType"],
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
        "competingMortalityApplied": False,
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
