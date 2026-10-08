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
PROGRESSION_CURVE_CONTRACT_VERSION = (
    "explicit_recent_remote_tbi_progression_curve_v1"
)

CURVE_CENTRAL_TIME_SINCE_INFECTION = (
    "central_piecewise_time_since_infection_progression_v1"
)
CURVE_CONSERVATIVE_TIME_SINCE_INFECTION = (
    "conservative_piecewise_time_since_infection_progression_v1"
)
CURVE_HIGHER_TIME_SINCE_INFECTION = (
    "higher_progression_piecewise_time_since_infection_progression_v1"
)

INTERPOLATION_PIECEWISE_LINEAR_CUMULATIVE_HAZARD = (
    "piecewise_linear_cumulative_hazard_v1"
)
POST_FINAL_ANCHOR_EXPLICIT_HAZARD = "explicit_post_final_anchor_hazard_v1"
POST_FINAL_ANCHOR_FROM_LAST_SEGMENT = (
    "post_final_anchor_hazard_from_10_to_25_year_slope_v1"
)

REINFECTION_POLICY_RESET_CLOCK = "recent_reinfection_resets_progression_clock_v1"
REINFECTION_POLICY_NO_RESET_CLOCK = (
    "recent_reinfection_does_not_reset_progression_clock_v1"
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
RECOMMENDED_PRODUCTION_RISK_POLICY = RISK_POLICY_NONE

STATE_BASELINE_ACTIVE_TB = "baseline_prevalent_active_tb"

COMPETING_MORTALITY_NOT_MODELLED = "not_modelled"
COMPETING_MORTALITY_CONSTANT_HAZARD = "external_constant_mortality_hazard"
COMPETING_MORTALITY_SURVIVAL_CURVE = "external_survival_curve"
COMPETING_MORTALITY_SCALAR_APPROXIMATION = (
    "scalar_horizon_survival_nonproduction_approximation"
)
INTEGRATION_ANALYTIC_PIECEWISE_CONSTANT = "analytic_piecewise_constant_tb_and_death_hazards"
INTEGRATION_NUMERICAL_TRAPEZOID = "numerical_trapezoid_external_survival_curve"

ROOT_TOLERANCE = 1e-12
MAX_ROOT_ITERATIONS = 200
MAX_HAZARD_PER_YEAR = 1e6
DEFAULT_INTEGRATION_TOLERANCE = 1e-10
DEFAULT_NUMERICAL_INTEGRATION_STEPS = 2048
PRACTICAL_IDENTIFIABILITY_CONDITION_NUMBER_THRESHOLD = 1.0e4
PRACTICAL_IDENTIFIABILITY_COMPOSITION_CONTRAST_THRESHOLD = 0.05

DEFAULT_REVIEW_WARNING_MULTIPLIER = 10.0
DEFAULT_REVIEW_BLOCK_MULTIPLIER = 100.0

DEFAULT_RECENT_WINDOW_YEARS = 5.0

CENTRAL_PROGRESSION_ANCHOR_TIMES = (1.0, 2.0, 5.0, 10.0, 25.0)
CENTRAL_PROGRESSION_CUMULATIVE_RISKS = (0.038, 0.050, 0.066, 0.072, 0.079)
CENTRAL_PROGRESSION_SOURCE = (
    "Time-since-infection natural-history synthesis for the United States; "
    "model-derived untreated cumulative risks cited in the branch evidence dossier."
)
CENTRAL_PROGRESSION_SOURCE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC7707158/"

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


def build_progression_curve_contract(
    *,
    curve_identifier: str,
    time_anchors: Iterable[float],
    cumulative_risk_anchors: Iterable[float],
    interpolation_method: str = INTERPOLATION_PIECEWISE_LINEAR_CUMULATIVE_HAZARD,
    post_final_anchor_hazard: float | None = None,
    post_final_anchor_extrapolation_method: str = POST_FINAL_ANCHOR_EXPLICIT_HAZARD,
    recent_window_years: float = DEFAULT_RECENT_WINDOW_YEARS,
    reinfection_policy: str = REINFECTION_POLICY_RESET_CLOCK,
    mortality_policy: str = COMPETING_MORTALITY_NOT_MODELLED,
    risk_factor_progression_policy: str = RISK_POLICY_NONE,
    sources: Iterable[Mapping[str, Any]] = (),
    review_status: str = "",
    notes: str = "",
    anchor_provenance: Iterable[Mapping[str, Any]] = (),
    version_identifier: str = PROGRESSION_CURVE_CONTRACT_VERSION,
) -> dict[str, Any]:
    source_times = tuple(
        _finite_nonnegative_float(value, "time_anchors") for value in time_anchors
    )
    source_risks = tuple(
        _probability_less_than_one(value, "cumulative_risk_anchors")
        for value in cumulative_risk_anchors
    )
    times = (
        (0.0,) + source_times
        if not source_times or source_times[0] != 0.0
        else source_times
    )
    risks = (
        (0.0,) + source_risks
        if not source_times or source_times[0] != 0.0
        else source_risks
    )
    if len(times) != len(risks):
        raise ValueError("time_anchors and cumulative_risk_anchors must have the same length.")
    if len(times) < 2:
        raise ValueError("At least one positive progression anchor is required.")
    for previous, current in zip(times[:-1], times[1:]):
        if current <= previous:
            raise ValueError("time_anchors must be strictly increasing.")
    for previous, current in zip(risks[:-1], risks[1:]):
        if current < previous:
            raise ValueError("cumulative_risk_anchors must be non-decreasing.")
    hazards = tuple(-math.log1p(-risk) for risk in risks)
    segment_hazards = tuple(
        (right_hazard - left_hazard) / (right_time - left_time)
        for left_time, right_time, left_hazard, right_hazard in zip(
            times[:-1],
            times[1:],
            hazards[:-1],
            hazards[1:],
        )
    )
    if post_final_anchor_hazard is None:
        if post_final_anchor_extrapolation_method != POST_FINAL_ANCHOR_FROM_LAST_SEGMENT:
            raise ValueError("post_final_anchor_hazard must be supplied explicitly.")
        post_final_anchor_hazard = segment_hazards[-1]
    payload = {
        "contractVersion": str(version_identifier),
        "curveIdentifier": str(curve_identifier),
        "timeAnchors": list(times),
        "cumulativeRiskAnchors": list(risks),
        "cumulativeHazardAnchors": list(hazards),
        "segmentHazards": list(segment_hazards),
        "interpolationMethod": str(interpolation_method),
        "postFinalAnchorExtrapolationMethod": str(
            post_final_anchor_extrapolation_method
        ),
        "postFinalAnchorHazard": post_final_anchor_hazard,
        "recentWindowYears": recent_window_years,
        "reinfectionPolicy": str(reinfection_policy),
        "mortalityPolicy": str(mortality_policy),
        "riskFactorProgressionPolicy": str(risk_factor_progression_policy),
        "sources": [dict(source) for source in sources],
        "anchorProvenance": [dict(item) for item in anchor_provenance],
        "reviewStatus": str(review_status),
        "notes": str(notes),
        "versionIdentifier": str(version_identifier),
        "transformation": (
            "Cumulative risk anchors F(t) are converted to cumulative hazard "
            "anchors with H(t)=-log(1-F(t)); segment hazards are slopes of "
            "piecewise-linear cumulative hazard and are not directly observed hazards."
        ),
    }
    return validate_progression_curve_contract(payload)


def validate_progression_curve_contract(payload: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("Progression curve contract must be a mapping.")
    contract = str(payload.get("contractVersion", PROGRESSION_CURVE_CONTRACT_VERSION))
    if contract != PROGRESSION_CURVE_CONTRACT_VERSION:
        raise ValueError(
            f"contractVersion must be {PROGRESSION_CURVE_CONTRACT_VERSION!r}."
        )
    curve_identifier = str(payload.get("curveIdentifier") or "")
    if not curve_identifier:
        raise ValueError("curveIdentifier must be supplied.")
    interpolation = str(payload.get("interpolationMethod") or "")
    if interpolation != INTERPOLATION_PIECEWISE_LINEAR_CUMULATIVE_HAZARD:
        raise ValueError("Unsupported progression-curve interpolationMethod.")
    post_method = str(payload.get("postFinalAnchorExtrapolationMethod") or "")
    if post_method not in {
        POST_FINAL_ANCHOR_EXPLICIT_HAZARD,
        POST_FINAL_ANCHOR_FROM_LAST_SEGMENT,
    }:
        raise ValueError("Unsupported postFinalAnchorExtrapolationMethod.")
    times = tuple(
        _finite_nonnegative_float(value, "timeAnchors")
        for value in payload.get("timeAnchors", ())
    )
    risks = tuple(
        _probability_less_than_one(value, "cumulativeRiskAnchors")
        for value in payload.get("cumulativeRiskAnchors", ())
    )
    if len(times) != len(risks):
        raise ValueError("timeAnchors and cumulativeRiskAnchors must have equal length.")
    if len(times) < 2:
        raise ValueError("At least origin and one positive anchor are required.")
    if times[0] != 0.0:
        raise ValueError("Progression curve must include origin time 0.")
    if risks[0] != 0.0:
        raise ValueError("Progression curve origin cumulative risk must be 0.")
    for previous, current in zip(times[:-1], times[1:]):
        if current <= previous:
            raise ValueError("timeAnchors must be strictly increasing.")
    for previous, current in zip(risks[:-1], risks[1:]):
        if current < previous:
            raise ValueError("cumulativeRiskAnchors must be non-decreasing.")
    hazards = tuple(-math.log1p(-risk) for risk in risks)
    for previous, current in zip(hazards[:-1], hazards[1:]):
        if current < previous - 1e-15:
            raise ValueError("cumulativeHazardAnchors must be non-decreasing.")
    segment_hazards = tuple(
        (right_hazard - left_hazard) / (right_time - left_time)
        for left_time, right_time, left_hazard, right_hazard in zip(
            times[:-1],
            times[1:],
            hazards[:-1],
            hazards[1:],
        )
    )
    for value in segment_hazards:
        if value < -1e-15:
            raise ValueError("segment hazards must be non-negative.")
    segment_hazards = tuple(max(0.0, value) for value in segment_hazards)
    supplied_hazards = payload.get("cumulativeHazardAnchors")
    if supplied_hazards not in (None, ""):
        supplied = tuple(float(value) for value in supplied_hazards)
        if len(supplied) != len(hazards):
            raise ValueError("cumulativeHazardAnchors length does not match anchors.")
        for expected, actual in zip(hazards, supplied):
            if abs(expected - actual) > 1e-12:
                raise ValueError("cumulativeHazardAnchors do not match transformed risks.")
    supplied_segments = payload.get("segmentHazards")
    if supplied_segments not in (None, ""):
        supplied = tuple(float(value) for value in supplied_segments)
        if len(supplied) != len(segment_hazards):
            raise ValueError("segmentHazards length does not match anchor intervals.")
        for expected, actual in zip(segment_hazards, supplied):
            if abs(expected - actual) > 1e-12:
                raise ValueError("segmentHazards do not match cumulative hazard slopes.")
    post_hazard = _finite_nonnegative_float(
        payload.get("postFinalAnchorHazard"),
        "postFinalAnchorHazard",
    )
    recent_window = _positive_float(
        payload.get("recentWindowYears", DEFAULT_RECENT_WINDOW_YEARS),
        "recentWindowYears",
    )
    reinfection_policy = str(payload.get("reinfectionPolicy") or "")
    if reinfection_policy not in {
        REINFECTION_POLICY_RESET_CLOCK,
        REINFECTION_POLICY_NO_RESET_CLOCK,
    }:
        raise ValueError("Unsupported reinfectionPolicy.")
    mortality_policy = str(payload.get("mortalityPolicy") or COMPETING_MORTALITY_NOT_MODELLED)
    if mortality_policy not in {
        COMPETING_MORTALITY_NOT_MODELLED,
        COMPETING_MORTALITY_CONSTANT_HAZARD,
        COMPETING_MORTALITY_SURVIVAL_CURVE,
    }:
        raise ValueError("Unsupported mortalityPolicy.")
    risk_policy = str(payload.get("riskFactorProgressionPolicy") or RISK_POLICY_NONE)
    if risk_policy != RISK_POLICY_NONE:
        raise ValueError(
            "Milestone 3A progression curves only support riskFactorProgressionPolicy=none."
        )
    return {
        "contractVersion": contract,
        "curveIdentifier": curve_identifier,
        "timeAnchors": list(times),
        "cumulativeRiskAnchors": list(risks),
        "cumulativeHazardAnchors": list(hazards),
        "segmentHazards": list(segment_hazards),
        "interpolationMethod": interpolation,
        "postFinalAnchorExtrapolationMethod": post_method,
        "postFinalAnchorHazard": post_hazard,
        "recentWindowYears": recent_window,
        "reinfectionPolicy": reinfection_policy,
        "mortalityPolicy": mortality_policy,
        "riskFactorProgressionPolicy": risk_policy,
        "sources": [dict(source) for source in payload.get("sources", ())],
        "anchorProvenance": [
            dict(item) for item in payload.get("anchorProvenance", ())
        ],
        "reviewStatus": str(payload.get("reviewStatus") or ""),
        "notes": str(payload.get("notes") or ""),
        "versionIdentifier": str(
            payload.get("versionIdentifier", PROGRESSION_CURVE_CONTRACT_VERSION)
        ),
        "transformation": str(payload.get("transformation") or ""),
        "validationDiagnostics": {
            "originIsZero": times[0] == 0.0 and risks[0] == 0.0 and hazards[0] == 0.0,
            "cumulativeHazardMonotonic": all(
                right >= left for left, right in zip(hazards[:-1], hazards[1:])
            ),
            "segmentHazardsNonNegative": all(value >= 0.0 for value in segment_hazards),
            "postFinalAnchorHazardExplicit": True,
        },
    }


def progression_curve_json(curve: Mapping[str, Any]) -> str:
    return _canonical_json(validate_progression_curve_contract(curve))


def progression_curve_hash(curve: Mapping[str, Any]) -> str:
    return hashlib.sha256(progression_curve_json(curve).encode("utf-8")).hexdigest()


def build_central_time_since_infection_progression_curve(
    *,
    reinfection_policy: str = REINFECTION_POLICY_RESET_CLOCK,
    mortality_policy: str = COMPETING_MORTALITY_NOT_MODELLED,
) -> dict[str, Any]:
    hazards = tuple(-math.log1p(-risk) for risk in CENTRAL_PROGRESSION_CUMULATIVE_RISKS)
    post_final = (hazards[-1] - hazards[-2]) / (
        CENTRAL_PROGRESSION_ANCHOR_TIMES[-1]
        - CENTRAL_PROGRESSION_ANCHOR_TIMES[-2]
    )
    return build_progression_curve_contract(
        curve_identifier=CURVE_CENTRAL_TIME_SINCE_INFECTION,
        time_anchors=CENTRAL_PROGRESSION_ANCHOR_TIMES,
        cumulative_risk_anchors=CENTRAL_PROGRESSION_CUMULATIVE_RISKS,
        post_final_anchor_hazard=post_final,
        post_final_anchor_extrapolation_method=POST_FINAL_ANCHOR_FROM_LAST_SEGMENT,
        recent_window_years=DEFAULT_RECENT_WINDOW_YEARS,
        reinfection_policy=reinfection_policy,
        mortality_policy=mortality_policy,
        risk_factor_progression_policy=RISK_POLICY_NONE,
        sources=(
            {
                "citation": CENTRAL_PROGRESSION_SOURCE,
                "url": CENTRAL_PROGRESSION_SOURCE_URL,
                "accessDate": "2026-10-08",
                "measure": "model-derived untreated cumulative progression risk",
                "role": "central cumulative-risk anchors",
            },
        ),
        review_status="approved_working_choice_subject_to_evidence_limitations",
        notes=(
            "Central Milestone 3A working curve. Segment hazards are transformed "
            "from cumulative risk anchors; post-25-year hazard is extrapolated "
            "from the 10-25 year cumulative-hazard slope."
        ),
        anchor_provenance=_central_anchor_provenance(),
    )


def build_conservative_time_since_infection_progression_curve() -> dict[str, Any]:
    early = 0.0034
    late = 0.00038
    times = (1.0, 2.0, 5.0, 10.0, 25.0)
    hazards = tuple(early * min(time, 5.0) + late * max(0.0, time - 5.0) for time in times)
    risks = tuple(1.0 - math.exp(-hazard) for hazard in hazards)
    return build_progression_curve_contract(
        curve_identifier=CURVE_CONSERVATIVE_TIME_SINCE_INFECTION,
        time_anchors=times,
        cumulative_risk_anchors=risks,
        post_final_anchor_hazard=late,
        recent_window_years=DEFAULT_RECENT_WINDOW_YEARS,
        reinfection_policy=REINFECTION_POLICY_RESET_CLOCK,
        mortality_policy=COMPETING_MORTALITY_NOT_MODELLED,
        risk_factor_progression_policy=RISK_POLICY_NONE,
        sources=(
            {
                "citation": "Explicit recent/remote TBI parameter decision dossier",
                "url": "docs/explicit_recent_remote_tbi_parameter_decision_dossier.md",
                "accessDate": "2026-10-08",
                "measure": "documented sensitivity segment hazards",
                "role": "conservative sensitivity curve",
            },
        ),
        review_status="sensitivity_only",
        notes=(
            "Conservative sensitivity encoded only from documented dossier "
            "segment hazards: 0.0034/year through five years and 0.00038/year after."
        ),
        anchor_provenance=(
            {
                "timeYears": time,
                "cumulativeRisk": risk,
                "source": "derived from documented conservative segment hazards",
                "transformation": "F(t)=1-exp[-H(t)]",
                "directlyObservedHazard": False,
            }
            for time, risk in zip(times, risks)
        ),
    )


def build_higher_progression_time_since_infection_curve() -> dict[str, Any]:
    first_year = 0.060
    five_year_hazard = -math.log1p(-0.145)
    years_one_to_five = (five_year_hazard - first_year) / 4.0
    late = 0.002
    times = (1.0, 2.0, 5.0, 10.0, 25.0)
    hazards = []
    for time in times:
        if time <= 1.0:
            hazard = first_year * time
        elif time <= 5.0:
            hazard = first_year + years_one_to_five * (time - 1.0)
        else:
            hazard = five_year_hazard + late * (time - 5.0)
        hazards.append(hazard)
    risks = tuple(1.0 - math.exp(-hazard) for hazard in hazards)
    return build_progression_curve_contract(
        curve_identifier=CURVE_HIGHER_TIME_SINCE_INFECTION,
        time_anchors=times,
        cumulative_risk_anchors=risks,
        post_final_anchor_hazard=late,
        recent_window_years=DEFAULT_RECENT_WINDOW_YEARS,
        reinfection_policy=REINFECTION_POLICY_RESET_CLOCK,
        mortality_policy=COMPETING_MORTALITY_NOT_MODELLED,
        risk_factor_progression_policy=RISK_POLICY_NONE,
        sources=(
            {
                "citation": "Explicit recent/remote TBI parameter decision dossier; Trauer five-year sensitivity anchor",
                "url": "docs/explicit_recent_remote_tbi_parameter_decision_dossier.md",
                "accessDate": "2026-10-08",
                "measure": "documented sensitivity segment hazards and five-year cumulative risk",
                "role": "higher-progression sensitivity curve",
            },
        ),
        review_status="sensitivity_only",
        notes=(
            "Higher-progression sensitivity encoded from documented dossier "
            "values: 0.060/year in year 0-1, 14.5% cumulative risk by year "
            "5, and 0.002/year after year 5."
        ),
        anchor_provenance=(
            {
                "timeYears": time,
                "cumulativeRisk": risk,
                "source": "derived from documented higher-progression sensitivity values",
                "transformation": "F(t)=1-exp[-H(t)]",
                "directlyObservedHazard": False,
            }
            for time, risk in zip(times, risks)
        ),
    )


def build_progression_curve_by_identifier(curve_identifier: str) -> dict[str, Any]:
    identifier = str(curve_identifier)
    if identifier == CURVE_CENTRAL_TIME_SINCE_INFECTION:
        return build_central_time_since_infection_progression_curve()
    if identifier == CURVE_CONSERVATIVE_TIME_SINCE_INFECTION:
        return build_conservative_time_since_infection_progression_curve()
    if identifier == CURVE_HIGHER_TIME_SINCE_INFECTION:
        return build_higher_progression_time_since_infection_curve()
    raise ValueError(
        f"Progression curve {identifier!r} is unavailable because complete numerical anchors were not documented."
    )


def progression_curve_availability(curve_identifier: str) -> dict[str, Any]:
    try:
        curve = build_progression_curve_by_identifier(curve_identifier)
    except ValueError as exc:
        return {
            "curveIdentifier": str(curve_identifier),
            "available": False,
            "reason": str(exc),
        }
    return {
        "curveIdentifier": curve["curveIdentifier"],
        "available": True,
        "reason": "Complete numerical anchors are available from the evidence dossier.",
    }


def time_since_curve_cumulative_hazard(
    curve: Mapping[str, Any],
    time_since_infection_years: float,
) -> float:
    validated = validate_progression_curve_contract(curve)
    time_since = _finite_nonnegative_float(
        time_since_infection_years,
        "time_since_infection_years",
    )
    times = tuple(float(value) for value in validated["timeAnchors"])
    hazards = tuple(float(value) for value in validated["cumulativeHazardAnchors"])
    segment_hazards = tuple(float(value) for value in validated["segmentHazards"])
    if time_since == 0.0:
        return 0.0
    for idx, right_time in enumerate(times[1:], start=1):
        left_time = times[idx - 1]
        if time_since <= right_time:
            return hazards[idx - 1] + segment_hazards[idx - 1] * (
                time_since - left_time
            )
    return hazards[-1] + float(validated["postFinalAnchorHazard"]) * (
        time_since - times[-1]
    )


def time_since_curve_cumulative_progression_risk(
    curve: Mapping[str, Any],
    time_since_infection_years: float,
) -> float:
    hazard = time_since_curve_cumulative_hazard(curve, time_since_infection_years)
    return 1.0 - math.exp(-hazard)


def time_since_curve_incremental_cumulative_hazard(
    curve: Mapping[str, Any],
    *,
    time_since_infection_at_baseline: float,
    horizon_years: float,
) -> float:
    time_since = _finite_nonnegative_float(
        time_since_infection_at_baseline,
        "time_since_infection_at_baseline",
    )
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    if horizon == 0.0:
        return 0.0
    start = time_since_curve_cumulative_hazard(curve, time_since)
    end = time_since_curve_cumulative_hazard(curve, time_since + horizon)
    return max(0.0, end - start)


def time_since_curve_conditional_future_progression_probability(
    curve: Mapping[str, Any],
    *,
    time_since_infection_at_baseline: float,
    horizon_years: float,
) -> float:
    increment = time_since_curve_incremental_cumulative_hazard(
        curve,
        time_since_infection_at_baseline=time_since_infection_at_baseline,
        horizon_years=horizon_years,
    )
    return 1.0 - math.exp(-increment)


def time_since_curve_instantaneous_hazard(
    curve: Mapping[str, Any],
    time_since_infection_years: float,
) -> float:
    validated = validate_progression_curve_contract(curve)
    time_since = _finite_nonnegative_float(
        time_since_infection_years,
        "time_since_infection_years",
    )
    times = tuple(float(value) for value in validated["timeAnchors"])
    segment_hazards = tuple(float(value) for value in validated["segmentHazards"])
    for idx, right_time in enumerate(times[1:], start=1):
        left_time = times[idx - 1]
        if left_time <= time_since < right_time:
            return segment_hazards[idx - 1]
    return float(validated["postFinalAnchorHazard"])


def time_since_curve_segment_exposure_times(
    curve: Mapping[str, Any],
    *,
    time_since_infection_at_baseline: float,
    horizon_years: float,
) -> tuple[dict[str, Any], ...]:
    validated = validate_progression_curve_contract(curve)
    start = _finite_nonnegative_float(
        time_since_infection_at_baseline,
        "time_since_infection_at_baseline",
    )
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    end = start + horizon
    if horizon == 0.0:
        return ()
    breakpoints = {start, end}
    for anchor in validated["timeAnchors"]:
        anchor_float = float(anchor)
        if start < anchor_float < end:
            breakpoints.add(anchor_float)
    points = sorted(breakpoints)
    rows = []
    for left, right in zip(points[:-1], points[1:]):
        if right <= left:
            continue
        rows.append(
            {
                "segmentStartSinceInfectionYears": left,
                "segmentEndSinceInfectionYears": right,
                "followUpStartYears": left - start,
                "followUpEndYears": right - start,
                "exposureYears": right - left,
                "segmentHazard": time_since_curve_instantaneous_hazard(curve, left),
                "incrementalCumulativeHazard": time_since_curve_cumulative_hazard(curve, right)
                - time_since_curve_cumulative_hazard(curve, left),
            }
        )
    return tuple(rows)


def time_since_curve_validation_diagnostics(curve: Mapping[str, Any]) -> dict[str, Any]:
    validated = validate_progression_curve_contract(curve)
    return {
        "contractVersion": PROGRESSION_CURVE_CONTRACT_VERSION,
        "curveIdentifier": validated["curveIdentifier"],
        "validationDiagnostics": dict(validated["validationDiagnostics"]),
        "timeAnchors": list(validated["timeAnchors"]),
        "cumulativeRiskAnchors": list(validated["cumulativeRiskAnchors"]),
        "cumulativeHazardAnchors": list(validated["cumulativeHazardAnchors"]),
        "segmentHazards": list(validated["segmentHazards"]),
        "postFinalAnchorHazard": validated["postFinalAnchorHazard"],
        "diagnosticMessages": [
            "Segment hazards are transformed slopes of cumulative hazard, not directly observed hazards.",
            "Post-final-anchor hazard is an explicit extrapolation parameter.",
        ],
    }


def resolve_reinfection_progression_clock(
    *,
    state: str,
    time_since_recent_infection: float | None = None,
    time_since_remote_infection: float | None = None,
    prior_remote_exposure: bool = False,
    reinfection_policy: str = REINFECTION_POLICY_RESET_CLOCK,
    recent_window_years: float = DEFAULT_RECENT_WINDOW_YEARS,
) -> dict[str, Any]:
    state_key = _state_key(state)
    policy = str(reinfection_policy)
    if policy not in {REINFECTION_POLICY_RESET_CLOCK, REINFECTION_POLICY_NO_RESET_CLOCK}:
        raise ValueError("Unsupported reinfection_policy.")
    window = _positive_float(recent_window_years, "recent_window_years")
    if state_key == STATE_UNINFECTED:
        return {
            "state": state_key,
            "progressionClockYears": None,
            "usedClock": "none",
            "priorRemoteExposure": bool(prior_remote_exposure),
            "reinfectionPolicy": policy,
        }
    if state_key == STATE_REMOTE_ONLY:
        if time_since_remote_infection is None:
            raise ValueError("remote_only state requires time_since_remote_infection.")
        remote_time = _finite_nonnegative_float(
            time_since_remote_infection,
            "time_since_remote_infection",
        )
        if remote_time < window:
            raise ValueError("remote_only progression clock must be at least recentWindowYears.")
        return {
            "state": state_key,
            "progressionClockYears": remote_time,
            "usedClock": "remote",
            "priorRemoteExposure": True,
            "reinfectionPolicy": policy,
        }
    if time_since_recent_infection is None:
        raise ValueError("recent state requires time_since_recent_infection.")
    recent_time = _finite_nonnegative_float(
        time_since_recent_infection,
        "time_since_recent_infection",
    )
    if recent_time > window + 1e-12:
        raise ValueError("recent progression clock cannot exceed recentWindowYears.")
    if bool(prior_remote_exposure) and policy == REINFECTION_POLICY_NO_RESET_CLOCK:
        if time_since_remote_infection is None:
            raise ValueError(
                "No-reset reinfection policy requires time_since_remote_infection."
            )
        remote_time = _finite_nonnegative_float(
            time_since_remote_infection,
            "time_since_remote_infection",
        )
        if remote_time < window:
            raise ValueError("remote reinfection clock must be at least recentWindowYears.")
        return {
            "state": state_key,
            "progressionClockYears": remote_time,
            "usedClock": "remote_no_reset_sensitivity",
            "priorRemoteExposure": True,
            "timeSinceRecentInfectionYears": recent_time,
            "timeSinceRemoteInfectionYears": remote_time,
            "reinfectionPolicy": policy,
        }
    return {
        "state": state_key,
        "progressionClockYears": recent_time,
        "usedClock": "recent_reset" if prior_remote_exposure else "recent",
        "priorRemoteExposure": bool(prior_remote_exposure),
        "timeSinceRecentInfectionYears": recent_time,
        "timeSinceRemoteInfectionYears": None
        if time_since_remote_infection is None
        else _finite_nonnegative_float(
            time_since_remote_infection,
            "time_since_remote_infection",
        ),
        "reinfectionPolicy": policy,
    }


def time_since_curve_progression_probability_for_state(
    row: Mapping[str, Any],
    curve: Mapping[str, Any],
    *,
    horizon_years: float,
) -> dict[str, Any]:
    validated = validate_progression_curve_contract(curve)
    state = _state_key(str(row.get("state", STATE_UNINFECTED)))
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    if state == STATE_UNINFECTED:
        return {
            "state": state,
            "timeSinceInfectionAtBaseline": None,
            "incrementalCumulativeHazard": 0.0,
            "conditionalFutureProgressionProbability": 0.0,
            "riskFactorMultiplierApplied": 1.0,
            "riskFactorProgressionPolicy": validated["riskFactorProgressionPolicy"],
        }
    clock = _progression_clock_from_row(row, validated)
    probability = time_since_curve_conditional_future_progression_probability(
        validated,
        time_since_infection_at_baseline=clock["progressionClockYears"],
        horizon_years=horizon,
    )
    increment = time_since_curve_incremental_cumulative_hazard(
        validated,
        time_since_infection_at_baseline=clock["progressionClockYears"],
        horizon_years=horizon,
    )
    return {
        "state": state,
        "timeSinceInfectionAtBaseline": clock["progressionClockYears"],
        "incrementalCumulativeHazard": increment,
        "conditionalFutureProgressionProbability": probability,
        "riskFactorMultiplierApplied": 1.0,
        "riskFactorProgressionPolicy": validated["riskFactorProgressionPolicy"],
        "clockResolution": clock,
    }


def competing_risk_progression_probability_for_curve(
    curve: Mapping[str, Any],
    *,
    time_since_infection_at_baseline: float,
    horizon_years: float,
    competing_survival_curve: Any = None,
    competing_mortality_hazard: Any = None,
    integration_tolerance: float = DEFAULT_INTEGRATION_TOLERANCE,
    integration_steps: int = DEFAULT_NUMERICAL_INTEGRATION_STEPS,
) -> dict[str, Any]:
    validated = validate_progression_curve_contract(curve)
    time_since = _finite_nonnegative_float(
        time_since_infection_at_baseline,
        "time_since_infection_at_baseline",
    )
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    tolerance = _positive_float(integration_tolerance, "integration_tolerance")
    if competing_survival_curve is not None and competing_mortality_hazard is not None:
        raise ValueError(
            "Provide either competing_survival_curve or competing_mortality_hazard, not both."
        )
    if competing_survival_curve is None and competing_mortality_hazard is None:
        probability = time_since_curve_conditional_future_progression_probability(
            validated,
            time_since_infection_at_baseline=time_since,
            horizon_years=horizon,
        )
        return {
            "cumulativeIncidence": probability,
            "competingMortality": COMPETING_MORTALITY_NOT_MODELLED,
            "integrationMethod": "closed_form_no_competing_mortality",
            "integrationTolerance": tolerance,
            "diagnosticMessages": [
                "Competing mortality is not modelled; long-horizon expected cases may be overstated."
            ],
        }
    if competing_mortality_hazard is not None:
        death_hazard = _finite_nonnegative_float(
            competing_mortality_hazard,
            "competing_mortality_hazard",
        )
        probability = _analytic_curve_competing_incidence(
            validated,
            time_since_infection_at_baseline=time_since,
            horizon_years=horizon,
            death_hazard=death_hazard,
        )
        return {
            "cumulativeIncidence": probability,
            "competingMortality": COMPETING_MORTALITY_CONSTANT_HAZARD,
            "competingMortalityHazard": death_hazard,
            "integrationMethod": INTEGRATION_ANALYTIC_PIECEWISE_CONSTANT,
            "integrationTolerance": tolerance,
            "diagnosticMessages": [
                "Competing mortality was integrated with the time-since-infection TB hazard."
            ],
        }
    survival = _validate_survival_curve(
        competing_survival_curve,
        horizon_years=horizon,
        integration_steps=integration_steps,
    )
    probability = _numerical_curve_competing_incidence_with_survival_curve(
        validated,
        time_since_infection_at_baseline=time_since,
        horizon_years=horizon,
        survival_at=survival,
        integration_steps=integration_steps,
    )
    return {
        "cumulativeIncidence": probability,
        "competingMortality": COMPETING_MORTALITY_SURVIVAL_CURVE,
        "integrationMethod": INTEGRATION_NUMERICAL_TRAPEZOID,
        "integrationTolerance": tolerance,
        "integrationSteps": int(integration_steps),
        "diagnosticMessages": [
            "Competing mortality was integrated from an externally supplied survival curve; no horizon-level multiplication was used."
        ],
    }


def expected_progression_events_for_curve(
    strata: Iterable[Mapping[str, Any]],
    curve: Mapping[str, Any],
    *,
    horizon_years: float,
    default_ascertainment_probability: float = 1.0,
    competing_survival_curve: Any = None,
    competing_mortality_hazard: Any = None,
    integration_steps: int = DEFAULT_NUMERICAL_INTEGRATION_STEPS,
    exclude_baseline_active_tb: bool = True,
) -> dict[str, Any]:
    validated = validate_progression_curve_contract(curve)
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    q_default = _probability(
        default_ascertainment_probability,
        "default_ascertainment_probability",
    )
    rows = []
    total = 0.0
    total_weight = 0.0
    included_weight = 0.0
    excluded_baseline = 0.0
    for idx, row in enumerate(strata):
        weight = _row_weight(row, idx)
        total_weight += weight
        raw_state = str(row.get("state", STATE_UNINFECTED))
        if exclude_baseline_active_tb and (
            raw_state == STATE_BASELINE_ACTIVE_TB
            or bool(row.get("baselineActiveTB", False))
        ):
            excluded_baseline += weight
            continue
        state = _state_key(raw_state)
        ascertainment = _probability(
            row.get("ascertainmentProbability", q_default),
            f"strata[{idx}].ascertainmentProbability",
        )
        if state == STATE_UNINFECTED:
            probability = 0.0
            clock = None
            mortality = {
                "competingMortality": COMPETING_MORTALITY_NOT_MODELLED,
                "integrationMethod": "no_tb_risk_uninfected",
            }
        else:
            clock = _progression_clock_from_row(row, validated)
            mortality = competing_risk_progression_probability_for_curve(
                validated,
                time_since_infection_at_baseline=clock["progressionClockYears"],
                horizon_years=horizon,
                competing_survival_curve=competing_survival_curve,
                competing_mortality_hazard=competing_mortality_hazard,
                integration_steps=integration_steps,
            )
            probability = mortality["cumulativeIncidence"]
        expected = weight * probability * ascertainment
        included_weight += weight
        total += expected
        rows.append(
            {
                "index": idx,
                "state": state,
                "weight": weight,
                "timeSinceInfectionAtBaseline": None
                if clock is None
                else clock["progressionClockYears"],
                "priorRemoteExposure": bool(row.get("priorRemoteExposure", False)),
                "ascertainmentProbability": ascertainment,
                "conditionalFutureProgressionProbability": probability,
                "expectedCases": expected,
                "riskFactorMultiplierApplied": 1.0,
                "riskFactorProgressionPolicy": validated["riskFactorProgressionPolicy"],
                "competingMortality": mortality["competingMortality"],
                "integrationMethod": mortality["integrationMethod"],
                "clockResolution": clock,
            }
        )
    return {
        "contractVersion": PROGRESSION_CONTRACT_VERSION,
        "progressionCurveContractVersion": PROGRESSION_CURVE_CONTRACT_VERSION,
        "curveIdentifier": validated["curveIdentifier"],
        "horizonYears": horizon,
        "expectedCases": total,
        "totalWeight": total_weight,
        "includedWeight": included_weight,
        "prospectiveAtRiskPopulation": included_weight,
        "tbiEligiblePopulation": included_weight,
        "baselineActiveTBCount": excluded_baseline,
        "excludedBaselineActiveTBWeight": excluded_baseline,
        "competingMortality": (
            COMPETING_MORTALITY_NOT_MODELLED
            if competing_survival_curve is None and competing_mortality_hazard is None
            else (
                COMPETING_MORTALITY_SURVIVAL_CURVE
                if competing_survival_curve is not None
                else COMPETING_MORTALITY_CONSTANT_HAZARD
            )
        ),
        "riskFactorProgressionPolicy": validated["riskFactorProgressionPolicy"],
        "diagnosticMessages": [
            "Risk-factor progression policy is none; inherited OR multipliers are not applied.",
            "Active-TB observations are validation-only under the Milestone 3A central policy.",
        ],
        "rows": rows,
    }


def validate_active_tb_observation_against_progression_curve(
    row: Mapping[str, Any],
    strata: Iterable[Mapping[str, Any]],
    curve: Mapping[str, Any],
    *,
    model_baseline_year: int,
    ascertainment_probability: float | None = None,
    ascertainment_source: str = "",
    ascertainment_review_status: str = "",
    observation_model_id: str | None = None,
    competing_survival_curve: Any = None,
    competing_mortality_hazard: Any = None,
) -> dict[str, Any]:
    validated_curve = validate_progression_curve_contract(curve)
    strata_tuple = tuple(strata)
    eligibility = assess_observation_prospective_calibration_eligibility(
        row,
        model_baseline_year=model_baseline_year,
        strata=strata_tuple,
        ascertainment_probability=ascertainment_probability,
        ascertainment_source=ascertainment_source,
        ascertainment_review_status=ascertainment_review_status,
    )
    if not eligibility["eligible"]:
        return {
            "contractVersion": PROGRESSION_CONTRACT_VERSION,
            "progressionCurveContractVersion": PROGRESSION_CURVE_CONTRACT_VERSION,
            "curveIdentifier": validated_curve["curveIdentifier"],
            "validationOnly": True,
            "eligibleForProspectiveValidation": False,
            "eligibility": eligibility,
            "diagnosticMessages": [
                "Observation is not eligible for prospective incident validation against the curve."
            ],
        }
    q = _strict_positive_probability(
        ascertainment_probability,
        "ascertainment_probability",
    )
    observed = eligibility["observation"]
    horizon = eligibility["horizonYears"]
    expected = expected_progression_events_for_curve(
        strata_tuple,
        validated_curve,
        horizon_years=horizon,
        default_ascertainment_probability=q,
        competing_survival_curve=competing_survival_curve,
        competing_mortality_hazard=competing_mortality_hazard,
    )
    observed_cases = float(observed["observedActiveTBCaseCount"])
    expected_cases = float(expected["expectedCases"])
    absolute_difference = expected_cases - observed_cases
    relative_difference = (
        math.inf
        if observed_cases == 0.0 and expected_cases > 0.0
        else (0.0 if observed_cases == 0.0 else absolute_difference / observed_cases)
    )
    likelihood = None
    if observation_model_id is not None:
        likelihood = active_tb_observation_log_likelihood(
            observed_cases=observed_cases,
            expected_cases=expected_cases,
            denominator=float(observed["populationDenominator"]),
            observation_model_id=observation_model_id,
        )
    denominators = eligibility["denominatorSummary"]
    return {
        "contractVersion": PROGRESSION_CONTRACT_VERSION,
        "progressionCurveContractVersion": PROGRESSION_CURVE_CONTRACT_VERSION,
        "curveIdentifier": validated_curve["curveIdentifier"],
        "validationOnly": True,
        "fittingPerformed": False,
        "eligibleForProspectiveValidation": True,
        "observationId": observed["observationId"],
        "observationStartYear": observed["startYear"],
        "observationEndYear": observed["endYear"],
        "observationHorizonYears": horizon,
        "observedActiveTBCaseCount": observed_cases,
        "sourcePopulationDenominator": denominators["sourcePopulationDenominator"],
        "populationDenominator": observed["populationDenominator"],
        "prospectiveAtRiskPopulation": denominators["prospectiveAtRiskPopulation"],
        "tbiEligiblePopulation": denominators["tbiEligiblePopulation"],
        "baselineActiveTBCount": denominators["baselineActiveTBCount"],
        "personTimeAtRisk": denominators["personTimeAtRisk"],
        "ascertainmentProbability": q,
        "ascertainmentSource": str(ascertainment_source),
        "ascertainmentReviewStatus": str(ascertainment_review_status),
        "expectedActiveTBCaseCount": expected_cases,
        "absoluteDifferenceExpectedMinusObserved": absolute_difference,
        "relativeDifferenceExpectedMinusObserved": relative_difference,
        "likelihood": likelihood,
        "expectedProgression": expected,
        "eligibility": eligibility,
        "diagnosticMessages": [
            "Validation-only comparison; active-TB observations do not fit or alter the progression curve."
        ],
    }


def worked_time_since_progression_diagnostic_table(
    curve: Mapping[str, Any] | None = None,
    *,
    horizons: Iterable[float] = (1.0, 2.0, 5.0, 10.0, 20.0),
    infection_times: Iterable[float] = (0.0, 0.5, 2.5, 4.9, 5.0, 10.0, 25.0, 50.0),
) -> tuple[dict[str, Any], ...]:
    validated = (
        build_central_time_since_infection_progression_curve()
        if curve is None
        else validate_progression_curve_contract(curve)
    )
    rows = []
    for time_since in infection_times:
        for horizon in horizons:
            conditioned = time_since_curve_conditional_future_progression_probability(
                validated,
                time_since_infection_at_baseline=float(time_since),
                horizon_years=float(horizon),
            )
            unconditioned = time_since_curve_cumulative_progression_risk(
                validated,
                float(time_since) + float(horizon),
            )
            rows.append(
                {
                    "curveIdentifier": validated["curveIdentifier"],
                    "timeSinceInfectionAtBaseline": float(time_since),
                    "horizonYears": float(horizon),
                    "conditionalFutureProgressionProbability": conditioned,
                    "incorrectUnconditionedProgressionProbability": unconditioned,
                    "biasAvoidedByConditioning": unconditioned - conditioned,
                    "incrementalCumulativeHazard": time_since_curve_incremental_cumulative_hazard(
                        validated,
                        time_since_infection_at_baseline=float(time_since),
                        horizon_years=float(horizon),
                    ),
                }
            )
    return tuple(rows)


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


def competing_risk_progression_probability(
    *,
    state: str,
    horizon_years: float,
    early_hazard: float,
    remote_hazard: float,
    multiplier: float = 1.0,
    remaining_early_risk_years: float = 0.0,
    competing_mortality_hazard: Any = None,
    competing_survival_curve: Any = None,
    integration_tolerance: float = DEFAULT_INTEGRATION_TOLERANCE,
    integration_steps: int = DEFAULT_NUMERICAL_INTEGRATION_STEPS,
) -> dict[str, Any]:
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    tolerance = _positive_float(integration_tolerance, "integration_tolerance")
    if competing_mortality_hazard is not None and competing_survival_curve is not None:
        raise ValueError(
            "Provide either competing_mortality_hazard or competing_survival_curve, not both."
        )
    if competing_mortality_hazard is None and competing_survival_curve is None:
        probability = progression_probability(
            state=state,
            horizon_years=horizon,
            early_hazard=early_hazard,
            remote_hazard=remote_hazard,
            multiplier=multiplier,
            remaining_early_risk_years=remaining_early_risk_years,
        )
        return {
            "cumulativeIncidence": probability,
            "competingMortality": COMPETING_MORTALITY_NOT_MODELLED,
            "integrationMethod": "closed_form_no_competing_mortality",
            "integrationTolerance": tolerance,
            "diagnosticMessages": [
                "Competing mortality is not modelled; long-horizon expected cases may be overstated."
            ],
        }
    if competing_mortality_hazard is not None:
        death_hazard = _finite_nonnegative_float(
            competing_mortality_hazard,
            "competing_mortality_hazard",
        )
        probability = _analytic_competing_incidence(
            state=state,
            horizon_years=horizon,
            early_hazard=early_hazard,
            remote_hazard=remote_hazard,
            multiplier=multiplier,
            remaining_early_risk_years=remaining_early_risk_years,
            death_hazard=death_hazard,
        )
        return {
            "cumulativeIncidence": probability,
            "competingMortality": COMPETING_MORTALITY_CONSTANT_HAZARD,
            "competingMortalityHazard": death_hazard,
            "integrationMethod": INTEGRATION_ANALYTIC_PIECEWISE_CONSTANT,
            "integrationTolerance": tolerance,
            "diagnosticMessages": [
                "Competing mortality was integrated as an external cause-specific death hazard."
            ],
        }
    survival = _validate_survival_curve(
        competing_survival_curve,
        horizon_years=horizon,
        integration_steps=integration_steps,
    )
    probability = _numerical_competing_incidence_with_survival_curve(
        state=state,
        horizon_years=horizon,
        early_hazard=early_hazard,
        remote_hazard=remote_hazard,
        multiplier=multiplier,
        remaining_early_risk_years=remaining_early_risk_years,
        survival_at=survival,
        integration_steps=integration_steps,
    )
    return {
        "cumulativeIncidence": probability,
        "competingMortality": COMPETING_MORTALITY_SURVIVAL_CURVE,
        "integrationMethod": INTEGRATION_NUMERICAL_TRAPEZOID,
        "integrationTolerance": tolerance,
        "integrationSteps": int(integration_steps),
        "diagnosticMessages": [
            "Competing mortality was integrated from an externally supplied survival curve."
        ],
    }


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
    competing_survival_curve: Any = None,
    competing_mortality_hazard: Any = None,
    allow_scalar_survival_approximation: bool = False,
    integration_tolerance: float = DEFAULT_INTEGRATION_TOLERANCE,
    integration_steps: int = DEFAULT_NUMERICAL_INTEGRATION_STEPS,
    exclude_baseline_active_tb: bool = True,
) -> dict[str, Any]:
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    ascertainment_default = _probability(
        default_ascertainment_probability,
        "default_ascertainment_probability",
    )
    mortality_spec = _resolve_competing_mortality_inputs(
        competing_survival_probability=competing_survival_probability,
        competing_survival_curve=competing_survival_curve,
        competing_mortality_hazard=competing_mortality_hazard,
        allow_scalar_survival_approximation=allow_scalar_survival_approximation,
    )
    rows = []
    total = 0.0
    total_weight = 0.0
    included_weight = 0.0
    excluded_baseline_weight = 0.0
    diagnostic_messages: list[str] = []
    if mortality_spec["mode"] == COMPETING_MORTALITY_NOT_MODELLED:
        diagnostic_messages.append(
            "Competing mortality is not modelled; long-horizon expected cases may be overstated."
        )
    if mortality_spec["mode"] == COMPETING_MORTALITY_SCALAR_APPROXIMATION:
        diagnostic_messages.append(
            "Scalar horizon survival is a non-production approximation and is not a competing-risk cumulative-incidence calculation."
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
        if mortality_spec["mode"] == COMPETING_MORTALITY_SCALAR_APPROXIMATION:
            base = progression_probability(
                state=state,
                horizon_years=horizon,
                early_hazard=early_hazard,
                remote_hazard=remote_hazard,
                multiplier=multiplier,
                remaining_early_risk_years=remaining,
            )
            competing_survival = _survival_probability_at(
                mortality_spec["scalarSurvival"],
                horizon,
                row,
                f"strata[{idx}].competingSurvivalProbability",
            )
            incidence = base * competing_survival
            mortality_diagnostics = {
                "cumulativeIncidence": incidence,
                "competingMortality": COMPETING_MORTALITY_SCALAR_APPROXIMATION,
                "integrationMethod": COMPETING_MORTALITY_SCALAR_APPROXIMATION,
                "diagnosticMessages": [
                    "Scalar horizon survival approximation used only because allow_scalar_survival_approximation=True."
                ],
            }
        else:
            mortality_diagnostics = competing_risk_progression_probability(
                state=state,
                horizon_years=horizon,
                early_hazard=early_hazard,
                remote_hazard=remote_hazard,
                multiplier=multiplier,
                remaining_early_risk_years=remaining,
                competing_mortality_hazard=mortality_spec["mortalityHazard"],
                competing_survival_curve=mortality_spec["survivalCurve"],
                integration_tolerance=integration_tolerance,
                integration_steps=integration_steps,
            )
            incidence = mortality_diagnostics["cumulativeIncidence"]
        expected = weight * incidence * ascertainment
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
                "progressionProbabilityNoCompetingMortality": progression_probability(
                    state=state,
                    horizon_years=horizon,
                    early_hazard=early_hazard,
                    remote_hazard=remote_hazard,
                    multiplier=multiplier,
                    remaining_early_risk_years=remaining,
                ),
                "competingRiskProgressionProbability": incidence,
                "competingMortality": mortality_diagnostics["competingMortality"],
                "integrationMethod": mortality_diagnostics["integrationMethod"],
                "expectedCases": expected,
            }
        )
    tbi_eligible_weight = included_weight
    return {
        "contractVersion": PROGRESSION_CONTRACT_VERSION,
        "horizonYears": horizon,
        "earlyHazard": _finite_nonnegative_float(early_hazard, "early_hazard"),
        "remoteHazard": _finite_nonnegative_float(remote_hazard, "remote_hazard"),
        "expectedCases": total,
        "totalWeight": total_weight,
        "includedWeight": included_weight,
        "prospectiveAtRiskPopulation": included_weight,
        "tbiEligiblePopulation": tbi_eligible_weight,
        "baselineActiveTBCount": excluded_baseline_weight,
        "excludedBaselineActiveTBWeight": excluded_baseline_weight,
        "competingMortality": mortality_spec["mode"],
        "integrationMethod": mortality_spec["integrationMethod"],
        "integrationTolerance": _positive_float(integration_tolerance, "integration_tolerance"),
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
    ascertainment_source: str = "",
    ascertainment_review_status: str = "",
) -> dict[str, Any]:
    classified = classify_active_tb_observation_target(
        row,
        model_baseline_year=model_baseline_year,
    )
    observed = classified["observation"]
    reasons: list[str] = []
    warnings: list[str] = []
    source_row = row
    strata_tuple = None if strata is None else tuple(strata)
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
        q = _strict_positive_probability(
            ascertainment_probability,
            "ascertainment_probability",
        )
        if not str(ascertainment_source).strip():
            reasons.append("Ascertainment probability requires a source.")
        if not str(ascertainment_review_status).strip():
            reasons.append("Ascertainment probability requires a review status.")
        if q == 1.0 and "complete" not in str(ascertainment_source).lower():
            reasons.append(
                "q=1 requires an explicit complete-ascertainment assumption in the source."
            )
        if q < 1.0:
            warnings.append("Incomplete ascertainment is modelled explicitly through q.")
    if strata_tuple is None:
        reasons.append("Compatible baseline population composition/strata are required.")
        included_weight = 0.0
        excluded_baseline = 0.0
    else:
        prepared = _prepare_progression_strata(strata_tuple)
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
    denominator_summary = _progression_denominator_summary(
        observed,
        strata=() if strata_tuple is None else strata_tuple,
        horizon_years=horizon,
        numerator_includes_baseline_active_tb=source_row.get(
            "numeratorIncludesBaselineActiveTB",
            observed.get("numeratorIncludesBaselineActiveTB"),
        ),
        numerator_includes_prevalent_cases=source_row.get(
            "numeratorIncludesPrevalentCases",
            observed.get("numeratorIncludesPrevalentCases"),
        ),
    )
    if denominator_summary["baselineActiveTBCount"] > 0.0:
        if denominator_summary["numeratorIncludesBaselineActiveTB"] is None:
            reasons.append(
                "Baseline active TB is present but numerator composition is ambiguous."
            )
        elif denominator_summary["numeratorIncludesBaselineActiveTB"]:
            reasons.append(
                "Prospective incident calibration cannot use a numerator that includes baseline active TB."
            )
    if denominator_summary["numeratorIncludesPrevalentCases"] is True:
        reasons.append(
            "Prospective incident calibration cannot use a numerator that includes prevalent cases."
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
        "denominatorSummary": denominator_summary,
        "ascertainment": {
            "probability": None
            if ascertainment_probability is None
            else _strict_positive_probability(
                ascertainment_probability,
                "ascertainment_probability",
            ),
            "source": str(ascertainment_source),
            "reviewStatus": str(ascertainment_review_status),
        },
    }


def expected_prospective_incident_cases_for_observation(
    row: Mapping[str, Any],
    strata: Iterable[Mapping[str, Any]],
    *,
    model_baseline_year: int,
    early_hazard: float,
    remote_hazard: float,
    ascertainment_probability: float | None = None,
    ascertainment_source: str = "",
    ascertainment_review_status: str = "",
    competing_survival_function: Any = None,
    competing_survival_probability: Any = None,
    competing_survival_curve: Any = None,
    competing_mortality_hazard: Any = None,
    allow_scalar_survival_approximation: bool = False,
    integration_tolerance: float = DEFAULT_INTEGRATION_TOLERANCE,
) -> dict[str, Any]:
    supplied_survival_specs = [
        value
        for value in (
            competing_survival_function,
            competing_survival_probability,
            competing_survival_curve,
        )
        if value is not None
    ]
    if len(supplied_survival_specs) > 1:
        raise ValueError(
            "Provide only one survival-curve input."
        )
    survival_spec = (
        competing_survival_function
        if competing_survival_function is not None
        else (
            competing_survival_curve
            if competing_survival_curve is not None
            else competing_survival_probability
        )
    )
    strata_tuple = tuple(strata)
    eligibility = assess_observation_prospective_calibration_eligibility(
        row,
        model_baseline_year=model_baseline_year,
        strata=strata_tuple,
        ascertainment_probability=ascertainment_probability,
        ascertainment_source=ascertainment_source,
        ascertainment_review_status=ascertainment_review_status,
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
        competing_survival_curve=survival_spec,
        competing_mortality_hazard=competing_mortality_hazard,
        allow_scalar_survival_approximation=allow_scalar_survival_approximation,
        integration_tolerance=integration_tolerance,
    )
    denominators = eligibility["denominatorSummary"]
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
        "sourcePopulationDenominator": denominators["sourcePopulationDenominator"],
        "baselineActiveTBCount": denominators["baselineActiveTBCount"],
        "prospectiveAtRiskPopulation": denominators["prospectiveAtRiskPopulation"],
        "tbiEligiblePopulation": denominators["tbiEligiblePopulation"],
        "personTimeAtRisk": denominators["personTimeAtRisk"],
        "denominatorSummary": denominators,
        "personYears": observed["personYears"],
        "ascertainmentProbability": _strict_positive_probability(
            ascertainment_probability,
            "ascertainment_probability",
        ),
        "ascertainmentSource": str(ascertainment_source),
        "ascertainmentReviewStatus": str(ascertainment_review_status),
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
    ascertainment_probability: float | None = None,
    ascertainment_source: str = "",
    ascertainment_review_status: str = "",
    competing_survival_probability: Any = None,
    competing_survival_curve: Any = None,
    competing_mortality_hazard: Any = None,
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
        ascertainment_source=ascertainment_source,
        ascertainment_review_status=ascertainment_review_status,
        competing_survival_probability=competing_survival_probability,
        competing_survival_curve=competing_survival_curve,
        competing_mortality_hazard=competing_mortality_hazard,
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
    ascertainment_probability: float | None = None,
    ascertainment_source: str = "",
    ascertainment_review_status: str = "",
    competing_survival_probability: Any = None,
    competing_survival_curve: Any = None,
    competing_mortality_hazard: Any = None,
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
            ascertainment_source=ascertainment_source,
            ascertainment_review_status=ascertainment_review_status,
            competing_survival_probability=competing_survival_probability,
            competing_survival_curve=competing_survival_curve,
            competing_mortality_hazard=competing_mortality_hazard,
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
    ascertainment_probability: float | None = None,
    ascertainment_source: str = "",
    ascertainment_review_status: str = "",
    competing_survival_probability: Any = None,
    competing_survival_curve: Any = None,
    competing_mortality_hazard: Any = None,
    tolerance: float = ROOT_TOLERANCE,
) -> dict[str, Any]:
    ratio = _positive_float(early_to_remote_ratio, "early_to_remote_ratio")
    if not str(ratio_source).strip():
        raise ValueError("ratio_source must be supplied.")
    if not str(ratio_review_status).strip():
        raise ValueError("ratio_review_status must be supplied.")
    if not str(ratio_provenance).strip():
        raise ValueError("ratio_provenance must be supplied.")
    try:
        ascertainment = _validate_ascertainment_assumption(
            ascertainment_probability,
            ascertainment_source=ascertainment_source,
            ascertainment_review_status=ascertainment_review_status,
        )
    except ValueError as exc:
        diagnostics = {
            "contractVersion": PROGRESSION_CALIBRATION_CONTRACT_VERSION,
            "policyId": POLICY_FIXED_RATIO_FIT_SCALE,
            "feasibilityStatus": "ineligible_ascertainment_identifiability",
            "diagnosticMessages": [str(exc)],
        }
        raise ProgressionCalibrationError(
            "Ascertainment must be fixed externally before fixed-ratio calibration.",
            diagnostics,
        ) from exc
    strata_tuple = tuple(strata)
    eligibility = assess_observation_prospective_calibration_eligibility(
        row,
        model_baseline_year=model_baseline_year,
        strata=strata_tuple,
        ascertainment_probability=ascertainment["probability"],
        ascertainment_source=ascertainment["source"],
        ascertainment_review_status=ascertainment["reviewStatus"],
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
            default_ascertainment_probability=ascertainment["probability"],
            competing_survival_probability=competing_survival_probability,
            competing_survival_curve=competing_survival_curve,
            competing_mortality_hazard=competing_mortality_hazard,
        )["expectedCases"]

    maximum = expected_for_scale(MAX_HAZARD_PER_YEAR)
    warnings = list(eligibility["warnings"])
    if requested > maximum + max(float(tolerance), 1e-12):
        diagnostics = _fixed_ratio_result(
            requested=requested,
            denominator=denominator,
            horizon=horizon,
            ascertainment_probability=ascertainment["probability"],
            ascertainment_source=ascertainment["source"],
            ascertainment_review_status=ascertainment["reviewStatus"],
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
                ascertainment_probability=ascertainment["probability"],
                ascertainment_source=ascertainment["source"],
                ascertainment_review_status=ascertainment["reviewStatus"],
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
        default_ascertainment_probability=ascertainment["probability"],
        competing_survival_probability=competing_survival_probability,
        competing_survival_curve=competing_survival_curve,
        competing_mortality_hazard=competing_mortality_hazard,
    )
    concentration = progression_case_concentration_diagnostics(expected["rows"])
    if concentration["shareExpectedCasesTop1PercentByWeight"] > 0.5:
        warnings.append("Highest 1% of eligible weight contributes more than half of expected cases.")
    return _fixed_ratio_result(
        requested=requested,
        denominator=denominator,
        horizon=horizon,
        ascertainment_probability=ascertainment["probability"],
        ascertainment_source=ascertainment["source"],
        ascertainment_review_status=ascertainment["reviewStatus"],
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
    ascertainment_source: str = "complete ascertainment assumption for structural diagnostic",
    ascertainment_review_status: str = "diagnostic_only",
    early_to_remote_ratio: float | None = None,
    condition_number_review_threshold: float = (
        PRACTICAL_IDENTIFIABILITY_CONDITION_NUMBER_THRESHOLD
    ),
    composition_contrast_review_threshold: float = (
        PRACTICAL_IDENTIFIABILITY_COMPOSITION_CONTRAST_THRESHOLD
    ),
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
    singular_diagnostics: dict[str, Any] = _sensitivity_singular_diagnostics(())
    if strata_groups is not None:
        for target, strata in zip(target_tuple, strata_groups):
            assessment = assess_observation_prospective_calibration_eligibility(
                target,
                model_baseline_year=model_baseline_year,
                strata=strata,
                ascertainment_probability=ascertainment_probability,
                ascertainment_source=ascertainment_source,
                ascertainment_review_status=ascertainment_review_status,
            )
            if assessment["eligible"]:
                sensitivity_vectors.append(
                    target_progression_sensitivity_vector(
                        strata,
                        horizon_years=assessment["horizonYears"],
                        ascertainment_probability=ascertainment_probability,
                    )
                )
        singular_diagnostics = _sensitivity_singular_diagnostics(
            sensitivity_vectors,
            rank_tolerance=1e-10,
        )
        rank = int(singular_diagnostics["rank"])
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
        condition_threshold = _positive_float(
            condition_number_review_threshold,
            "condition_number_review_threshold",
        )
        contrast_threshold = _probability(
            composition_contrast_review_threshold,
            "composition_contrast_review_threshold",
        )
        condition_number = singular_diagnostics["conditionNumber"]
        composition_contrast = singular_diagnostics["compositionContrast"]
        structural = len(eligible) >= 2 and rank >= 2
        practical = (
            structural
            and math.isfinite(condition_number)
            and condition_number <= condition_threshold
            and composition_contrast >= contrast_threshold
        )
        identified = practical
        status = "available" if identified else "unavailable_nonidentifiable"
        if len(eligible) < 2:
            messages.append(
                "One aggregate target cannot identify both lambda_E and lambda_L."
            )
        if len(eligible) >= 2 and rank < 2:
            messages.append(
                "Targets with identical or near-identical recent/remote sensitivity do not add independent information."
            )
        if structural and not practical:
            messages.append(
                "Sensitivity rank is not sufficient; practical identifiability review failed under numerical conditioning/composition rules."
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
        "singularValues": singular_diagnostics["singularValues"],
        "conditionNumber": singular_diagnostics["conditionNumber"],
        "compositionContrast": singular_diagnostics["compositionContrast"],
        "practicalIdentifiability": {
            "conditionNumberReviewThreshold": condition_number_review_threshold,
            "compositionContrastReviewThreshold": composition_contrast_review_threshold,
            "thresholdMeaning": "numerical review rule, not biological evidence",
            "profileLikelihoodStatus": (
                "not_implemented; required before production joint estimation"
            ),
            "parameterBoundaryStatus": (
                "not_evaluated; boundary effects must be reviewed before production"
            ),
        },
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
    if competing_survival_probability is not None:
        raise ValueError(
            "Competing survival is not supported in this linear sensitivity diagnostic; use expected-case competing-risk functions for mortality-aware calculations."
        )
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
        if state == STATE_RECENT:
            early += weight * multiplier * min(horizon, remaining) * ascertainment
            remote += weight * multiplier * max(0.0, horizon - remaining) * ascertainment
        elif state == STATE_REMOTE_ONLY:
            remote += weight * multiplier * horizon * ascertainment
    return (early, remote)


def build_risk_factor_application_policy(
    *,
    policy_id: str,
    reviewed_factor_names: Iterable[str] = (),
    review_warning_threshold: float = DEFAULT_REVIEW_WARNING_MULTIPLIER,
    review_block_threshold: float = DEFAULT_REVIEW_BLOCK_MULTIPLIER,
    production_default: bool = False,
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
    if production_default and policy == RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC:
        raise ValueError(
            "legacy_or_as_hazard_diagnostic_only cannot be selected as a production/default runner policy."
        )
    reviewed = tuple(str(name) for name in reviewed_factor_names)
    return {
        "policyId": policy,
        "reviewedFactorNames": list(reviewed),
        "reviewWarningThreshold": warning,
        "reviewBlockThreshold": blocking,
        "productionDefault": bool(production_default),
        "recommendedProductionPolicy": RECOMMENDED_PRODUCTION_RISK_POLICY,
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
        production_default=bool(policy.get("productionDefault", False)),
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
    ascertainment_probability: float | None = None,
    ascertainment_source: str = "",
    ascertainment_review_status: str = "",
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
            ascertainment_source=ascertainment_source,
            ascertainment_review_status=ascertainment_review_status,
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
        ascertainment_source="synthetic complete ascertainment assumption",
        ascertainment_review_status="unreviewed_synthetic",
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
        ascertainment_source="synthetic incomplete ascertainment assumption",
        ascertainment_review_status="unreviewed_synthetic",
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
        ascertainment_source="synthetic complete ascertainment assumption",
        ascertainment_review_status="unreviewed_synthetic",
        competing_mortality_hazard=0.01,
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
        ascertainment_source="synthetic complete ascertainment assumption",
        ascertainment_review_status="unreviewed_synthetic",
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
        ascertainment_probability=1.0,
        ascertainment_source="synthetic complete ascertainment assumption",
        ascertainment_review_status="unreviewed_synthetic",
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
        ascertainment_probability=1.0,
        ascertainment_source="synthetic complete ascertainment assumption",
        ascertainment_review_status="unreviewed_synthetic",
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
            ascertainment_source="synthetic complete ascertainment assumption",
            ascertainment_review_status="unreviewed_synthetic",
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
                ascertainment_probability=1.0,
                ascertainment_source="synthetic complete ascertainment assumption",
                ascertainment_review_status="unreviewed_synthetic",
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


def ascertainment_progression_scale_confounding_diagnostic(
    *,
    observed_cases: float,
    expected_true_cases_at_scale: float,
    ascertainment_probability: float,
) -> dict[str, Any]:
    observed = _finite_nonnegative_float(observed_cases, "observed_cases")
    expected_true = _positive_float(
        expected_true_cases_at_scale,
        "expected_true_cases_at_scale",
    )
    q = _strict_positive_probability(
        ascertainment_probability,
        "ascertainment_probability",
    )
    required_scale_multiplier = observed / (q * expected_true)
    return {
        "observedCases": observed,
        "expectedTrueCasesAtUnitScale": expected_true,
        "ascertainmentProbability": q,
        "requiredProgressionScaleMultiplier": required_scale_multiplier,
        "identifiabilityStatus": "q_and_progression_scale_are_confounding_without_external_q",
        "diagnosticMessage": (
            "For E[C]=q E[C_true(k)], lower q can be offset by higher progression scale."
        ),
    }


def _central_anchor_provenance() -> tuple[dict[str, Any], ...]:
    return tuple(
        {
            "timeYears": time,
            "cumulativeRisk": risk,
            "source": CENTRAL_PROGRESSION_SOURCE,
            "url": CENTRAL_PROGRESSION_SOURCE_URL,
            "measure": "model-derived untreated cumulative progression risk",
            "transformation": "H(t)=-log(1-F(t))",
            "directlyObservedHazard": False,
            "limitations": (
                "Cumulative risk anchor from evidence synthesis; transformed "
                "segment hazards are not directly observed hazards."
            ),
        }
        for time, risk in zip(
            CENTRAL_PROGRESSION_ANCHOR_TIMES,
            CENTRAL_PROGRESSION_CUMULATIVE_RISKS,
        )
    )


def _progression_clock_from_row(
    row: Mapping[str, Any],
    curve: Mapping[str, Any],
) -> dict[str, Any]:
    state = _state_key(str(row.get("state", STATE_UNINFECTED)))
    if state == STATE_UNINFECTED:
        return resolve_reinfection_progression_clock(
            state=state,
            prior_remote_exposure=bool(row.get("priorRemoteExposure", False)),
            reinfection_policy=str(curve["reinfectionPolicy"]),
            recent_window_years=float(curve["recentWindowYears"]),
        )
    if "progressionClockYears" in row:
        clock = _finite_nonnegative_float(
            row["progressionClockYears"],
            "progressionClockYears",
        )
        if state == STATE_REMOTE_ONLY and clock < float(curve["recentWindowYears"]):
            raise ValueError("remote_only progressionClockYears must be at least recentWindowYears.")
        if state == STATE_RECENT and clock > float(curve["recentWindowYears"]) + 1e-12:
            raise ValueError("recent progressionClockYears cannot exceed recentWindowYears.")
        return {
            "state": state,
            "progressionClockYears": clock,
            "usedClock": "explicit_progression_clock",
            "priorRemoteExposure": bool(row.get("priorRemoteExposure", False)),
            "reinfectionPolicy": str(curve["reinfectionPolicy"]),
        }
    recent_time = row.get(
        "timeSinceRecentInfectionYears",
        row.get("timeSinceMostRecentInfection", row.get("timeSinceInfectionYears")),
    )
    remote_time = row.get(
        "timeSinceRemoteInfectionYears",
        row.get("timeSinceRemoteInfection"),
    )
    if state == STATE_REMOTE_ONLY and remote_time is None:
        remote_time = row.get("timeSinceInfectionYears")
    return resolve_reinfection_progression_clock(
        state=state,
        time_since_recent_infection=None if recent_time is None else float(recent_time),
        time_since_remote_infection=None if remote_time is None else float(remote_time),
        prior_remote_exposure=bool(row.get("priorRemoteExposure", False)),
        reinfection_policy=str(curve["reinfectionPolicy"]),
        recent_window_years=float(curve["recentWindowYears"]),
    )


def _analytic_curve_competing_incidence(
    curve: Mapping[str, Any],
    *,
    time_since_infection_at_baseline: float,
    horizon_years: float,
    death_hazard: float,
) -> float:
    start = _finite_nonnegative_float(
        time_since_infection_at_baseline,
        "time_since_infection_at_baseline",
    )
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    death = _finite_nonnegative_float(death_hazard, "death_hazard")
    if horizon == 0.0:
        return 0.0
    probability = 0.0
    baseline_hazard = time_since_curve_cumulative_hazard(curve, start)
    for segment in time_since_curve_segment_exposure_times(
        curve,
        time_since_infection_at_baseline=start,
        horizon_years=horizon,
    ):
        left_follow = float(segment["followUpStartYears"])
        duration = float(segment["exposureYears"])
        tb_hazard = float(segment["segmentHazard"])
        if duration <= 0.0 or tb_hazard <= 0.0:
            continue
        cumulative_tb_to_left = (
            time_since_curve_cumulative_hazard(
                curve,
                float(segment["segmentStartSinceInfectionYears"]),
            )
            - baseline_hazard
        )
        total_hazard = tb_hazard + death
        survival_to_left = math.exp(-(cumulative_tb_to_left + death * left_follow))
        probability += (
            survival_to_left
            * tb_hazard
            * (-math.expm1(-total_hazard * duration))
            / total_hazard
        )
    return min(max(probability, 0.0), 1.0)


def _numerical_curve_competing_incidence_with_survival_curve(
    curve: Mapping[str, Any],
    *,
    time_since_infection_at_baseline: float,
    horizon_years: float,
    survival_at,
    integration_steps: int,
) -> float:
    start = _finite_nonnegative_float(
        time_since_infection_at_baseline,
        "time_since_infection_at_baseline",
    )
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    if horizon == 0.0:
        return 0.0
    steps = int(integration_steps)
    if steps < 16:
        raise ValueError("integration_steps must be at least 16.")
    breakpoints = {0.0, horizon}
    for anchor in validate_progression_curve_contract(curve)["timeAnchors"]:
        follow_time = float(anchor) - start
        if 0.0 < follow_time < horizon:
            breakpoints.add(follow_time)
    probability = 0.0
    for left, right in zip(sorted(breakpoints)[:-1], sorted(breakpoints)[1:]):
        width = right - left
        if width <= 0.0:
            continue
        local_steps = max(16, int(math.ceil(steps * width / horizon)))
        dt = width / local_steps
        previous = _curve_competing_integrand(
            curve,
            time_since_infection_at_baseline=start,
            follow_up_time=left,
            death_survival=survival_at(left),
        )
        segment = 0.0
        for step in range(1, local_steps + 1):
            follow_time = left + step * dt
            current = _curve_competing_integrand(
                curve,
                time_since_infection_at_baseline=start,
                follow_up_time=follow_time,
                death_survival=survival_at(follow_time),
            )
            segment += 0.5 * (previous + current) * dt
            previous = current
        probability += segment
    return min(max(probability, 0.0), 1.0)


def _curve_competing_integrand(
    curve: Mapping[str, Any],
    *,
    time_since_infection_at_baseline: float,
    follow_up_time: float,
    death_survival: float,
) -> float:
    start = _finite_nonnegative_float(
        time_since_infection_at_baseline,
        "time_since_infection_at_baseline",
    )
    follow = _finite_nonnegative_float(follow_up_time, "follow_up_time")
    death_survival_value = _probability(death_survival, "death_survival")
    incremental_hazard = time_since_curve_incremental_cumulative_hazard(
        curve,
        time_since_infection_at_baseline=start,
        horizon_years=follow,
    )
    tb_survival = math.exp(-incremental_hazard)
    tb_hazard = time_since_curve_instantaneous_hazard(curve, start + follow)
    return tb_survival * death_survival_value * tb_hazard


def _validate_ascertainment_assumption(
    ascertainment_probability: float | None,
    *,
    ascertainment_source: str,
    ascertainment_review_status: str,
) -> dict[str, Any]:
    if ascertainment_probability is None:
        raise ValueError(
            "ascertainment_probability must be fixed externally before Policy B can fit."
        )
    q = _strict_positive_probability(
        ascertainment_probability,
        "ascertainment_probability",
    )
    source = str(ascertainment_source or "").strip()
    review_status = str(ascertainment_review_status or "").strip()
    if not source:
        raise ValueError("ascertainment_source must be supplied.")
    if not review_status:
        raise ValueError("ascertainment_review_status must be supplied.")
    if q == 1.0 and "complete" not in source.lower():
        raise ValueError(
            "q=1 requires an explicit complete-ascertainment assumption in ascertainment_source."
        )
    return {
        "probability": q,
        "source": source,
        "reviewStatus": review_status,
    }


def _progression_denominator_summary(
    observed: Mapping[str, Any],
    *,
    strata: Iterable[Mapping[str, Any]],
    horizon_years: float,
    numerator_includes_baseline_active_tb: Any = None,
    numerator_includes_prevalent_cases: Any = None,
) -> dict[str, Any]:
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    prepared = _prepare_progression_strata(strata)
    prospective_at_risk = sum(row["weight"] for row in prepared["includedRows"])
    baseline_active = float(prepared["excludedBaselineActiveTBWeight"])
    source_denominator = float(observed["populationDenominator"])
    supplied_person_years = observed.get("personYears")
    if supplied_person_years is None:
        person_time = prospective_at_risk * horizon
        person_time_source = "calculated_from_prospective_at_risk_population_and_horizon"
    else:
        person_time = float(supplied_person_years)
        person_time_source = "source_observation_person_years"
    return {
        "sourcePopulationDenominator": source_denominator,
        "sourceObservedActiveTBCaseCount": float(observed["observedActiveTBCaseCount"]),
        "baselineActiveTBCount": baseline_active,
        "prospectiveAtRiskPopulation": prospective_at_risk,
        "personTimeAtRisk": person_time,
        "personTimeAtRiskSource": person_time_source,
        "tbiEligiblePopulation": prospective_at_risk,
        "sourceDenominatorPreserved": True,
        "numeratorIncludesBaselineActiveTB": _optional_bool_from_any(
            numerator_includes_baseline_active_tb,
            "numeratorIncludesBaselineActiveTB",
        ),
        "numeratorIncludesPrevalentCases": _optional_bool_from_any(
            numerator_includes_prevalent_cases,
            "numeratorIncludesPrevalentCases",
        ),
        "denominatorRelationship": (
            "source denominator is preserved separately from modelled at-risk and TBI-eligible counts"
        ),
    }


def _resolve_competing_mortality_inputs(
    *,
    competing_survival_probability: Any,
    competing_survival_curve: Any,
    competing_mortality_hazard: Any,
    allow_scalar_survival_approximation: bool,
) -> dict[str, Any]:
    supplied = [
        value
        for value in (
            competing_survival_probability,
            competing_survival_curve,
            competing_mortality_hazard,
        )
        if value is not None
    ]
    if len(supplied) > 1:
        raise ValueError(
            "Provide only one competing-mortality input: survival probability, survival curve or mortality hazard."
        )
    if competing_mortality_hazard is not None:
        return {
            "mode": COMPETING_MORTALITY_CONSTANT_HAZARD,
            "integrationMethod": INTEGRATION_ANALYTIC_PIECEWISE_CONSTANT,
            "mortalityHazard": competing_mortality_hazard,
            "survivalCurve": None,
        }
    if competing_survival_curve is not None:
        return {
            "mode": COMPETING_MORTALITY_SURVIVAL_CURVE,
            "integrationMethod": INTEGRATION_NUMERICAL_TRAPEZOID,
            "mortalityHazard": None,
            "survivalCurve": competing_survival_curve,
        }
    if competing_survival_probability is not None:
        if callable(competing_survival_probability) or isinstance(
            competing_survival_probability,
            Mapping,
        ):
            return {
                "mode": COMPETING_MORTALITY_SURVIVAL_CURVE,
                "integrationMethod": INTEGRATION_NUMERICAL_TRAPEZOID,
                "mortalityHazard": None,
                "survivalCurve": competing_survival_probability,
            }
        if not allow_scalar_survival_approximation:
            raise ValueError(
                "A scalar horizon survival probability is insufficient for production competing-risk cumulative incidence."
            )
        return {
            "mode": COMPETING_MORTALITY_SCALAR_APPROXIMATION,
            "integrationMethod": COMPETING_MORTALITY_SCALAR_APPROXIMATION,
            "mortalityHazard": None,
            "survivalCurve": None,
            "scalarSurvival": competing_survival_probability,
        }
    return {
        "mode": COMPETING_MORTALITY_NOT_MODELLED,
        "integrationMethod": "closed_form_no_competing_mortality",
        "mortalityHazard": None,
        "survivalCurve": None,
    }


def _analytic_competing_incidence(
    *,
    state: str,
    horizon_years: float,
    early_hazard: float,
    remote_hazard: float,
    multiplier: float,
    remaining_early_risk_years: float,
    death_hazard: float,
) -> float:
    segments = _tb_hazard_segments(
        state=state,
        horizon_years=horizon_years,
        early_hazard=early_hazard,
        remote_hazard=remote_hazard,
        multiplier=multiplier,
        remaining_early_risk_years=remaining_early_risk_years,
    )
    cumulative_tb = 0.0
    probability = 0.0
    for start, end, tb_hazard in segments:
        duration = end - start
        if duration <= 0.0:
            continue
        if tb_hazard <= 0.0:
            cumulative_tb += tb_hazard * duration
            continue
        total_hazard = tb_hazard + death_hazard
        survival_to_start = math.exp(-(cumulative_tb + death_hazard * start))
        if total_hazard == 0.0:
            contribution = survival_to_start * tb_hazard * duration
        else:
            contribution = (
                survival_to_start
                * tb_hazard
                * (-math.expm1(-total_hazard * duration))
                / total_hazard
            )
        probability += contribution
        cumulative_tb += tb_hazard * duration
    return min(max(probability, 0.0), 1.0)


def _numerical_competing_incidence_with_survival_curve(
    *,
    state: str,
    horizon_years: float,
    early_hazard: float,
    remote_hazard: float,
    multiplier: float,
    remaining_early_risk_years: float,
    survival_at,
    integration_steps: int,
) -> float:
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    if horizon == 0.0:
        return 0.0
    steps = int(integration_steps)
    if steps < 16:
        raise ValueError("integration_steps must be at least 16.")
    breakpoints = {0.0, horizon}
    remaining = _finite_nonnegative_float(
        remaining_early_risk_years,
        "remaining_early_risk_years",
    )
    if 0.0 < remaining < horizon:
        breakpoints.add(remaining)
    probability = 0.0
    sorted_points = sorted(breakpoints)
    for left, right in zip(sorted_points[:-1], sorted_points[1:]):
        width = right - left
        if width <= 0.0:
            continue
        local_steps = max(16, int(math.ceil(steps * width / horizon)))
        dt = width / local_steps
        previous = _competing_integrand(
            left,
            state=state,
            early_hazard=early_hazard,
            remote_hazard=remote_hazard,
            multiplier=multiplier,
            remaining_early_risk_years=remaining,
            death_survival=survival_at(left),
        )
        segment = 0.0
        for step in range(1, local_steps + 1):
            t = left + step * dt
            current = _competing_integrand(
                t,
                state=state,
                early_hazard=early_hazard,
                remote_hazard=remote_hazard,
                multiplier=multiplier,
                remaining_early_risk_years=remaining,
                death_survival=survival_at(t),
            )
            segment += 0.5 * (previous + current) * dt
            previous = current
        probability += segment
    return min(max(probability, 0.0), 1.0)


def _competing_integrand(
    time_years: float,
    *,
    state: str,
    early_hazard: float,
    remote_hazard: float,
    multiplier: float,
    remaining_early_risk_years: float,
    death_survival: float,
) -> float:
    tb_survival = progression_survival_probability(
        state=state,
        horizon_years=time_years,
        early_hazard=early_hazard,
        remote_hazard=remote_hazard,
        multiplier=multiplier,
        remaining_early_risk_years=remaining_early_risk_years,
    )
    tb_hazard = progression_piecewise_hazard(
        state=state,
        time_years=time_years,
        early_hazard=early_hazard,
        remote_hazard=remote_hazard,
        multiplier=multiplier,
        remaining_early_risk_years=remaining_early_risk_years,
    )
    return tb_survival * death_survival * tb_hazard


def _tb_hazard_segments(
    *,
    state: str,
    horizon_years: float,
    early_hazard: float,
    remote_hazard: float,
    multiplier: float,
    remaining_early_risk_years: float,
) -> tuple[tuple[float, float, float], ...]:
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    state_key = _state_key(state)
    early = _finite_nonnegative_float(early_hazard, "early_hazard")
    late = _finite_nonnegative_float(remote_hazard, "remote_hazard")
    mult = _finite_nonnegative_float(multiplier, "multiplier")
    remaining = _finite_nonnegative_float(
        remaining_early_risk_years,
        "remaining_early_risk_years",
    )
    if horizon == 0.0 or state_key == STATE_UNINFECTED:
        return ()
    if state_key == STATE_REMOTE_ONLY:
        return ((0.0, horizon, mult * late),)
    early_end = min(horizon, remaining)
    segments = []
    if early_end > 0.0:
        segments.append((0.0, early_end, mult * early))
    if horizon > early_end:
        segments.append((early_end, horizon, mult * late))
    return tuple(segments)


def _validate_survival_curve(
    survival_spec: Any,
    *,
    horizon_years: float,
    integration_steps: int,
):
    horizon = _finite_nonnegative_float(horizon_years, "horizon_years")
    steps = int(integration_steps)
    if steps < 16:
        raise ValueError("integration_steps must be at least 16.")
    if callable(survival_spec):
        times = [horizon * idx / min(steps, 128) for idx in range(min(steps, 128) + 1)]

        def survival_at_callable(t: float) -> float:
            try:
                value = survival_spec(t)
            except TypeError:
                value = survival_spec(t, {})
            return _probability(value, "competing_survival_curve")

        _validate_survival_values(times, [survival_at_callable(t) for t in times])
        return survival_at_callable
    if isinstance(survival_spec, Mapping):
        raw_points = sorted((float(key), value) for key, value in survival_spec.items())
    else:
        raw_points = []
        for item in survival_spec:
            if isinstance(item, Mapping):
                raw_points.append((float(item["timeYears"]), item["survivalProbability"]))
            else:
                time_value, survival_value = item
                raw_points.append((float(time_value), survival_value))
        raw_points.sort()
    if not raw_points:
        raise ValueError("competing_survival_curve must include at least time 0.")
    times = [point[0] for point in raw_points]
    values = [_probability(point[1], "competing_survival_curve") for point in raw_points]
    _validate_survival_values(times, values)
    if times[-1] < horizon - 1e-12:
        raise ValueError("competing_survival_curve must cover the observation horizon.")

    def survival_at_points(t: float) -> float:
        if t <= times[0]:
            return values[0]
        for idx in range(1, len(times)):
            if t <= times[idx]:
                left_t = times[idx - 1]
                right_t = times[idx]
                left_s = values[idx - 1]
                right_s = values[idx]
                if right_t == left_t:
                    return right_s
                fraction = (t - left_t) / (right_t - left_t)
                return left_s + fraction * (right_s - left_s)
        return values[-1]

    return survival_at_points


def _validate_survival_values(times: Iterable[float], values: Iterable[float]) -> None:
    time_tuple = tuple(float(value) for value in times)
    value_tuple = tuple(float(value) for value in values)
    if len(time_tuple) != len(value_tuple):
        raise ValueError("survival times and values must have the same length.")
    if not time_tuple or abs(time_tuple[0]) > 1e-12:
        raise ValueError("competing survival must begin at time 0.")
    if abs(value_tuple[0] - 1.0) > 1e-10:
        raise ValueError("competing survival must begin at one.")
    previous_time = time_tuple[0]
    previous_value = value_tuple[0]
    for time_value, survival_value in zip(time_tuple, value_tuple):
        if not math.isfinite(time_value) or time_value < 0.0:
            raise ValueError("survival times must be finite and non-negative.")
        if time_value < previous_time - 1e-12:
            raise ValueError("survival times must be non-decreasing.")
        if not math.isfinite(survival_value) or survival_value < 0.0 or survival_value > 1.0:
            raise ValueError("survival probabilities must be finite and in [0,1].")
        if survival_value > previous_value + 1e-10:
            raise ValueError("competing survival must be non-increasing.")
        previous_time = time_value
        previous_value = survival_value


def _sensitivity_singular_diagnostics(
    vectors: Iterable[tuple[float, float]],
    *,
    rank_tolerance: float = 1e-10,
) -> dict[str, Any]:
    vector_tuple = tuple(vectors)
    nonzero = [
        vector for vector in vector_tuple if abs(vector[0]) > rank_tolerance or abs(vector[1]) > rank_tolerance
    ]
    if not nonzero:
        return {
            "rank": 0,
            "singularValues": (0.0, 0.0),
            "conditionNumber": math.inf,
            "compositionContrast": 0.0,
        }
    a = sum(vector[0] * vector[0] for vector in nonzero)
    b = sum(vector[0] * vector[1] for vector in nonzero)
    c = sum(vector[1] * vector[1] for vector in nonzero)
    trace = a + c
    disc = max((a - c) * (a - c) + 4.0 * b * b, 0.0)
    lambda_1 = max((trace + math.sqrt(disc)) / 2.0, 0.0)
    lambda_2 = max((trace - math.sqrt(disc)) / 2.0, 0.0)
    s1 = math.sqrt(lambda_1)
    s2 = math.sqrt(lambda_2)
    if s1 <= rank_tolerance:
        rank = 0
    elif s2 <= rank_tolerance:
        rank = 1
    else:
        rank = 2
    condition = math.inf if s2 <= rank_tolerance else s1 / s2
    composition_contrast = 0.0
    for i, left in enumerate(nonzero):
        left_norm = math.hypot(left[0], left[1])
        for right in nonzero[i + 1 :]:
            right_norm = math.hypot(right[0], right[1])
            if left_norm == 0.0 or right_norm == 0.0:
                continue
            determinant = abs(left[0] * right[1] - left[1] * right[0])
            composition_contrast = max(
                composition_contrast,
                determinant / (left_norm * right_norm),
            )
    return {
        "rank": rank,
        "singularValues": (s1, s2),
        "conditionNumber": condition,
        "compositionContrast": composition_contrast,
    }


def _optional_bool_from_any(value: Any, label: str) -> bool | None:
    if value in (None, ""):
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "yes", "1"}:
            return True
        if normalized in {"false", "no", "0"}:
            return False
    raise ValueError(f"{label} must be true, false or omitted.")


def _fixed_ratio_result(
    *,
    requested: float,
    denominator: float,
    horizon: float,
    ascertainment_probability: float,
    ascertainment_source: str,
    ascertainment_review_status: str,
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
        "sourcePopulationDenominator": float(
            eligibility.get("denominatorSummary", {}).get(
                "sourcePopulationDenominator",
                denominator,
            )
        ),
        "baselineActiveTBCount": float(
            eligibility.get("denominatorSummary", {}).get("baselineActiveTBCount", 0.0)
        ),
        "prospectiveAtRiskPopulation": float(
            eligibility.get("denominatorSummary", {}).get(
                "prospectiveAtRiskPopulation",
                eligibility.get("includedPopulationWeight", 0.0),
            )
        ),
        "tbiEligiblePopulation": float(
            eligibility.get("denominatorSummary", {}).get(
                "tbiEligiblePopulation",
                eligibility.get("includedPopulationWeight", 0.0),
            )
        ),
        "personTimeAtRisk": eligibility.get("denominatorSummary", {}).get(
            "personTimeAtRisk"
        ),
        "denominatorSummary": dict(eligibility.get("denominatorSummary", {})),
        "horizonYears": float(horizon),
        "ascertainmentProbability": _probability(
            ascertainment_probability,
            "ascertainment_probability",
        ),
        "ascertainmentSource": str(ascertainment_source),
        "ascertainmentReviewStatus": str(ascertainment_review_status),
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


def _probability_less_than_one(value: Any, label: str) -> float:
    number = _probability(value, label)
    if number >= 1.0:
        raise ValueError(f"{label} must be in [0,1).")
    return number


def _strict_positive_probability(value: Any, label: str) -> float:
    number = _probability(value, label)
    if number <= 0.0:
        raise ValueError(f"{label} must be in (0,1].")
    return number
