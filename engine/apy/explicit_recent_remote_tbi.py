from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime
import math
from typing import Any, Iterable, Mapping

import numpy as np


CALIBRATION_CONTRACT_VERSION = "explicit_recent_remote_tbi_calibration_v1"
ANALYSIS_BASIS = "explicit_recent_remote_tbi_foundation_v1"
NATURAL_HISTORY_SEMANTICS = "explicit_recent_remote_tbi_history_v1"
CALIBRATION_METHOD = "constant_hazard_age_window_recent_remote_v1"
INFECTION_HISTORY_CONTRACT_VERSION = "recent_remote_tbi_history_contract_v1"
ACTIVE_TB_OBSERVATION_SCHEMA_VERSION = "active_tb_observation_targets_v1"

RECENT_WINDOW_YEARS = 5.0
REMOTE_HISTORY_CAP_YEARS = 100.0
AGE_PROPORTION_TOLERANCE = 1e-9
PROBABILITY_TOLERANCE = 1e-10
ROOT_TOLERANCE = 1e-12
MAX_ROOT_ITERATIONS = 200
MAX_HAZARD_PER_YEAR = 1e6

RECENT_LABEL = "Recently infected within 5 years"
REMOTE_ONLY_LABEL = (
    "Remote infection only — infected 5 or more years ago with no infection "
    "in the last 5 years"
)

OBSERVATION_WINDOW_MEANINGS = {
    "baseline_prevalent",
    "screen_detected_prevalent",
    "follow_up_incident",
}
CASE_CLASSIFICATIONS = {
    "prevalent_baseline",
    "screen_detected",
    "incident_follow_up",
}
DENOMINATOR_TYPES = {
    "census_population",
    "eligible_population",
    "screened_population",
    "study_population",
    "person_years",
}
POPULATION_SCOPES = {"whole_population", "screened_subgroup"}
ASCERTAINMENT_METHODS = {
    "passive_notification",
    "active_screening",
    "prevalence_survey",
    "combined",
}
ACTIVE_TB_CLASSIFICATIONS = {"all_active_tb", "pulmonary", "unknown"}


class CalibrationError(ValueError):
    def __init__(self, message: str, assessment: Mapping[str, Any] | None = None) -> None:
        super().__init__(message)
        self.assessment = dict(assessment or {})


@dataclass(frozen=True)
class AgeDistribution:
    ages: tuple[float, ...]
    proportions: tuple[float, ...]
    tolerance: float = AGE_PROPORTION_TOLERANCE

    def as_dict(self) -> dict[str, Any]:
        return {
            "ages": list(self.ages),
            "proportions": list(self.proportions),
            "tolerance": self.tolerance,
        }


@dataclass(frozen=True)
class ExposureDurations:
    ages: tuple[float, ...]
    recent_years: tuple[float, ...]
    remote_years: tuple[float, ...]
    recent_window_years: float = RECENT_WINDOW_YEARS
    remote_history_cap_years: float = REMOTE_HISTORY_CAP_YEARS

    def as_dict(self) -> dict[str, Any]:
        return {
            "ages": list(self.ages),
            "recentYears": list(self.recent_years),
            "remoteYears": list(self.remote_years),
            "recentWindowYears": self.recent_window_years,
            "remoteHistoryCapYears": self.remote_history_cap_years,
        }


@dataclass(frozen=True)
class PopulationPrevalences:
    recent: float
    remote_only: float
    uninfected: float
    total_tbi: float
    prior_remote_exposure: float
    recent_with_prior_remote: float

    def as_dict(self) -> dict[str, float]:
        return {
            "recent": self.recent,
            "remoteOnly": self.remote_only,
            "uninfected": self.uninfected,
            "totalTBI": self.total_tbi,
            "priorRemoteExposure": self.prior_remote_exposure,
            "recentWithPriorRemote": self.recent_with_prior_remote,
        }


@dataclass(frozen=True)
class CalibrationResult:
    requested_recent_prevalence: float
    requested_remote_only_prevalence: float
    fitted_recent_hazard: float
    fitted_remote_hazard: float
    achieved_recent_prevalence: float
    achieved_remote_only_prevalence: float
    achieved_total_tbi_prevalence: float
    achieved_uninfected_prevalence: float
    absolute_residuals: Mapping[str, float]
    convergence_status: str
    feasibility_status: str
    diagnostic_messages: tuple[str, ...]
    max_recent_prevalence: float
    max_remote_only_prevalence_given_recent: float
    calibration_contract_version: str = CALIBRATION_CONTRACT_VERSION
    analysis_basis: str = ANALYSIS_BASIS
    natural_history_semantics: str = NATURAL_HISTORY_SEMANTICS
    calibration_method: str = CALIBRATION_METHOD
    infection_history_contract_version: str = INFECTION_HISTORY_CONTRACT_VERSION

    def as_dict(self) -> dict[str, Any]:
        return {
            "calibrationContractVersion": self.calibration_contract_version,
            "analysisBasis": self.analysis_basis,
            "naturalHistorySemantics": self.natural_history_semantics,
            "calibrationMethod": self.calibration_method,
            "infectionHistoryContractVersion": self.infection_history_contract_version,
            "targetDefinitions": {
                "recent": RECENT_LABEL,
                "remoteOnly": REMOTE_ONLY_LABEL,
            },
            "requestedRecentPrevalence": self.requested_recent_prevalence,
            "requestedRemoteOnlyPrevalence": self.requested_remote_only_prevalence,
            "fittedRecentHazard": self.fitted_recent_hazard,
            "fittedRemoteHazard": self.fitted_remote_hazard,
            "achievedRecentPrevalence": self.achieved_recent_prevalence,
            "achievedRemoteOnlyPrevalence": self.achieved_remote_only_prevalence,
            "achievedTotalTBIPrevalence": self.achieved_total_tbi_prevalence,
            "achievedUninfectedPrevalence": self.achieved_uninfected_prevalence,
            "absoluteResiduals": dict(self.absolute_residuals),
            "convergenceStatus": self.convergence_status,
            "feasibilityStatus": self.feasibility_status,
            "diagnosticMessages": list(self.diagnostic_messages),
            "maxRecentPrevalence": self.max_recent_prevalence,
            "maxRemoteOnlyPrevalenceGivenRecent": self.max_remote_only_prevalence_given_recent,
        }


def calibration_result_from_dict(payload: Mapping[str, Any]) -> CalibrationResult:
    return CalibrationResult(
        requested_recent_prevalence=float(payload["requestedRecentPrevalence"]),
        requested_remote_only_prevalence=float(payload["requestedRemoteOnlyPrevalence"]),
        fitted_recent_hazard=float(payload["fittedRecentHazard"]),
        fitted_remote_hazard=float(payload["fittedRemoteHazard"]),
        achieved_recent_prevalence=float(payload["achievedRecentPrevalence"]),
        achieved_remote_only_prevalence=float(payload["achievedRemoteOnlyPrevalence"]),
        achieved_total_tbi_prevalence=float(payload["achievedTotalTBIPrevalence"]),
        achieved_uninfected_prevalence=float(payload["achievedUninfectedPrevalence"]),
        absolute_residuals={
            str(key): float(value)
            for key, value in dict(payload["absoluteResiduals"]).items()
        },
        convergence_status=str(payload["convergenceStatus"]),
        feasibility_status=str(payload["feasibilityStatus"]),
        diagnostic_messages=tuple(str(item) for item in payload["diagnosticMessages"]),
        max_recent_prevalence=float(payload["maxRecentPrevalence"]),
        max_remote_only_prevalence_given_recent=float(
            payload["maxRemoteOnlyPrevalenceGivenRecent"]
        ),
        calibration_contract_version=str(
            payload.get("calibrationContractVersion", CALIBRATION_CONTRACT_VERSION)
        ),
        analysis_basis=str(payload.get("analysisBasis", ANALYSIS_BASIS)),
        natural_history_semantics=str(
            payload.get("naturalHistorySemantics", NATURAL_HISTORY_SEMANTICS)
        ),
        calibration_method=str(payload.get("calibrationMethod", CALIBRATION_METHOD)),
        infection_history_contract_version=str(
            payload.get(
                "infectionHistoryContractVersion",
                INFECTION_HISTORY_CONTRACT_VERSION,
            )
        ),
    )


def exposure_durations_for_ages(
    ages: Iterable[float],
    *,
    recent_window_years: float = RECENT_WINDOW_YEARS,
    remote_history_cap_years: float = REMOTE_HISTORY_CAP_YEARS,
) -> ExposureDurations:
    age_tuple = _finite_nonnegative_tuple(ages, "ages")
    recent_window = _positive_float(recent_window_years, "recent_window_years")
    remote_cap = _positive_float(remote_history_cap_years, "remote_history_cap_years")
    if remote_cap < recent_window:
        raise ValueError("remote_history_cap_years must be at least recent_window_years.")
    recent = tuple(min(recent_window, age) for age in age_tuple)
    remote = tuple(max(0.0, min(remote_cap, age) - recent_window) for age in age_tuple)
    return ExposureDurations(
        ages=age_tuple,
        recent_years=recent,
        remote_years=remote,
        recent_window_years=recent_window,
        remote_history_cap_years=remote_cap,
    )


def cumulative_hazards_for_ages(
    ages: Iterable[float],
    recent_hazard: float,
    remote_hazard: float,
    *,
    recent_window_years: float = RECENT_WINDOW_YEARS,
    remote_history_cap_years: float = REMOTE_HISTORY_CAP_YEARS,
) -> dict[str, tuple[float, ...]]:
    recent_rate = _finite_nonnegative_float(recent_hazard, "recent_hazard")
    remote_rate = _finite_nonnegative_float(remote_hazard, "remote_hazard")
    durations = exposure_durations_for_ages(
        ages,
        recent_window_years=recent_window_years,
        remote_history_cap_years=remote_history_cap_years,
    )
    return {
        "recent": tuple(recent_rate * duration for duration in durations.recent_years),
        "remote": tuple(remote_rate * duration for duration in durations.remote_years),
    }


def mutually_exclusive_state_probabilities(
    ages: Iterable[float],
    recent_hazard: float,
    remote_hazard: float,
    *,
    recent_window_years: float = RECENT_WINDOW_YEARS,
    remote_history_cap_years: float = REMOTE_HISTORY_CAP_YEARS,
    tolerance: float = PROBABILITY_TOLERANCE,
) -> dict[str, tuple[float, ...]]:
    hazards = cumulative_hazards_for_ages(
        ages,
        recent_hazard,
        remote_hazard,
        recent_window_years=recent_window_years,
        remote_history_cap_years=remote_history_cap_years,
    )
    h_recent = np.asarray(hazards["recent"], dtype=float)
    h_remote = np.asarray(hazards["remote"], dtype=float)
    exp_recent = np.exp(-h_recent)
    exp_remote = np.exp(-h_remote)
    recent = 1.0 - exp_recent
    prior_remote_exposure = 1.0 - exp_remote
    remote_only = exp_recent * prior_remote_exposure
    uninfected = exp_recent * exp_remote
    recent_with_prior_remote = recent * prior_remote_exposure
    recent_without_prior_remote = recent * exp_remote
    total = recent + remote_only + uninfected
    _validate_probability_vector(recent, "recent", tolerance)
    _validate_probability_vector(remote_only, "remote_only", tolerance)
    _validate_probability_vector(uninfected, "uninfected", tolerance)
    if not np.allclose(total, 1.0, atol=tolerance, rtol=0.0):
        raise ValueError("State probabilities must sum to one within tolerance.")
    return {
        "recent": _tuple(recent),
        "remote_only": _tuple(remote_only),
        "uninfected": _tuple(uninfected),
        "total_tbi": _tuple(recent + remote_only),
        "prior_remote_exposure": _tuple(prior_remote_exposure),
        "recent_with_prior_remote": _tuple(recent_with_prior_remote),
        "recent_without_prior_remote": _tuple(recent_without_prior_remote),
    }


def validate_age_distribution(
    ages: Iterable[float],
    proportions: Iterable[float],
    *,
    tolerance: float = AGE_PROPORTION_TOLERANCE,
) -> AgeDistribution:
    age_tuple = _finite_nonnegative_tuple(ages, "ages")
    prop_tuple = _finite_nonnegative_tuple(proportions, "proportions")
    if len(age_tuple) == 0:
        raise ValueError("Age distribution must contain at least one age stratum.")
    if len(age_tuple) != len(prop_tuple):
        raise ValueError("Age and proportion arrays must have the same length.")
    total = sum(prop_tuple)
    if abs(total - 1.0) > float(tolerance):
        raise ValueError(
            f"Age proportions must sum to one within tolerance {tolerance:g}; "
            f"received {total:.17g}."
        )
    return AgeDistribution(ages=age_tuple, proportions=prop_tuple, tolerance=tolerance)


def population_weighted_prevalences(
    ages: Iterable[float],
    proportions: Iterable[float],
    recent_hazard: float,
    remote_hazard: float,
) -> PopulationPrevalences:
    distribution = validate_age_distribution(ages, proportions)
    probabilities = mutually_exclusive_state_probabilities(
        distribution.ages,
        recent_hazard,
        remote_hazard,
    )
    weights = np.asarray(distribution.proportions, dtype=float)

    def weighted(name: str) -> float:
        return float(np.sum(weights * np.asarray(probabilities[name], dtype=float)))

    recent = weighted("recent")
    remote_only = weighted("remote_only")
    uninfected = weighted("uninfected")
    return PopulationPrevalences(
        recent=recent,
        remote_only=remote_only,
        uninfected=uninfected,
        total_tbi=recent + remote_only,
        prior_remote_exposure=weighted("prior_remote_exposure"),
        recent_with_prior_remote=weighted("recent_with_prior_remote"),
    )


def assess_calibration_feasibility(
    ages: Iterable[float],
    proportions: Iterable[float],
    requested_recent_prevalence: float,
    requested_remote_only_prevalence: float,
) -> dict[str, Any]:
    distribution = validate_age_distribution(ages, proportions)
    target_recent = _probability(requested_recent_prevalence, "requested_recent_prevalence")
    target_remote = _probability(
        requested_remote_only_prevalence,
        "requested_remote_only_prevalence",
    )
    durations = exposure_durations_for_ages(distribution.ages)
    weights = np.asarray(distribution.proportions, dtype=float)
    recent_durations = np.asarray(durations.recent_years, dtype=float)
    remote_durations = np.asarray(durations.remote_years, dtype=float)
    max_recent = float(np.sum(weights[recent_durations > 0.0]))
    messages: list[str] = []
    feasible = True
    recent_hazard_for_remote_assessment = 0.0

    if target_recent + target_remote > 1.0 + PROBABILITY_TOLERANCE:
        feasible = False
        messages.append("Requested recent and remote-only prevalences exceed one.")
    if target_recent > max_recent + PROBABILITY_TOLERANCE:
        feasible = False
        messages.append(
            "Requested recent prevalence exceeds the maximum possible value "
            f"{max_recent:.12g} for this age distribution."
        )
    if feasible and target_recent > 0.0:
        recent_hazard_for_remote_assessment = _solve_hazard_for_target(
            recent_durations,
            weights,
            target_recent,
            max_recent,
            "recent",
        )
    no_recent = np.exp(-recent_hazard_for_remote_assessment * recent_durations)
    max_remote = float(np.sum(weights[remote_durations > 0.0] * no_recent[remote_durations > 0.0]))
    if target_remote > max_remote + PROBABILITY_TOLERANCE:
        feasible = False
        messages.append(
            "Requested remote-only prevalence exceeds the maximum possible value "
            f"{max_remote:.12g} after accounting for no recent infection."
        )
    return {
        "isFeasible": feasible,
        "feasibilityStatus": "feasible" if feasible else "infeasible",
        "requestedRecentPrevalence": target_recent,
        "requestedRemoteOnlyPrevalence": target_remote,
        "maxRecentPrevalence": max_recent,
        "maxRemoteOnlyPrevalenceGivenRecent": max_remote,
        "diagnosticMessages": messages,
        "calibrationContractVersion": CALIBRATION_CONTRACT_VERSION,
    }


def calibrate_recent_remote_hazards(
    ages: Iterable[float],
    proportions: Iterable[float],
    requested_recent_prevalence: float,
    requested_remote_only_prevalence: float,
    *,
    tolerance: float = ROOT_TOLERANCE,
) -> CalibrationResult:
    distribution = validate_age_distribution(ages, proportions)
    target_recent = _probability(requested_recent_prevalence, "requested_recent_prevalence")
    target_remote = _probability(
        requested_remote_only_prevalence,
        "requested_remote_only_prevalence",
    )
    assessment = assess_calibration_feasibility(
        distribution.ages,
        distribution.proportions,
        target_recent,
        target_remote,
    )
    if not assessment["isFeasible"]:
        raise CalibrationError(
            "Explicit recent/remote TBI calibration target is infeasible: "
            + "; ".join(assessment["diagnosticMessages"]),
            assessment,
        )

    durations = exposure_durations_for_ages(distribution.ages)
    weights = np.asarray(distribution.proportions, dtype=float)
    recent_durations = np.asarray(durations.recent_years, dtype=float)
    remote_durations = np.asarray(durations.remote_years, dtype=float)
    max_recent = float(assessment["maxRecentPrevalence"])
    messages: list[str] = []

    if target_recent <= tolerance:
        recent_hazard = 0.0
        messages.append("Requested recent prevalence is zero; fitted recent hazard set to zero.")
    else:
        recent_hazard = _solve_hazard_for_target(
            recent_durations,
            weights,
            target_recent,
            max_recent,
            "recent",
            tolerance=tolerance,
        )

    no_recent = np.exp(-recent_hazard * recent_durations)
    max_remote = float(np.sum(weights[remote_durations > 0.0] * no_recent[remote_durations > 0.0]))
    if target_remote > max_remote + PROBABILITY_TOLERANCE:
        raise CalibrationError(
            "Explicit recent/remote TBI calibration target is infeasible: "
            "requested remote-only prevalence exceeds the maximum possible value "
            f"{max_remote:.12g} after fitting recent infection.",
            {**assessment, "maxRemoteOnlyPrevalenceGivenRecent": max_remote},
        )
    if target_remote <= tolerance:
        remote_hazard = 0.0
        messages.append("Requested remote-only prevalence is zero; fitted remote hazard set to zero.")
    else:
        remote_hazard = _solve_remote_hazard_for_target(
            remote_durations,
            weights,
            no_recent,
            target_remote,
            max_remote,
            tolerance=tolerance,
        )

    achieved = population_weighted_prevalences(
        distribution.ages,
        distribution.proportions,
        recent_hazard,
        remote_hazard,
    )
    residuals = {
        "recent": abs(achieved.recent - target_recent),
        "remoteOnly": abs(achieved.remote_only - target_remote),
    }
    if any(value > max(10.0 * tolerance, PROBABILITY_TOLERANCE) for value in residuals.values()):
        messages.append("Calibration converged numerically but residuals exceed tolerance.")
        convergence = "residual_tolerance_not_met"
    else:
        convergence = "converged"
    return CalibrationResult(
        requested_recent_prevalence=target_recent,
        requested_remote_only_prevalence=target_remote,
        fitted_recent_hazard=recent_hazard,
        fitted_remote_hazard=remote_hazard,
        achieved_recent_prevalence=achieved.recent,
        achieved_remote_only_prevalence=achieved.remote_only,
        achieved_total_tbi_prevalence=achieved.total_tbi,
        achieved_uninfected_prevalence=achieved.uninfected,
        absolute_residuals=residuals,
        convergence_status=convergence,
        feasibility_status="feasible",
        diagnostic_messages=tuple(messages),
        max_recent_prevalence=max_recent,
        max_remote_only_prevalence_given_recent=max_remote,
    )


def infection_time_quantile(
    *,
    window_start_years_before_baseline: float,
    window_duration_years: float,
    hazard_per_year: float,
    quantile: float,
) -> float:
    start = _finite_nonnegative_float(
        window_start_years_before_baseline,
        "window_start_years_before_baseline",
    )
    duration = _positive_float(window_duration_years, "window_duration_years")
    hazard = _positive_float(hazard_per_year, "hazard_per_year")
    q = _probability(quantile, "quantile")
    event_probability = -math.expm1(-hazard * duration)
    if event_probability <= 0:
        raise ValueError("Cannot condition on an event with zero event probability.")
    years_after_start = -math.log1p(-q * event_probability) / hazard
    return start + min(max(years_after_start, 0.0), duration)


def infection_timing_specification() -> dict[str, Any]:
    return {
        "contractVersion": INFECTION_HISTORY_CONTRACT_VERSION,
        "effectiveStateRule": (
            "The most recent infection event determines the effective progression "
            "state. A recent infection supersedes earlier remote exposure."
        ),
        "recentWindow": {
            "state": "recent",
            "lookbackYears": "[0, min(5, age)]",
            "conditionalTiming": (
                "Draw the most recent event in the recent window from the "
                "constant-hazard distribution truncated to the available recent "
                "window."
            ),
        },
        "remoteOnlyWindow": {
            "state": "remote_only",
            "lookbackYears": "[5, min(100, age)]",
            "conditionalTiming": (
                "Conditional on at least one remote-window event and no recent "
                "event, draw the most recent remote-window event from the "
                "constant-hazard distribution truncated to the remote window."
            ),
        },
        "bothRemoteAndRecent": {
            "effectiveState": "recent",
            "retainedInformation": "priorRemoteExposure=true",
            "progressionImplication": (
                "Recent reinfection resets the higher-progression-risk state in "
                "future integration unless a later scientific decision chooses a "
                "different mechanism."
            ),
        },
        "unresolvedIntegrationQuestion": (
            "The inherited runner uses a Markov recent-to-remote progression "
            "compartment and does not currently distinguish first infection, most "
            "recent infection, and reinfection reset mechanisms."
        ),
    }


def validate_active_tb_observation(row: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(row, Mapping):
        raise ValueError("Active-TB observation must be a mapping.")
    observed_count = _finite_nonnegative_float(
        _required(row, "observedActiveTBCaseCount"),
        "observedActiveTBCaseCount",
    )
    population_denominator = _positive_float(
        _required(row, "populationDenominator"),
        "populationDenominator",
    )
    person_years_value = row.get("personYears")
    person_years = None
    if person_years_value not in (None, ""):
        person_years = _positive_float(person_years_value, "personYears")
    start = _parse_observation_endpoint(row, "start")
    end = _parse_observation_endpoint(row, "end")
    if _endpoint_sort_key(start) > _endpoint_sort_key(end):
        raise ValueError("Observation start must be on or before observation end.")
    window_meaning = _one_of(row, "observationWindowMeaning", OBSERVATION_WINDOW_MEANINGS)
    case_classification = _one_of(row, "caseClassification", CASE_CLASSIFICATIONS)
    denominator_type = _one_of(row, "denominatorType", DENOMINATOR_TYPES)
    population_scope = _one_of(row, "populationScope", POPULATION_SCOPES)
    ascertainment = _one_of(row, "ascertainmentMethod", ASCERTAINMENT_METHODS)
    tb_classification = _one_of(row, "activeTBClassification", ACTIVE_TB_CLASSIFICATIONS)
    uncertainty = row.get("uncertainty", {})
    if uncertainty in (None, ""):
        uncertainty = {}
    if not isinstance(uncertainty, Mapping):
        raise ValueError("uncertainty must be a mapping when supplied.")
    out = {
        "schemaVersion": ACTIVE_TB_OBSERVATION_SCHEMA_VERSION,
        "observationId": str(row.get("observationId") or ""),
        "startDate": start.get("date"),
        "startYear": start["year"],
        "endDate": end.get("date"),
        "endYear": end["year"],
        "observedActiveTBCaseCount": observed_count,
        "populationDenominator": population_denominator,
        "personYears": person_years,
        "denominatorType": denominator_type,
        "populationScope": population_scope,
        "caseClassification": case_classification,
        "observationWindowMeaning": window_meaning,
        "ascertainmentMethod": ascertainment,
        "activeTBClassification": tb_classification,
        "source": str(_required(row, "source")).strip(),
        "reviewStatus": str(_required(row, "reviewStatus")).strip(),
        "notes": str(row.get("notes") or ""),
        "uncertainty": {str(key): value for key, value in uncertainty.items()},
        "observedRatePer100000Population": observed_count / population_denominator * 100000.0,
        "observedRatePer100000PersonYears": (
            None if person_years is None else observed_count / person_years * 100000.0
        ),
    }
    if not out["source"]:
        raise ValueError("source must not be empty.")
    if not out["reviewStatus"]:
        raise ValueError("reviewStatus must not be empty.")
    return out


def validate_active_tb_observations(rows: Iterable[Mapping[str, Any]]) -> tuple[dict[str, Any], ...]:
    return tuple(validate_active_tb_observation(row) for row in rows)


def _solve_hazard_for_target(
    durations: np.ndarray,
    weights: np.ndarray,
    target: float,
    max_value: float,
    label: str,
    *,
    tolerance: float = ROOT_TOLERANCE,
) -> float:
    if target > max_value + PROBABILITY_TOLERANCE:
        raise CalibrationError(
            f"Requested {label} prevalence exceeds its maximum possible value."
        )
    if target <= tolerance:
        return 0.0

    def value(hazard: float) -> float:
        return float(np.sum(weights * (1.0 - np.exp(-hazard * durations))))

    return _bisect_nonnegative_monotone(value, target, label, tolerance=tolerance)


def _solve_remote_hazard_for_target(
    remote_durations: np.ndarray,
    weights: np.ndarray,
    no_recent: np.ndarray,
    target: float,
    max_value: float,
    *,
    tolerance: float = ROOT_TOLERANCE,
) -> float:
    if target > max_value + PROBABILITY_TOLERANCE:
        raise CalibrationError(
            "Requested remote-only prevalence exceeds its maximum possible value."
        )
    if target <= tolerance:
        return 0.0

    def value(hazard: float) -> float:
        return float(
            np.sum(weights * no_recent * (1.0 - np.exp(-hazard * remote_durations)))
        )

    return _bisect_nonnegative_monotone(value, target, "remote-only", tolerance=tolerance)


def _bisect_nonnegative_monotone(
    value,
    target: float,
    label: str,
    *,
    tolerance: float,
) -> float:
    lo = 0.0
    hi = 1.0
    while value(hi) < target and hi < MAX_HAZARD_PER_YEAR:
        hi *= 2.0
    if value(hi) < target:
        raise CalibrationError(
            f"Could not bracket {label} hazard below {MAX_HAZARD_PER_YEAR:g} per year."
        )
    for _ in range(MAX_ROOT_ITERATIONS):
        mid = (lo + hi) / 2.0
        mid_value = value(mid)
        if abs(mid_value - target) <= tolerance or (hi - lo) / 2.0 <= tolerance:
            return mid
        if mid_value < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def _parse_observation_endpoint(row: Mapping[str, Any], prefix: str) -> dict[str, Any]:
    date_key = f"{prefix}Date"
    year_key = f"{prefix}Year"
    has_date = row.get(date_key) not in (None, "")
    has_year = row.get(year_key) not in (None, "")
    if not has_date and not has_year:
        raise ValueError(f"{date_key} or {year_key} is required.")
    if has_date:
        parsed = _parse_iso_date(row[date_key], date_key)
        if has_year and int(row[year_key]) != parsed.year:
            raise ValueError(f"{date_key} and {year_key} disagree.")
        return {"date": parsed.isoformat(), "year": parsed.year, "month": parsed.month, "day": parsed.day}
    year = int(row[year_key])
    if year < 0:
        raise ValueError(f"{year_key} must be non-negative.")
    return {
        "date": None,
        "year": year,
        "month": 1 if prefix == "start" else 12,
        "day": 1 if prefix == "start" else 31,
    }


def _parse_iso_date(value: Any, field: str) -> date:
    try:
        return datetime.strptime(str(value), "%Y-%m-%d").date()
    except ValueError as exc:
        raise ValueError(f"{field} must use YYYY-MM-DD format.") from exc


def _endpoint_sort_key(endpoint: Mapping[str, Any]) -> tuple[int, int, int]:
    return (int(endpoint["year"]), int(endpoint["month"]), int(endpoint["day"]))


def _required(row: Mapping[str, Any], field: str) -> Any:
    value = row.get(field)
    if value in (None, ""):
        raise ValueError(f"{field} is required.")
    return value


def _one_of(row: Mapping[str, Any], field: str, allowed: set[str]) -> str:
    value = str(_required(row, field)).strip()
    if value not in allowed:
        raise ValueError(f"{field} must be one of: {', '.join(sorted(allowed))}.")
    return value


def _finite_nonnegative_tuple(values: Iterable[float], label: str) -> tuple[float, ...]:
    try:
        out = tuple(float(value) for value in values)
    except TypeError as exc:
        raise ValueError(f"{label} must be an iterable of numbers.") from exc
    if any(not math.isfinite(value) for value in out):
        raise ValueError(f"{label} must contain only finite values.")
    if any(value < 0 for value in out):
        raise ValueError(f"{label} must contain only non-negative values.")
    return out


def _finite_nonnegative_float(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"{label} must be finite and non-negative.")
    return number


def _positive_float(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{label} must be finite and positive.")
    return number


def _probability(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0 or number > 1.0:
        raise ValueError(f"{label} must be a finite probability in [0,1].")
    return number


def _validate_probability_vector(values: np.ndarray, label: str, tolerance: float) -> None:
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{label} probabilities must be finite.")
    if np.any(values < -tolerance) or np.any(values > 1.0 + tolerance):
        raise ValueError(f"{label} probabilities must lie between zero and one.")


def _tuple(values: np.ndarray) -> tuple[float, ...]:
    return tuple(float(value) for value in values)
