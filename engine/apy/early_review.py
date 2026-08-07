from __future__ import annotations

from copy import deepcopy
import math
from typing import Any

import numpy as np
import pandas as pd

from engine.apy.decision_analysis import run_scenario_comparison
from engine.apy.event_ledger_economics import PRIMARY_DISCOUNT_RATE


EARLY_REVIEW_CONTRACT_VERSION = "apy_early_screening_review_v1"


def run_early_screening_review(
    base_config: dict[str, Any],
    economics_config: dict[str, Any] | None,
    review_input: dict[str, Any],
) -> dict[str, Any]:
    validation = _validate_review_input(base_config, review_input)
    if validation["errors"]:
        return _empty_result(review_input, validation)
    prior = _prior_weights(review_input)
    if prior["errors"]:
        validation["errors"].extend(prior["errors"])
        return _empty_result(review_input, validation)

    grid_rows = []
    projection_rows = []
    log_likelihoods = []
    projection_cache: dict[tuple[float, float], dict[str, Any]] = {}
    screened_to_date = float(review_input["screenedToDate"])
    positives = float(review_input["testPositiveToDate"])
    planned_total = float(review_input["plannedTotalScreened"])
    eligible = float(review_input.get("eligiblePopulation", base_config.get("N")))
    completed_coverage = 0.0 if eligible <= 0 else screened_to_date / eligible
    planned_coverage = 0.0 if eligible <= 0 else planned_total / eligible

    for prevalence, prior_weight in zip(prior["prevalenceGrid"], prior["weights"]):
        completed = _cached_projection(
            projection_cache,
            base_config,
            economics_config,
            prevalence,
            completed_coverage,
            "stop_after_review",
        )
        planned = _cached_projection(
            projection_cache,
            base_config,
            economics_config,
            prevalence,
            planned_coverage,
            "continue_to_planned_screening",
        )
        p_positive = _positive_probability(completed["summary"])
        log_likelihood = _binomial_log_likelihood(screened_to_date, positives, p_positive)
        log_likelihoods.append(log_likelihood)
        grid_rows.append(
            {
                "prevalence": prevalence,
                "priorWeight": prior_weight,
                "predictedPositiveProbability": p_positive,
                "logLikelihood": log_likelihood,
                "observedPositiveYield": None if screened_to_date == 0 else positives / screened_to_date,
            }
        )
        projection_rows.append(
            _projection_delta_row(
                prevalence,
                completed,
                planned,
                completed_coverage,
                planned_coverage,
            )
        )

    posterior_weights = _posterior_weights(prior["weights"], log_likelihoods)
    for row, weight in zip(grid_rows, posterior_weights):
        row["posteriorWeight"] = weight
    for row, weight in zip(projection_rows, posterior_weights):
        row["posteriorWeight"] = weight

    prior_summary = _weighted_summary(prior["prevalenceGrid"], prior["weights"])
    posterior_summary = _weighted_summary(prior["prevalenceGrid"], posterior_weights)
    posterior_projection = _posterior_projection_summary(projection_rows, posterior_weights)
    informativeness_warnings = _test_informativeness_warnings(base_config)
    validation["warnings"].extend(informativeness_warnings)
    direction = (
        "Selected test characteristics are non-informative or direction-reversed; standard yield-to-prevalence interpretation is suppressed."
        if informativeness_warnings
        else _yield_direction(prior_summary.get("mean"), posterior_summary.get("mean"))
    )
    timing_approximation = True
    timing_reasons = _economic_timing_approximation_reasons()
    return {
        "contractVersion": EARLY_REVIEW_CONTRACT_VERSION,
        "reviewId": review_input.get("reviewId", "early_review"),
        "scenarioId": review_input.get("scenarioId", base_config.get("scenarioLabel")),
        "inputs": deepcopy(review_input),
        "prior": {
            "type": review_input.get("prior", {}).get("type"),
            "summary": prior_summary,
            "provenance": review_input.get("prior", {}).get("source", ""),
            "discretisationMethod": prior.get("discretisationMethod"),
        },
        "posterior": {
            "summary": posterior_summary,
            "probabilityBelowThresholds": _probability_below_thresholds(
                prior["prevalenceGrid"],
                posterior_weights,
                review_input.get("prevalenceThresholds") or [],
            ),
            "interpretation": direction,
        },
        "likelihoodMethod": "aggregate_binomial_approximation",
        "likelihoodNotes": (
            "positives ~ Binomial(screened_to_date, model-predicted positive probability "
            "for the screened tranche); this is an aggregate approximation for targeted "
            "heterogeneous screening."
        ),
        "priorPosteriorTable": grid_rows,
        "projection": projection_rows,
        "posteriorProjectionSummary": posterior_projection,
        "stopStrategy": "stop_after_review",
        "continueStrategy": "continue_to_planned_screening",
        "timingApproximation": timing_approximation,
        "timingApproximationDescription": (
            "Stop-versus-continue is represented by completed versus planned total "
            "coverage using the same deterministic targeting order. Completed screening "
            "is not erased; remaining programme consequences are planned minus completed. "
            "Within-window review-time scheduling is not yet exact."
        ),
        "economicTimingStatus": "approximate_not_decision_grade",
        "economicDecisionComplete": False,
        "economicTimingApproximation": True,
        "economicTimingApproximationReasons": timing_reasons,
        "conclusionPermitted": False,
        "posteriorProbabilityPositiveNMB": None,
        "economicConclusion": None,
        "approximateEconomicComponentsLabel": "Approximate economic components - timing not fully represented",
        "validation": validation,
        "calibrationPolicy": "infection_intercept_only",
        "referenceCalibrationHash": _first_projection_field(projection_rows, "referenceCalibrationHash"),
        "activeTBSurveillanceJointUpdate": False,
        "activeTBSurveillanceJointUpdateNotes": (
            "Screening-yield updating changes the LTBI-prevalence distribution only; "
            "active-TB surveillance data are not jointly updated in this workflow."
        ),
    }


def _run_coverage_projection(
    base_config: dict[str, Any],
    economics_config: dict[str, Any] | None,
    prevalence: float,
    coverage: float,
    scenario_id: str,
) -> dict[str, Any]:
    comparison = run_scenario_comparison(
        base_config,
        economics_config,
        [
            {
                "scenarioId": scenario_id,
                "changes": {
                    "ltbiPrevalence": float(max(prevalence, 0.0)),
                    "screenCoverage": float(max(0.0, min(1.0, coverage))),
                },
                "calibrationPolicy": "infection_intercept_only",
            }
        ],
        model_type="expected_value",
    )
    return {
        "comparison": comparison,
        "summary": comparison["scenarioSummaries"][0],
    }


def _cached_projection(
    cache: dict[tuple[float, float], dict[str, Any]],
    base_config: dict[str, Any],
    economics_config: dict[str, Any] | None,
    prevalence: float,
    coverage: float,
    scenario_id: str,
) -> dict[str, Any]:
    key = (round(float(prevalence), 12), round(float(coverage), 12))
    if key not in cache:
        cache[key] = _run_coverage_projection(
            base_config,
            economics_config,
            prevalence,
            coverage,
            scenario_id,
        )
    return cache[key]


def _projection_delta_row(
    prevalence: float,
    completed: dict[str, Any],
    planned: dict[str, Any],
    completed_coverage: float,
    planned_coverage: float,
) -> dict[str, Any]:
    completed_summary = completed["summary"]
    planned_summary = planned["summary"]
    row = {
        "prevalence": prevalence,
        "completedCoverage": completed_coverage,
        "plannedCoverage": planned_coverage,
        "additionalCoverage": planned_coverage - completed_coverage,
    }
    metric_map = {
        "additionalPeopleScreened": "screened",
        "additionalTruePositiveRecent": "true_positive_recent",
        "additionalTruePositiveRemote": "true_positive_remote",
        "additionalFalsePositives": "false_positive",
        "additionalTPTStarts": "tpt_started_total",
        "additionalTPTCompletions": "tpt_completed_total",
        "additionalInfectionsEffectivelyTreated": "infection_effectively_treated_total",
        "additionalActiveTBCasesPrevented": "active_tb_cases_prevented",
        "plannedComparatorActiveTB": "comparator_active_tb",
        "plannedInterventionActiveTB": "intervention_active_tb",
        "plannedActiveTBCasesPrevented": "active_tb_cases_prevented",
        "additionalIncrementalCostOfContinuing": "incrementalCost",
        "additionalDALYsAverted": "dalysAverted",
    }
    for out_key, source_key in metric_map.items():
        if out_key.startswith("planned"):
            row[out_key] = planned_summary.get(source_key)
        else:
            row[out_key] = _subtract(planned_summary.get(source_key), completed_summary.get(source_key))
    row["economicsComplete"] = planned_summary.get("incrementalCost") is not None and completed_summary.get("incrementalCost") is not None
    row["economicTimingStatus"] = "approximate_not_decision_grade"
    row["economicDecisionComplete"] = False
    row["economicTimingApproximation"] = True
    row["economicTimingApproximationReasons"] = _economic_timing_approximation_reasons()
    row["conclusionPermitted"] = False
    row["additionalNMB"] = None
    row["posteriorProbabilityPositiveNMB"] = None
    row["economicConclusion"] = None
    row["deprecatedAdditionalProgrammeCost"] = row.get("additionalIncrementalCostOfContinuing")
    row["deprecatedAdditionalProgrammeCostNotes"] = (
        "Deprecated compatibility alias. Authoritative field is additionalIncrementalCostOfContinuing; "
        "the value is planned incremental cost minus completed-coverage incremental cost and may include "
        "non-programme health-system consequences."
    )
    row["additionalProgrammeSetupCost"] = 0.0 if completed_coverage > 0 else None
    row["additionalProgrammeSetupCostNotes"] = (
        "Setup expenditure is treated as already incurred after screening has begun and is not charged again as a future continuation cost."
        if completed_coverage > 0
        else "Future setup timing is unavailable under the current approximate review-time allocation."
    )
    row["calibrationPolicy"] = planned_summary.get("calibrationPolicy")
    row["referenceCalibrationHash"] = planned_summary.get("referenceCalibrationHash")
    row["timingApproximation"] = True
    row["fixedProgrammeCostTimingApproximate"] = True
    row["variableProgrammeCostTimingApproximate"] = True
    return row


def _validate_review_input(base_config: dict[str, Any], review: dict[str, Any]) -> dict[str, Any]:
    errors = []
    warnings = []
    required = ["screenedToDate", "testPositiveToDate", "plannedTotalScreened"]
    for field in required:
        if field not in review:
            errors.append({"field": field, "message": f"{field} is required."})
    if errors:
        return {"isValid": False, "errors": errors, "warnings": warnings}
    screened = _count_value(review["screenedToDate"])
    positives = _count_value(review["testPositiveToDate"])
    planned = _count_value(review["plannedTotalScreened"])
    eligible = _count_value(review.get("eligiblePopulation", base_config.get("N", planned)))
    if any(value is None for value in [screened, positives, planned, eligible]):
        errors.append(
            {
                "field": "screeningYield",
                "message": "Screened, positive, planned and eligible values must be integer-valued counts.",
            }
        )
        return {"isValid": False, "errors": errors, "warnings": warnings}
    if not 0 <= positives <= screened <= planned <= eligible:
        errors.append(
            {
                "field": "screeningYield",
                "message": "Require 0 <= positives <= screened_to_date <= planned_total <= eligible_population.",
            }
        )
    if "prior" not in review:
        errors.append({"field": "prior", "message": "An explicit prior is required; no hidden prior strength is supplied."})
    return {"isValid": not errors, "errors": errors, "warnings": warnings}


def _prior_weights(review: dict[str, Any]) -> dict[str, Any]:
    prior = review.get("prior") or {}
    grid = np.asarray(review.get("prevalenceGrid") or prior.get("prevalenceGrid") or [], dtype=float)
    errors = []
    if grid.size == 0:
        errors.append({"field": "prevalenceGrid", "message": "An explicit prevalence grid is required."})
        return {"errors": errors}
    if np.any(grid < 0) or np.any(grid > 1):
        errors.append({"field": "prevalenceGrid", "message": "Prevalence grid values must be in [0,1]."})
    if grid.size > 1 and np.any(np.diff(grid) <= 0):
        errors.append({"field": "prevalenceGrid", "message": "Prevalence grid must be strictly increasing."})
    prior_type = prior.get("type")
    if prior_type == "beta":
        if grid.size > 2 and not np.allclose(np.diff(grid), np.diff(grid)[0], rtol=1e-8, atol=1e-12):
            errors.append(
                {
                    "field": "prevalenceGrid",
                    "message": "Beta-prior prevalence grid must be evenly spaced for point-mass discretisation.",
                }
            )
        alpha = prior.get("alpha")
        beta = prior.get("beta")
        if alpha in (None, "") or beta in (None, ""):
            mean = prior.get("mean")
            ess = prior.get("effectiveSampleSize")
            if mean in (None, "") or ess in (None, ""):
                errors.append({"field": "prior", "message": "Beta prior requires alpha/beta or mean/effective sample size."})
            else:
                mean = float(mean)
                ess = float(ess)
                if not 0 < mean < 1:
                    errors.append({"field": "prior.mean", "message": "Beta prior mean must be strictly between 0 and 1."})
                if ess <= 0:
                    errors.append({"field": "prior.effectiveSampleSize", "message": "Beta prior effective sample size must be > 0."})
                alpha = mean * ess
                beta = (1.0 - mean) * ess
        if not errors and (float(alpha) <= 0 or float(beta) <= 0):
            errors.append({"field": "prior", "message": "Beta prior alpha and beta must be > 0."})
        if not errors:
            weights = np.asarray([_beta_pdf(x, float(alpha), float(beta)) for x in grid], dtype=float)
    elif prior_type == "discrete_grid":
        weights = np.asarray(prior.get("weights") or [], dtype=float)
        if weights.size != grid.size:
            errors.append({"field": "prior.weights", "message": "Discrete-grid prior weights must match the prevalence grid."})
    else:
        errors.append({"field": "prior.type", "message": "Supported priors are beta and discrete_grid."})
    if errors:
        return {"errors": errors}
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        errors.append({"field": "prior.weights", "message": "Prior weights must be finite and sum to a positive value."})
        return {"errors": errors}
    weights = weights / weights.sum()
    method = (
        "explicit_discrete_grid_weights"
        if prior_type == "discrete_grid"
        else "evenly_spaced_point_mass_approximation"
    )
    return {
        "errors": [],
        "prevalenceGrid": grid.tolist(),
        "weights": weights.tolist(),
        "discretisationMethod": method,
    }


def _posterior_weights(prior_weights: list[float], log_likelihoods: list[float]) -> list[float]:
    log_prior = np.log(np.maximum(np.asarray(prior_weights, dtype=float), np.finfo(float).tiny))
    log_like = np.asarray(log_likelihoods, dtype=float)
    log_post = log_prior + log_like
    log_post -= np.nanmax(log_post)
    weights = np.exp(log_post)
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        weights = np.asarray(prior_weights, dtype=float)
    return (weights / weights.sum()).tolist()


def _binomial_log_likelihood(n: float, k: float, p: float | None) -> float:
    if n == 0:
        return 0.0
    if p is None:
        return -math.inf
    p = min(max(float(p), 1e-15), 1 - 1e-15)
    return k * math.log(p) + (n - k) * math.log1p(-p)


def _positive_probability(summary: dict[str, Any]) -> float | None:
    screened = summary.get("screened")
    positives = summary.get("test_positive_total")
    if screened in (None, 0) or positives is None:
        return None
    return max(0.0, min(1.0, float(positives) / float(screened)))


def _weighted_summary(values: list[float], weights: list[float]) -> dict[str, Any]:
    values_arr = np.asarray(values, dtype=float)
    weights_arr = np.asarray(weights, dtype=float)
    order = np.argsort(values_arr)
    values_arr = values_arr[order]
    weights_arr = weights_arr[order] / weights_arr.sum()
    cdf = np.cumsum(weights_arr)
    return {
        "mean": float(np.sum(values_arr * weights_arr)),
        "median": float(values_arr[np.searchsorted(cdf, 0.5, side="left")]),
        "p2_5": float(values_arr[np.searchsorted(cdf, 0.025, side="left")]),
        "p97_5": float(values_arr[np.searchsorted(cdf, 0.975, side="left")]),
    }


def _posterior_projection_summary(rows: list[dict[str, Any]], weights: list[float]) -> list[dict[str, Any]]:
    metrics = [
        "additionalPeopleScreened",
        "additionalTruePositiveRecent",
        "additionalTruePositiveRemote",
        "additionalFalsePositives",
        "additionalTPTStarts",
        "additionalTPTCompletions",
        "additionalInfectionsEffectivelyTreated",
        "additionalActiveTBCasesPrevented",
        "plannedComparatorActiveTB",
        "plannedInterventionActiveTB",
        "plannedActiveTBCasesPrevented",
        "additionalIncrementalCostOfContinuing",
        "additionalDALYsAverted",
    ]
    out = []
    for metric in metrics:
        values = [row.get(metric) for row in rows]
        if any(value is None or pd.isna(value) for value in values):
            out.append({"metric": metric, "mean": None, "p2_5": None, "p97_5": None, "complete": False})
            continue
        summary = _weighted_summary([float(value) for value in values], weights)
        out.append({"metric": metric, **summary, "complete": True})
    return out


def _probability_below_thresholds(values: list[float], weights: list[float], thresholds: list[float]) -> list[dict[str, Any]]:
    out = []
    for threshold in thresholds:
        out.append(
            {
                "threshold": float(threshold),
                "probability": float(
                    sum(weight for value, weight in zip(values, weights) if float(value) < float(threshold))
                ),
            }
        )
    return out


def _test_informativeness_warnings(config: dict[str, Any]) -> list[dict[str, str]]:
    if str(config.get("testType", "IGRA")).upper() == "TST":
        sens = float(config.get("tstSensitivity", 0.0))
        spec = min(float(config.get("tstSpecificityBCG", 0.0)), float(config.get("tstSpecificityNoBCG", 0.0)))
    else:
        sens = float(config.get("testSensitivity", 0.0))
        spec = float(config.get("testSpecificity", 0.0))
    if sens + spec <= 1:
        return [{"field": "testCharacteristics", "message": "Sensitivity + specificity <= 1; directional prevalence updating is not informative."}]
    return []


def _yield_direction(prior_mean: float | None, posterior_mean: float | None) -> str:
    if prior_mean is None or posterior_mean is None:
        return ""
    if posterior_mean < prior_mean:
        return "Observed positive yield shifts the prevalence distribution downward under the selected prior and test assumptions."
    if posterior_mean > prior_mean:
        return "Observed positive yield shifts the prevalence distribution upward under the selected prior and test assumptions."
    return "Observed positive yield leaves the posterior mean unchanged at the displayed precision."


def _economic_timing_approximation_reasons() -> list[str]:
    return [
        "Exact review-time allocation of programme setup expenditure is not represented.",
        "Programme-running costs already incurred versus remaining after review are not separately scheduled.",
        "Already-started versus future TPT courses are approximated by completed versus planned total coverage.",
        "Future active-TB and discount timing are inherited from coverage scenarios rather than an exact review-time schedule.",
    ]


def _count_value(value: Any) -> int | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number) or not number.is_integer():
        return None
    return int(number)


def _first_projection_field(rows: list[dict[str, Any]], field: str) -> Any:
    for row in rows:
        value = row.get(field)
        if value not in (None, ""):
            return value
    return None


def _beta_pdf(x: float, alpha: float, beta: float) -> float:
    if x <= 0 or x >= 1 or alpha <= 0 or beta <= 0:
        if x == 0 and alpha < 1:
            return math.inf
        if x == 1 and beta < 1:
            return math.inf
        if alpha <= 0 or beta <= 0:
            return 0.0
    log_pdf = (
        (alpha - 1.0) * math.log(max(x, 1e-15))
        + (beta - 1.0) * math.log(max(1.0 - x, 1e-15))
        + math.lgamma(alpha + beta)
        - math.lgamma(alpha)
        - math.lgamma(beta)
    )
    return math.exp(min(log_pdf, 700))


def _subtract(a: Any, b: Any) -> float | None:
    if a is None or b is None:
        return None
    try:
        if pd.isna(a) or pd.isna(b):
            return None
        return float(a) - float(b)
    except (TypeError, ValueError):
        return None


def _empty_result(review_input: dict[str, Any], validation: dict[str, Any]) -> dict[str, Any]:
    return {
        "contractVersion": EARLY_REVIEW_CONTRACT_VERSION,
        "reviewId": review_input.get("reviewId", "early_review"),
        "inputs": deepcopy(review_input),
        "priorPosteriorTable": [],
        "projection": [],
        "posteriorProjectionSummary": [],
        "validation": {**validation, "isValid": False},
    }
