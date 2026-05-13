from __future__ import annotations

import math
from datetime import datetime, timezone
from typing import Any


ALIGNED_METRICS = [
    ("projection_horizon", "horizon", "horizon", "Projection/follow-up horizon."),
    ("population", "population", "population", "Population size if available."),
    (
        "cumulative_baseline_active_tb_cases",
        "cumulative_baseline_active_tb_cases",
        "cumulative baseline active TB cases",
        "Only comparable when both models expose a baseline/no-intervention cumulative count.",
    ),
    (
        "cumulative_intervention_active_tb_cases",
        "cumulative_intervention_active_tb_cases",
        "cumulative intervention active TB cases",
        "Only comparable when both models expose an intervention cumulative count.",
    ),
    ("cumulative_cases_averted", "cumulative_cases_averted", "cumulative cases averted", "Closest aligned benefit metric."),
    (
        "relative_reduction_cumulative_active_tb_cases",
        "relative_reduction",
        "relative reduction",
        "Relative reduction in cumulative active TB cases where available.",
    ),
]

STRUCTURALLY_NON_COMPARABLE_METRICS = [
    "NNS/NNT: APY ABM screening efficiency outputs are not directly comparable to dynamic incidence projections.",
    "False positives treated: APY ABM pathway outputs have no direct dynamic model equivalent.",
    "Treatment starts/completions: APY ABM cascade outputs have no direct dynamic model equivalent.",
    "Risk-factor attributable outputs: APY ABM attributable outputs are not interchangeable with dynamic risk-factor inputs.",
    "ABM economics outputs: APY ABM cost outputs are optional and are not directly comparable to dynamic epidemiologic projections.",
]


def number_or_none(value: Any) -> float | int | None:
    try:
        if value is None or value == "":
            return None
        result = float(value)
        if not math.isfinite(result):
            return None
        if result.is_integer():
            return int(result)
        return result
    except (TypeError, ValueError):
        return None


def metric_rows(bundle: dict) -> list[dict[str, Any]]:
    headline = bundle.get("headline") or {}
    rows = headline.get("keyMetricsRows") or headline.get("summaryRows") or []
    return rows if isinstance(rows, list) else []


def row_lookup(bundle: dict) -> dict[str, Any]:
    lookup: dict[str, Any] = {}
    for row in metric_rows(bundle):
        if not isinstance(row, dict):
            continue
        metric = row.get("Metric") or row.get("metric")
        if metric in (None, ""):
            continue
        value = row.get("Value", row.get("Median", row.get("value")))
        lookup[str(metric)] = value
    return lookup


def first_metric(lookup: dict[str, Any], names: list[str]) -> Any:
    for name in names:
        if name in lookup:
            return lookup[name]
    return None


def extract_dynamic_key_metrics(bundle: dict) -> dict[str, Any]:
    lookup = row_lookup(bundle or {})
    metadata = (bundle or {}).get("metadata") or {}
    return {
        "projection_horizon": first_metric(lookup, ["projection_horizon"]),
        "population": metadata.get("population"),
        "cumulative_baseline_active_tb_cases": first_metric(lookup, ["cumulative_baseline_active_tb_cases"]),
        "cumulative_intervention_active_tb_cases": first_metric(lookup, ["cumulative_intervention_active_tb_cases"]),
        "cumulative_cases_averted": first_metric(lookup, ["cumulative_cases_averted"]),
        "relative_reduction_cumulative_active_tb_cases": first_metric(
            lookup,
            ["relative_reduction_cumulative_active_tb_cases", "relative_reduction"],
        ),
    }


def extract_abm_key_metrics(bundle: dict) -> dict[str, Any]:
    lookup = row_lookup(bundle or {})
    technical = (bundle or {}).get("technical") or {}
    interface_config = technical.get("interfaceConfig") or {}
    return {
        "horizon": interface_config.get("followHorizon"),
        "population": interface_config.get("N"),
        "cumulative_baseline_active_tb_cases": first_metric(
            lookup,
            ["cumulative_baseline_active_tb_cases", "baselineActiveTB", "nBaselineActiveTB"],
        ),
        "cumulative_intervention_active_tb_cases": first_metric(
            lookup,
            ["cumulative_intervention_active_tb_cases", "interventionActiveTB", "nInterventionActiveTB"],
        ),
        "cumulative_cases_averted": first_metric(
            lookup,
            ["cumulative_cases_averted", "nPreventedActiveTB", "preventedActiveTB"],
        ),
        "relative_reduction": first_metric(
            lookup,
            ["relative_reduction", "relativeReduction", "relative_reduction_cumulative_active_tb_cases"],
        ),
    }


def comparison_row(metric: str, dynamic_value: Any, abm_value: Any, notes: str) -> dict[str, Any]:
    dynamic_number = number_or_none(dynamic_value)
    abm_number = number_or_none(abm_value)
    absolute_difference = None
    relative_difference = None
    comparable = dynamic_number is not None and abm_number is not None
    if comparable:
        absolute_difference = dynamic_number - abm_number
        if abm_number != 0:
            relative_difference = absolute_difference / abm_number
    return {
        "metric": metric,
        "dynamic_value": dynamic_value,
        "abm_value": abm_value,
        "absolute_difference": absolute_difference,
        "relative_difference": relative_difference,
        "comparable": comparable,
        "notes": notes if comparable else f"{notes} Missing or non-numeric value in one model.",
    }


def compare_dynamic_abm_v9(dynamic_bundle: dict, abm_bundle: dict) -> dict[str, Any]:
    dynamic_metrics = extract_dynamic_key_metrics(dynamic_bundle or {})
    abm_metrics = extract_abm_key_metrics(abm_bundle or {})
    rows = []
    missing_dynamic_metrics = []
    missing_abm_metrics = []
    for dynamic_key, abm_key, label, notes in ALIGNED_METRICS:
        dynamic_value = dynamic_metrics.get(dynamic_key)
        abm_value = abm_metrics.get(abm_key)
        rows.append(comparison_row(label, dynamic_value, abm_value, notes))
        if number_or_none(dynamic_value) is None:
            missing_dynamic_metrics.append(f"{label} ({dynamic_key})")
        if number_or_none(abm_value) is None:
            missing_abm_metrics.append(f"{label} ({abm_key})")

    warnings = []
    if missing_dynamic_metrics:
        warnings.append("Missing or non-numeric dynamic metrics: " + ", ".join(missing_dynamic_metrics) + ".")
    if missing_abm_metrics:
        warnings.append("Missing or non-numeric APY ABM metrics: " + ", ".join(missing_abm_metrics) + ".")
    warnings.extend(STRUCTURALLY_NON_COMPARABLE_METRICS)
    return {
        "comparisonRows": rows,
        "warnings": warnings,
        "missing_dynamic_metrics": missing_dynamic_metrics,
        "missing_abm_metrics": missing_abm_metrics,
        "structurally_non_comparable_metrics": list(STRUCTURALLY_NON_COMPARABLE_METRICS),
        "metadata": {
            "contractVersion": "dynamic_abm_compare_v9",
            "createdAt": datetime.now(timezone.utc).isoformat(),
            "dynamicModelVersion": (dynamic_bundle or {}).get("modelVersion"),
            "abmModelVersion": ((abm_bundle or {}).get("metadata") or {}).get("modelVersion"),
        },
    }
