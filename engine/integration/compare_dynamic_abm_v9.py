from __future__ import annotations

import math
from statistics import median
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


def rows_lookup(rows: list[dict[str, Any]]) -> dict[str, Any]:
    lookup: dict[str, Any] = {}
    for row in rows if isinstance(rows, list) else []:
        if not isinstance(row, dict):
            continue
        metric = row.get("Metric") or row.get("metric")
        if metric in (None, ""):
            continue
        lookup[str(metric)] = row.get("Value", row.get("Median", row.get("value")))
    return lookup


def first_metric(lookup: dict[str, Any], names: list[str]) -> Any:
    for name in names:
        if name in lookup:
            return lookup[name]
    return None


def nested_dict(bundle: dict, *keys: str) -> dict[str, Any]:
    value: Any = bundle
    for key in keys:
        if not isinstance(value, dict):
            return {}
        value = value.get(key)
    return value if isinstance(value, dict) else {}


def numeric_median(values: list[Any]) -> float | int | None:
    nums = [number_or_none(value) for value in values]
    nums = [value for value in nums if value is not None]
    if not nums:
        return None
    result = float(median(nums))
    if result.is_integer():
        return int(result)
    return result


def median_column(rows: list[dict[str, Any]], column: str) -> float | int | None:
    return numeric_median([row.get(column) for row in rows if isinstance(row, dict)])


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
    technical = nested_dict(bundle or {}, "technical")
    interface_config = technical.get("interfaceConfig") or {}
    dynamic_comparison = technical.get("dynamicComparison") if isinstance(technical.get("dynamicComparison"), dict) else {}
    dynamic_comparison_rows = rows_lookup(dynamic_comparison.get("metricRows") or [])
    do_nothing = nested_dict(bundle or {}, "doNothing")
    derived_rows = do_nothing.get("derivedRows") or []
    summary_rows = rows_lookup(do_nothing.get("summaryRows") or [])

    baseline = first_metric(
        dynamic_comparison_rows,
        ["cumulative_baseline_active_tb_cases"],
    )
    baseline = dynamic_comparison.get("cumulative_baseline_active_tb_cases", baseline)
    baseline = baseline if number_or_none(baseline) is not None else first_metric(
        lookup,
        ["cumulative_baseline_active_tb_cases", "baselineActiveTB", "nBaselineActiveTB"],
    )
    baseline = baseline if number_or_none(baseline) is not None else median_column(derived_rows, "nActiveBy20y_DoNothing")

    intervention = first_metric(
        dynamic_comparison_rows,
        ["cumulative_intervention_active_tb_cases"],
    )
    intervention = dynamic_comparison.get("cumulative_intervention_active_tb_cases", intervention)
    intervention = intervention if number_or_none(intervention) is not None else first_metric(
        lookup,
        ["cumulative_intervention_active_tb_cases", "interventionActiveTB", "nInterventionActiveTB"],
    )
    intervention = intervention if number_or_none(intervention) is not None else median_column(derived_rows, "nActiveBy20y_AfterStrategy")

    averted = first_metric(dynamic_comparison_rows, ["cumulative_cases_averted"])
    averted = dynamic_comparison.get("cumulative_cases_averted", averted)
    averted = averted if number_or_none(averted) is not None else first_metric(
        lookup,
        ["cumulative_cases_averted", "preventedActiveTB"],
    )
    averted = averted if number_or_none(averted) is not None else median_column(derived_rows, "nActiveBy20y_Prevented")
    averted = averted if number_or_none(averted) is not None else first_metric(
        lookup,
        ["nPreventedActiveTB"],
    )

    rel_reduction = first_metric(dynamic_comparison_rows, ["relative_reduction_cumulative_active_tb_cases"])
    rel_reduction = dynamic_comparison.get("relative_reduction_cumulative_active_tb_cases", rel_reduction)
    rel_reduction = rel_reduction if number_or_none(rel_reduction) is not None else first_metric(
        lookup,
        ["relative_reduction", "relativeReduction", "relative_reduction_cumulative_active_tb_cases"],
    )
    rel_reduction = rel_reduction if number_or_none(rel_reduction) is not None else first_metric(
        summary_rows,
        ["Relative reduction in 20-year active TB burden"],
    )
    rel_reduction = rel_reduction if number_or_none(rel_reduction) is not None else median_column(derived_rows, "relReduction20y")

    return {
        "horizon": interface_config.get("followHorizon"),
        "population": interface_config.get("N"),
        "cumulative_baseline_active_tb_cases": baseline,
        "cumulative_intervention_active_tb_cases": intervention,
        "cumulative_cases_averted": averted,
        "relative_reduction": rel_reduction,
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
    dynamic_horizon = number_or_none(dynamic_metrics.get("projection_horizon"))
    abm_horizon = number_or_none(abm_metrics.get("horizon"))
    dynamic_population = number_or_none(dynamic_metrics.get("population"))
    abm_population = number_or_none(abm_metrics.get("population"))
    if dynamic_horizon is not None and abm_horizon is not None and dynamic_horizon != abm_horizon:
        warnings.append(f"Dynamic horizon ({dynamic_horizon}) differs from APY followHorizon ({abm_horizon}).")
    if dynamic_population is not None and abm_population is not None and dynamic_population != abm_population:
        warnings.append(f"Dynamic population ({dynamic_population}) differs from APY population N ({abm_population}).")
        comparable_count_metrics = [
            row["metric"]
            for row in rows
            if row["comparable"] and row["metric"] in {
                "cumulative baseline active TB cases",
                "cumulative intervention active TB cases",
                "cumulative cases averted",
            }
        ]
        if comparable_count_metrics:
            warnings.append(
                "Count metrics are compared despite different population sizes: "
                + ", ".join(comparable_count_metrics)
                + "."
            )
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
