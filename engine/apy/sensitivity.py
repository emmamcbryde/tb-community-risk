from __future__ import annotations

from copy import deepcopy
import csv
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from engine.apy.decision_analysis import (
    DECISION_ANALYSIS_CONTRACT_VERSION,
    build_scenario_config,
    run_scenario_comparison,
)


SENSITIVITY_CONTRACT_VERSION = "apy_one_way_sensitivity_v1"
THRESHOLD_CONTRACT_VERSION = "apy_threshold_analysis_v1"
DEFAULT_SENSITIVITY_SPEC_PATH = Path(__file__).resolve().parents[2] / "data" / "apy_sensitivity_spec.csv"


def load_sensitivity_specs(path: str | Path | None = None) -> list[dict[str, Any]]:
    source = Path(path) if path is not None else DEFAULT_SENSITIVITY_SPEC_PATH
    with source.open("r", encoding="utf-8-sig", newline="") as f:
        return [_normalise_spec(row) for row in csv.DictReader(f)]


def run_one_way_sensitivity(
    base_config: dict[str, Any],
    economics_config: dict[str, Any] | None,
    sensitivity_specs: list[dict[str, Any]],
    outcomes: list[str],
) -> dict[str, Any]:
    rows = []
    raw_outputs = []
    for spec in sensitivity_specs:
        parameter_id = spec.get("parameterId") or spec.get("id")
        if not _has_explicit_range(spec):
            rows.append(_unresolved_parameter_row(spec, outcomes))
            continue
        values = {
            "low": spec["lowValue"],
            "base": _base_value(base_config, spec, economics_config),
            "high": spec["highValue"],
        }
        outcome_rows = {}
        for level, value in values.items():
            try:
                comparison = _run_parameter_value(base_config, economics_config, spec, value, level)
                summary = comparison["scenarioSummaries"][0]
                raw_outputs.append(
                    {
                        "parameterId": parameter_id,
                        "level": level,
                        "value": value,
                        "comparison": comparison,
                    }
                )
                outcome_rows[level] = {outcome: summary.get(outcome) for outcome in outcomes}
            except Exception as exc:
                outcome_rows[level] = {outcome: None for outcome in outcomes}
                outcome_rows[level]["failure"] = str(exc)
        for outcome in outcomes:
            rows.append(
                {
                    "parameterId": parameter_id,
                    "label": spec.get("label", parameter_id),
                    "outcome": outcome,
                    "lowValue": values["low"],
                    "baseValue": values["base"],
                    "highValue": values["high"],
                    "lowOutcome": outcome_rows.get("low", {}).get(outcome),
                    "baseOutcome": outcome_rows.get("base", {}).get(outcome),
                    "highOutcome": outcome_rows.get("high", {}).get(outcome),
                    "lowMinusBase": _subtract(outcome_rows.get("low", {}).get(outcome), outcome_rows.get("base", {}).get(outcome)),
                    "highMinusBase": _subtract(outcome_rows.get("high", {}).get(outcome), outcome_rows.get("base", {}).get(outcome)),
                    "complete": all(outcome_rows.get(level, {}).get(outcome) is not None for level in ["low", "base", "high"]),
                    "reviewStatus": spec.get("reviewStatus"),
                    "provisional": spec.get("provisional"),
                    "source": spec.get("source"),
                    "modelRecalibrationRule": spec.get("modelRecalibrationRule", ""),
                }
            )
    tornado_rows = _tornado_rows(rows)
    validation = {
        "isValid": True,
        "warnings": [
            {"parameterId": row["parameterId"], "message": row["unresolvedReason"]}
            for row in rows
            if row.get("unresolvedReason")
        ],
        "errors": [],
    }
    return {
        "contractVersion": SENSITIVITY_CONTRACT_VERSION,
        "decisionAnalysisContractVersion": DECISION_ANALYSIS_CONTRACT_VERSION,
        "analysisType": "one_way_sensitivity",
        "outcomes": outcomes,
        "specifications": sensitivity_specs,
        "results": rows,
        "tornadoRows": tornado_rows,
        "rawOutputs": raw_outputs,
        "validation": validation,
    }


def run_threshold_analysis(
    base_config: dict[str, Any],
    economics_config: dict[str, Any] | None,
    parameter_spec: dict[str, Any],
    decision_metric: str,
    search_bounds: dict[str, Any],
) -> dict[str, Any]:
    if decision_metric == "nmb" and not _valid_threshold(economics_config):
        return _threshold_unavailable(parameter_spec, decision_metric, "WTP threshold is unresolved.")
    low = float(search_bounds["low"])
    high = float(search_bounds["high"])
    grid_points = int(search_bounds.get("gridPoints", 11))
    if grid_points < 3:
        raise ValueError("Threshold analysis requires at least three grid points.")
    values = np.linspace(low, high, grid_points)
    target = float(search_bounds.get("target", 0.0))
    grid = []
    for value in values:
        metric = _evaluate_threshold_metric(
            base_config,
            economics_config,
            parameter_spec,
            decision_metric,
            float(value),
        )
        grid.append(
            {
                "parameterValue": float(value),
                "metric": metric,
                "target": target,
                "difference": None if metric is None else float(metric) - target,
                "complete": metric is not None,
            }
        )
    crossings = _detect_crossings(grid)
    roots = []
    for crossing in crossings:
        if not crossing.get("bracketed"):
            continue
        root = _bisect_root(
            base_config,
            economics_config,
            parameter_spec,
            decision_metric,
            target,
            crossing["low"],
            crossing["high"],
            float(search_bounds.get("tolerance", 1e-5)),
            int(search_bounds.get("maxIterations", 40)),
        )
        roots.append(root)
    return {
        "contractVersion": THRESHOLD_CONTRACT_VERSION,
        "analysisType": "threshold_analysis",
        "parameterSpec": parameter_spec,
        "decisionMetric": decision_metric,
        "searchBounds": search_bounds,
        "grid": grid,
        "monotonicity": monotonicity_diagnostic(grid),
        "crossingCount": len(crossings),
        "crossings": crossings,
        "roots": roots,
        "validation": {"isValid": True, "errors": [], "warnings": []},
    }


def monotonicity_diagnostic(grid: list[dict[str, Any]]) -> dict[str, Any]:
    if any(not row.get("complete") for row in grid):
        return {
            "isMonotonic": None,
            "direction": "gapped_grid",
            "hasGaps": True,
            "message": "Monotonicity was not assessed across incomplete grid points.",
        }
    values = [row["difference"] for row in grid if row.get("complete")]
    if len(values) < 2:
        return {"isMonotonic": None, "direction": "insufficient_complete_grid"}
    diffs = np.diff(values)
    nondec = bool(np.all(diffs >= -1e-10))
    noninc = bool(np.all(diffs <= 1e-10))
    if nondec and not noninc:
        return {"isMonotonic": True, "direction": "nondecreasing"}
    if noninc and not nondec:
        return {"isMonotonic": True, "direction": "nonincreasing"}
    if nondec and noninc:
        return {"isMonotonic": True, "direction": "flat"}
    return {"isMonotonic": False, "direction": "nonmonotonic"}


def crossings_from_evaluated_grid(grid: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return _detect_crossings(grid)


def _evaluate_threshold_metric(
    base_config: dict[str, Any],
    economics_config: dict[str, Any] | None,
    spec: dict[str, Any],
    metric: str,
    value: float,
) -> float | None:
    comparison = _run_parameter_value(base_config, economics_config, spec, value, f"{value:g}")
    summary = comparison["scenarioSummaries"][0]
    metric_map = {
        "positive_nmb": "nmb",
        "nmb": "nmb",
        "active_tb_cases_prevented": "active_tb_cases_prevented",
        "infections_effectively_treated": "infection_effectively_treated_total",
        "dalys_averted": "dalysAverted",
        "incremental_cost": "incrementalCost",
    }
    value = summary.get(metric_map.get(metric, metric))
    return None if value is None or pd.isna(value) else float(value)


def _bisect_root(
    base_config: dict[str, Any],
    economics_config: dict[str, Any] | None,
    spec: dict[str, Any],
    metric: str,
    target: float,
    low: float,
    high: float,
    tolerance: float,
    max_iterations: int,
) -> dict[str, Any]:
    f_low = _evaluate_threshold_metric(base_config, economics_config, spec, metric, low)
    f_high = _evaluate_threshold_metric(base_config, economics_config, spec, metric, high)
    if f_low is None or f_high is None:
        return {"root": None, "low": low, "high": high, "status": "incomplete_bracket"}
    a = low
    b = high
    fa = f_low - target
    fb = f_high - target
    if fa == 0:
        return {"root": a, "low": low, "high": high, "status": "exact_grid_point"}
    if fb == 0:
        return {"root": b, "low": low, "high": high, "status": "exact_grid_point"}
    for _ in range(max_iterations):
        mid = (a + b) / 2.0
        f_mid_raw = _evaluate_threshold_metric(base_config, economics_config, spec, metric, mid)
        if f_mid_raw is None:
            return {"root": None, "low": low, "high": high, "status": "incomplete_midpoint"}
        fm = f_mid_raw - target
        if abs(fm) <= tolerance or abs(b - a) <= tolerance:
            return {"root": mid, "low": low, "high": high, "status": "converged"}
        if fa * fm <= 0:
            b = mid
            fb = fm
        else:
            a = mid
            fa = fm
    return {"root": (a + b) / 2.0, "low": low, "high": high, "status": "max_iterations"}


def _detect_crossings(grid: list[dict[str, Any]]) -> list[dict[str, Any]]:
    crossings = []
    for left, right in zip(grid, grid[1:]):
        if not left.get("complete") or not right.get("complete"):
            continue
        a = left.get("difference")
        b = right.get("difference")
        if a is None or b is None:
            continue
        if a == 0:
            crossings.append({"type": "exact", "root": left["parameterValue"], "bracketed": False})
        if a * b < 0:
            crossings.append(
                {
                    "type": "bracket",
                    "low": left["parameterValue"],
                    "high": right["parameterValue"],
                    "bracketed": True,
                }
            )
    if grid and grid[-1].get("complete") and grid[-1].get("difference") == 0:
        crossings.append({"type": "exact", "root": grid[-1]["parameterValue"], "bracketed": False})
    return crossings


def _tornado_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for row in rows:
        if not row.get("complete"):
            continue
        lows = [row.get("lowOutcome"), row.get("highOutcome")]
        numeric = [float(v) for v in lows if v is not None and not pd.isna(v)]
        if not numeric:
            continue
        out.append(
            {
                "parameterId": row["parameterId"],
                "label": row["label"],
                "outcome": row["outcome"],
                "baseOutcome": row["baseOutcome"],
                "lowOutcome": row["lowOutcome"],
                "highOutcome": row["highOutcome"],
                "swing": max(numeric) - min(numeric),
            }
        )
    return sorted(out, key=lambda item: abs(float(item["swing"])), reverse=True)


def _run_parameter_value(
    base_config: dict[str, Any],
    economics_config: dict[str, Any] | None,
    spec: dict[str, Any],
    value: Any,
    level: str,
) -> dict[str, Any]:
    adapter = str(spec.get("adapter") or spec.get("parameterId"))
    econ = deepcopy(economics_config) if economics_config is not None else None
    scenario_changes: dict[str, Any] = {}
    if adapter in _ECONOMIC_COST_ADAPTERS:
        if econ is None:
            econ = {}
        _ECONOMIC_COST_ADAPTERS[adapter](econ, value, spec)
    else:
        scenario_changes[adapter] = value
    scenario = {
        "scenarioId": f"{spec.get('parameterId', 'parameter')}_{level}",
        "label": f"{spec.get('label', spec.get('parameterId', 'parameter'))} {level}",
        "changes": scenario_changes,
    }
    return run_scenario_comparison(
        base_config,
        econ,
        [scenario],
        model_type="expected_value",
    )


def _set_cost_sensitivity(
    econ: dict[str, Any],
    cost_item_id: str,
    value: Any,
    spec: dict[str, Any],
) -> None:
    basis = str(spec.get("monetaryValueBasis") or "target_year_converted_cost")
    for item in econ.get("costItems") or []:
        if item.get("costItemId") == cost_item_id:
            if basis == "target_year_converted_cost":
                item["targetYearCostSensitivityOverride"] = {
                    "active": True,
                    "value": float(value),
                    "unit": item.get("targetCurrency"),
                    "priceYear": item.get("targetPriceYear"),
                    "basis": basis,
                    "source": spec.get("source", ""),
                    "notes": "Sensitivity override of authoritative target-year cost; original source record preserved.",
                }
            elif basis == "original_source_year_cost":
                factor = item.get("inflationFactor")
                if factor in (None, ""):
                    raise ValueError(f"Cost item {cost_item_id} lacks a valid stored conversion factor.")
                item["sourceYearCostSensitivityOverride"] = {
                    "active": True,
                    "value": float(value),
                    "currency": item.get("originalCurrency"),
                    "priceYear": item.get("originalPriceYear"),
                    "basis": basis,
                    "source": spec.get("source", ""),
                }
                item["originalCost"] = float(value)
                item["convertedTargetYearCost"] = float(value) * float(factor)
                item["conversionStatus"] = "valid"
            else:
                raise ValueError(f"Unsupported monetaryValueBasis for {cost_item_id}: {basis}")
            return
    raise ValueError(f"Cost item {cost_item_id} is not available for sensitivity analysis.")


_ECONOMIC_COST_ADAPTERS: dict[str, Callable[[dict[str, Any], Any, dict[str, Any]], None]] = {
    "programRunningCost": lambda econ, value, spec: _set_cost_sensitivity(econ, "program_running", value, spec),
    "programSetupCost": lambda econ, value, spec: _set_cost_sensitivity(econ, "program_setup", value, spec),
    "testIGRACost": lambda econ, value, spec: _set_cost_sensitivity(econ, "test_igra", value, spec),
    "testTSTCost": lambda econ, value, spec: _set_cost_sensitivity(econ, "test_tst", value, spec),
    "activeTBDiseaseCost": lambda econ, value, spec: _set_cost_sensitivity(econ, "active_tb_disease", value, spec),
}


def _base_value(base_config: dict[str, Any], spec: dict[str, Any], economics_config: dict[str, Any] | None = None) -> Any:
    if spec.get("baseValue") not in (None, ""):
        return spec["baseValue"]
    adapter = str(spec.get("adapter") or spec.get("parameterId"))
    if adapter in _ECONOMIC_COST_ADAPTERS and economics_config:
        cost_item_id = {
            "programRunningCost": "program_running",
            "programSetupCost": "program_setup",
            "testIGRACost": "test_igra",
            "testTSTCost": "test_tst",
            "activeTBDiseaseCost": "active_tb_disease",
        }[adapter]
        for item in economics_config.get("costItems") or []:
            if item.get("costItemId") == cost_item_id:
                basis = str(spec.get("monetaryValueBasis") or "target_year_converted_cost")
                if basis == "original_source_year_cost":
                    return item.get("originalCost")
                return item.get("convertedTargetYearCost")
    if adapter == "test":
        return base_config.get("testType")
    if adapter == "regimen":
        return base_config.get("regimen")
    return base_config.get(adapter)


def _has_explicit_range(spec: dict[str, Any]) -> bool:
    return spec.get("lowValue") not in (None, "") and spec.get("highValue") not in (None, "")


def _unresolved_parameter_row(spec: dict[str, Any], outcomes: list[str]) -> dict[str, Any]:
    return {
        "parameterId": spec.get("parameterId") or spec.get("id"),
        "label": spec.get("label", spec.get("parameterId")),
        "outcome": ";".join(outcomes),
        "complete": False,
        "reviewStatus": spec.get("reviewStatus"),
        "provisional": spec.get("provisional"),
        "source": spec.get("source"),
        "unresolvedReason": spec.get("unresolvedReason") or "Explicit low/high values are required.",
    }


def _valid_threshold(economics_config: dict[str, Any] | None) -> bool:
    threshold = (economics_config or {}).get("threshold") or {}
    return threshold.get("value") not in (None, "", []) and threshold.get("source") and threshold.get("status") in {
        "configured_reviewed",
        "model_derived_reviewed",
    }


def _threshold_unavailable(spec: dict[str, Any], metric: str, reason: str) -> dict[str, Any]:
    return {
        "contractVersion": THRESHOLD_CONTRACT_VERSION,
        "analysisType": "threshold_analysis",
        "parameterSpec": spec,
        "decisionMetric": metric,
        "grid": [],
        "crossingCount": 0,
        "crossings": [],
        "roots": [],
        "validation": {"isValid": False, "errors": [], "warnings": [{"message": reason}]},
    }


def _normalise_spec(row: dict[str, Any]) -> dict[str, Any]:
    out = dict(row)
    for field in ["baseValue", "lowValue", "highValue"]:
        out[field] = _parse_optional_float(out.get(field))
    out["provisional"] = str(out.get("provisional", "")).lower() == "true"
    return out


def _parse_optional_float(value: Any) -> Any:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def _subtract(a: Any, b: Any) -> float | None:
    if a is None or b is None:
        return None
    try:
        if pd.isna(a) or pd.isna(b):
            return None
        return float(a) - float(b)
    except (TypeError, ValueError):
        return None
