from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd


CONTRACT_VERSION = "dynamic_output_contract_v9"
MODEL_VERSION = "dynamic_python_v1"


def json_like(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, dict):
        return {str(key): json_like(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_like(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_like(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, (str, int, bool)):
        return value
    return str(value)


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


def rows_from_dataframe(df_future: Any) -> list[dict[str, Any]]:
    if df_future is None:
        return []
    frame = pd.DataFrame(df_future).copy()
    if frame.empty:
        return []
    return [json_like(row) for row in frame.to_dict(orient="records")]


def last_value(frame: pd.DataFrame, column: str) -> float | int | None:
    if column not in frame.columns or frame.empty:
        return None
    series = frame[column].dropna()
    if series.empty:
        return None
    return number_or_none(series.iloc[-1])


def projection_horizon(frame: pd.DataFrame) -> float | int | None:
    if "Year" in frame.columns:
        return last_value(frame, "Year")
    if frame.empty:
        return None
    return int(len(frame))


def relative_reduction(baseline: Any, intervention: Any) -> float | None:
    baseline_value = number_or_none(baseline)
    intervention_value = number_or_none(intervention)
    if baseline_value in (None, 0) or intervention_value is None:
        return None
    return float((baseline_value - intervention_value) / baseline_value)


def metric_row(metric: str, value: Any, notes: str = "") -> dict[str, Any]:
    return {
        "Metric": metric,
        "Value": json_like(value),
        "Notes": notes,
    }


def build_dynamic_results_bundle_v9(
    df_future: Any,
    params_base: dict[str, Any] | None = None,
    params_intervention: dict[str, Any] | None = None,
    calibration: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    frame = pd.DataFrame(df_future).copy() if df_future is not None else pd.DataFrame()
    annual_rows = rows_from_dataframe(frame)

    horizon = projection_horizon(frame)
    baseline_cum = last_value(frame, "Baseline_cum_count")
    intervention_cum = last_value(frame, "Intervention_cum_count")
    cases_averted_cum = last_value(frame, "Cases_averted_cum_count")
    if cases_averted_cum is None and baseline_cum is not None and intervention_cum is not None:
        cases_averted_cum = baseline_cum - intervention_cum

    summary_rows = [
        metric_row("projection_horizon", horizon),
        metric_row("final_year_baseline_incidence_per100k", last_value(frame, "Baseline_inc_per100k")),
        metric_row("final_year_intervention_incidence_per100k", last_value(frame, "Intervention_inc_per100k")),
        metric_row("cumulative_baseline_active_tb_cases", baseline_cum),
        metric_row("cumulative_intervention_active_tb_cases", intervention_cum),
        metric_row("cumulative_cases_averted", cases_averted_cum),
        metric_row("relative_reduction_cumulative_active_tb_cases", relative_reduction(baseline_cum, intervention_cum)),
    ]
    key_metrics = [row for row in summary_rows if row["Value"] is not None]

    return {
        "model": "Dynamic TB model",
        "modelVersion": MODEL_VERSION,
        "contractVersion": CONTRACT_VERSION,
        "metadata": json_like(metadata or {}),
        "headline": {
            "summaryRows": summary_rows,
            "keyMetricsRows": key_metrics,
        },
        "projection": {
            "annualRows": annual_rows,
        },
        "technical": {
            "paramsBase": json_like(params_base or {}),
            "paramsIntervention": json_like(params_intervention or {}),
            "calibration": json_like(calibration or {}),
        },
    }
