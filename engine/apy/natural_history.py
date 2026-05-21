from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from engine.apy.attributable_risk import run_attributable_risk
from engine.apy.summary import empirical_quantile


SUMMARY_METRICS = [
    (
        "LTBI prevalence at baseline (do nothing)",
        "ltbiPrev_DoNothing",
        "proportion",
    ),
    (
        "Active TB cases by 2 years (do nothing)",
        "nActiveBy2y_DoNothing",
        "count",
    ),
    (
        "Active TB prevalence by 2 years (do nothing)",
        "activePrev2y_DoNothing",
        "proportion",
    ),
    (
        "Active TB cases by 20 years (do nothing)",
        "nActiveBy20y_DoNothing",
        "count",
    ),
    (
        "Active TB prevalence by 20 years (do nothing)",
        "activePrev20y_DoNothing",
        "proportion",
    ),
    (
        "Active TB cases by 20 years after current strategy",
        "nActiveBy20y_AfterStrategy",
        "count",
    ),
    (
        "Active TB prevalence by 20 years after current strategy",
        "activePrev20y_AfterStrategy",
        "proportion",
    ),
    (
        "Active TB cases prevented by 20 years",
        "nActiveBy20y_Prevented",
        "count",
    ),
    (
        "Relative reduction in 20-year active TB burden",
        "relReduction20y",
        "proportion",
    ),
    (
        "NNS to prevent one active TB case (current strategy)",
        "NNS_preventActiveTB",
        "count",
    ),
    (
        "NNS to cure one infection (current strategy)",
        "NNS_cureInfection",
        "count",
    ),
]


REQUIRED_RAW_COLUMNS = [
    "nInfected",
    "nActiveBy2y",
    "nActiveBy20y",
    "nPreventedActiveTB",
    "NNS_preventActiveTB",
    "NNS_cureInfection",
]

_ADDON_REPORT_SOURCE = "apy_python_natural_history_addon_report"


def run_do_nothing_summary(results: dict[str, Any]) -> dict[str, Any]:
    raw = _require_raw_dataframe(results)
    n = _population_size(results)

    post_active_horizon = raw["nActiveBy20y"] - raw["nPreventedActiveTB"]
    rel_reduction_horizon = safe_fraction_vector(
        raw["nPreventedActiveTB"],
        raw["nActiveBy20y"],
    )
    ltbi_prev = safe_fraction_vector(raw["nInfected"], n)
    active_prev_2y = safe_fraction_vector(raw["nActiveBy2y"], n)
    active_prev_horizon = safe_fraction_vector(raw["nActiveBy20y"], n)
    post_active_prev_horizon = safe_fraction_vector(post_active_horizon, n)

    derived = pd.DataFrame(
        {
            "nInfected": raw["nInfected"],
            "nActiveBy2y_DoNothing": raw["nActiveBy2y"],
            "nActiveBy20y_DoNothing": raw["nActiveBy20y"],
            "nActiveBy20y_AfterStrategy": post_active_horizon,
            "nActiveBy20y_Prevented": raw["nPreventedActiveTB"],
            "relReduction20y": rel_reduction_horizon,
            "ltbiPrev_DoNothing": ltbi_prev,
            "activePrev2y_DoNothing": active_prev_2y,
            "activePrev20y_DoNothing": active_prev_horizon,
            "activePrev20y_AfterStrategy": post_active_prev_horizon,
            "nActiveByHorizon_DoNothing": raw["nActiveBy20y"],
            "nActiveByHorizon_AfterStrategy": post_active_horizon,
            "nActiveByHorizon_Prevented": raw["nPreventedActiveTB"],
            "relReductionHorizon": rel_reduction_horizon,
            "NNS_preventActiveTB": raw["NNS_preventActiveTB"],
            "NNS_cureInfection": raw["NNS_cureInfection"],
        }
    )

    summary = pd.DataFrame(
        [
            summary_row(metric, derived[column], display_scale)
            for metric, column, display_scale in SUMMARY_METRICS
        ]
    )
    return {
        "summary": summary,
        "derived": derived,
        "resultsUsed": results,
    }


def build_natural_history_addon_report(
    results: Mapping[str, Any],
    do_nothing: Mapping[str, Any] | None = None,
    requested_attributable_metrics: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Build a JSON-safe natural-history add-on report contract."""
    missing_inputs: list[dict[str, str]] = []
    messages: list[str] = []

    if do_nothing is None:
        summary_rows: list[dict[str, Any]] = []
        do_nothing_payload = {
            "status": "missing-input",
            "source": None,
            "summaryRows": [],
            "derivedRows": [],
        }
        missing_inputs.append(
            {
                "field": "doNothing",
                "message": (
                    "Do-nothing natural-history summary was not provided; "
                    "run_do_nothing_summary output is required for these rows."
                ),
            }
        )
        messages.append("Missing do-nothing natural-history summary.")
    else:
        summary_rows = _rows(do_nothing.get("summary"))
        do_nothing_payload = {
            "status": "available",
            "source": "run_do_nothing_summary",
            "summaryRows": summary_rows,
            "derivedRows": _rows(do_nothing.get("derived")),
        }
        messages.append("Do-nothing natural-history rows were included.")

    attributable_risk = run_attributable_risk(
        results,
        requested_metrics=requested_attributable_metrics,
    )
    attributable_missing = _json_like(attributable_risk.get("missingInputs", []))
    unsupported_metrics = _json_like(attributable_risk.get("unsupportedMetrics", []))
    missing_inputs.extend(attributable_missing)
    messages.extend(_json_like(attributable_risk.get("messages", [])))

    if unsupported_metrics:
        messages.append("Some attributable-risk metrics are unsupported.")

    return {
        "status": _addon_report_status(missing_inputs, unsupported_metrics),
        "source": _ADDON_REPORT_SOURCE,
        "summaryRows": summary_rows,
        "doNothing": do_nothing_payload,
        "attributableRisk": _json_like(attributable_risk),
        "missingInputs": missing_inputs,
        "unsupportedMetrics": unsupported_metrics,
        "messages": messages,
    }


def safe_fraction_vector(num, den):
    numerator = pd.Series(num, dtype="float64")
    if np.isscalar(den):
        denominator = pd.Series(float(den), index=numerator.index, dtype="float64")
    else:
        denominator = pd.Series(den, dtype="float64").reset_index(drop=True)
        if len(denominator) != len(numerator):
            denominator = denominator.reindex(numerator.index)

    with np.errstate(divide="ignore", invalid="ignore"):
        out = numerator / denominator
    out = out.mask(denominator == 0, np.nan)
    return out.to_numpy()


def summarise_vector(values) -> dict[str, float]:
    x = pd.Series(values, dtype="float64")
    x = x[np.isfinite(x)]
    if x.empty:
        return {"Median": np.nan, "Low95": np.nan, "High95": np.nan}
    return {
        "Median": float(x.median()),
        "Low95": float(empirical_quantile(x.to_numpy(), 0.025)),
        "High95": float(empirical_quantile(x.to_numpy(), 0.975)),
    }


def summary_row(metric: str, values, display_scale: str) -> dict[str, Any]:
    row = {"Metric": metric, **summarise_vector(values)}
    row["DisplayScale"] = display_scale
    return row


def _require_raw_dataframe(results: dict[str, Any]) -> pd.DataFrame:
    raw = results.get("raw")
    if not isinstance(raw, pd.DataFrame):
        raise ValueError("Results must contain a raw pandas DataFrame.")

    missing = [name for name in REQUIRED_RAW_COLUMNS if name not in raw.columns]
    if missing:
        raise ValueError(
            "Results raw table is missing required columns: "
            + ", ".join(missing)
        )
    return raw


def _population_size(results: dict[str, Any]) -> float:
    cfg = results.get("interfaceConfig", {})
    if "N" not in cfg:
        raise ValueError("Results must contain interfaceConfig.N.")
    return float(cfg["N"])


def _rows(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, pd.DataFrame):
        return _json_like(value.to_dict(orient="records"))
    if isinstance(value, list):
        return _json_like(value)
    return []


def _json_like(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_like(item) for key, item in value.items()}
    if isinstance(value, pd.DataFrame):
        return _json_like(value.to_dict(orient="records"))
    if isinstance(value, pd.Series):
        return _json_like(value.tolist())
    if isinstance(value, tuple):
        return [_json_like(item) for item in value]
    if isinstance(value, list):
        return [_json_like(item) for item in value]
    if isinstance(value, np.generic):
        return _json_like(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _addon_report_status(
    missing_inputs: list[dict[str, str]],
    unsupported_metrics: list[Any],
) -> str:
    if missing_inputs:
        return "missing-input"
    if unsupported_metrics:
        return "unsupported"
    return "available"
