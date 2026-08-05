from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

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


def run_do_nothing_summary(results: dict[str, Any]) -> dict[str, Any]:
    raw = _require_raw_dataframe(results)
    n = _population_size(results)
    ledger = _ledger_active_tb_frame(results)
    if ledger is not None:
        comparator_horizon = ledger["comparator_active_tb_cases"]
        post_active_horizon = ledger["intervention_active_tb_cases"]
        prevented_horizon = ledger["active_tb_cases_prevented"]
    else:
        comparator_horizon = raw["nActiveBy20y"]
        post_active_horizon = raw["nActiveBy20y"] - raw["nPreventedActiveTB"]
        prevented_horizon = raw["nPreventedActiveTB"]

    rel_reduction_horizon = safe_fraction_vector(
        prevented_horizon,
        comparator_horizon,
    )
    ltbi_prev = safe_fraction_vector(raw["nInfected"], n)
    active_prev_2y = safe_fraction_vector(raw["nActiveBy2y"], n)
    active_prev_horizon = safe_fraction_vector(comparator_horizon, n)
    post_active_prev_horizon = safe_fraction_vector(post_active_horizon, n)

    derived = pd.DataFrame(
        {
            "nInfected": raw["nInfected"],
            "nActiveBy2y_DoNothing": raw["nActiveBy2y"],
            "nActiveBy20y_DoNothing": comparator_horizon,
            "nActiveBy20y_AfterStrategy": post_active_horizon,
            "nActiveBy20y_Prevented": prevented_horizon,
            "relReduction20y": rel_reduction_horizon,
            "ltbiPrev_DoNothing": ltbi_prev,
            "activePrev2y_DoNothing": active_prev_2y,
            "activePrev20y_DoNothing": active_prev_horizon,
            "activePrev20y_AfterStrategy": post_active_prev_horizon,
            "nActiveByHorizon_DoNothing": comparator_horizon,
            "nActiveByHorizon_AfterStrategy": post_active_horizon,
            "nActiveByHorizon_Prevented": prevented_horizon,
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


def _ledger_active_tb_frame(results: dict[str, Any]) -> pd.DataFrame | None:
    ledger = results.get("eventLedger") or {}
    totals = ledger.get("replicateTotals")
    if not isinstance(totals, pd.DataFrame) or totals.empty:
        return None
    subset = totals[totals["eventName"].isin(["active_tb_cases", "active_tb_cases_prevented"])]
    if subset.empty:
        return None
    wide = subset.pivot_table(
        index=["replicateId", "arm"],
        columns="eventName",
        values="value",
        aggfunc="first",
    ).reset_index()
    comparator = wide[wide["arm"] == "comparator"][["replicateId", "active_tb_cases"]].rename(
        columns={"active_tb_cases": "comparator_active_tb_cases"}
    )
    intervention = wide[wide["arm"] == "intervention"][
        ["replicateId", "active_tb_cases", "active_tb_cases_prevented"]
    ].rename(columns={"active_tb_cases": "intervention_active_tb_cases"})
    merged = comparator.merge(intervention, on="replicateId", how="inner")
    if merged.empty:
        return None
    return merged
