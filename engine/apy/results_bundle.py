from __future__ import annotations

from typing import Any

import pandas as pd

from engine.apy.summary import empirical_quantile


KEY_METRICS = [
    "nScreened",
    "nTestPositiveNonActive",
    "nFalsePositiveTreated",
    "nTotalCoursesStarted",
    "nTotalCoursesCompleted",
    "nCuredInfection",
    "nPreventedActiveTB",
    "NNS_cureInfection",
    "NNS_preventActiveTB",
    "NNT_started_cureInfection",
    "NNT_started_preventActiveTB",
]


def build_results_bundle(
    results: dict[str, Any],
    do_nothing: dict[str, Any] | None = None,
) -> dict[str, Any]:
    summary = results["summary"]
    key_metrics = summary[summary["Metric"].isin(KEY_METRICS)].copy()
    dynamic_comparison = _build_dynamic_comparison(results, do_nothing=do_nothing)
    event_ledger = results.get("eventLedger") or {}
    return {
        "metadata": {
            "available": True,
            "modelVersion": results.get("modelVersion", "python_apy_v9_port"),
            "backend": results.get("backend", "python"),
            "scenarioLabel": results.get("interfaceConfig", {}).get("scenarioLabel"),
            "contractVersion": "apy_results_bundle_v9_python_port",
        },
        "headline": {
            "available": True,
            "strategy": results["strategy"],
            "calibration": results["calibration"],
            "keyMetricsRows": _rows(key_metrics),
            "summaryRows": _rows(summary),
        },
        "technical": {
            "available": True,
            "interfaceConfig": results["interfaceConfig"],
            "calibration": results["calibration"],
            "tableMetadata": {
                "rawRows": int(len(results["raw"])),
                "rawColumns": list(results["raw"].columns),
                "summaryRows": int(len(summary)),
                "summaryColumns": list(summary.columns),
                "exampleCohortRows": _table_len(results.get("exampleCohort")),
            },
            "dynamicComparison": dynamic_comparison,
            "eventLedger": event_ledger,
        },
        "downloads": {
            "eventLedgerTotals": _rows(event_ledger.get("replicateTotals")) if event_ledger else [],
            "eventLedgerAnnual": _rows(event_ledger.get("annualEvents")) if event_ledger else [],
            "eventDefinitions": _rows(event_ledger.get("definitions")) if event_ledger else [],
        },
    }


def _build_dynamic_comparison(
    results: dict[str, Any],
    do_nothing: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if do_nothing is not None and isinstance(do_nothing.get("derived"), pd.DataFrame):
        return _build_complete_dynamic_comparison(results, do_nothing["derived"])
    return _build_partial_dynamic_comparison(results)


def _build_complete_dynamic_comparison(
    results: dict[str, Any],
    derived: pd.DataFrame,
) -> dict[str, Any]:
    cfg = results.get("interfaceConfig", {})
    metric_map = {
        "cumulative_baseline_active_tb_cases": "nActiveBy20y_DoNothing",
        "cumulative_intervention_active_tb_cases": "nActiveBy20y_AfterStrategy",
        "cumulative_cases_averted": "nActiveBy20y_Prevented",
        "relative_reduction_cumulative_active_tb_cases": "relReduction20y",
    }
    missing = [
        column for column in metric_map.values()
        if column not in derived.columns
    ]
    if missing:
        partial = _build_partial_dynamic_comparison(results)
        partial["source"] = "doNothing.derived incomplete"
        partial["missingFields"] = missing
        partial["notes"] = (
            "Do-nothing derived table was supplied but lacks required "
            "dynamic-comparison fields."
        )
        return partial

    metric_rows = [
        _dynamic_metric_row(metric, derived[column])
        for metric, column in metric_map.items()
    ]
    rows_by_metric = {row["Metric"]: row for row in metric_rows}
    return {
        "available": True,
        "source": "doNothing.derived",
        "population": cfg.get("N"),
        "followHorizon": cfg.get("followHorizon"),
        "cumulative_baseline_active_tb_cases": rows_by_metric[
            "cumulative_baseline_active_tb_cases"
        ]["Median"],
        "cumulative_intervention_active_tb_cases": rows_by_metric[
            "cumulative_intervention_active_tb_cases"
        ]["Median"],
        "cumulative_cases_averted": rows_by_metric[
            "cumulative_cases_averted"
        ]["Median"],
        "relative_reduction_cumulative_active_tb_cases": rows_by_metric[
            "relative_reduction_cumulative_active_tb_cases"
        ]["Median"],
        "metricRows": metric_rows,
        "missingFields": [],
        "notes": "",
    }


def _build_partial_dynamic_comparison(results: dict[str, Any]) -> dict[str, Any]:
    cfg = results.get("interfaceConfig", {})
    ledger_dynamic = _dynamic_from_event_ledger(results)
    if ledger_dynamic is not None:
        return ledger_dynamic
    summary_by_metric = {
        row["Metric"]: row for row in results["summary"].to_dict(orient="records")
    }
    missing = [
        "cumulative_baseline_active_tb_cases",
        "relative_reduction_cumulative_active_tb_cases",
    ]
    metric_rows = []
    if "nPreventedActiveTB" in summary_by_metric:
        row = summary_by_metric["nPreventedActiveTB"]
        metric_rows.append(
            {
                "Metric": "cumulative_cases_averted",
                "Median": row.get("Median"),
                "Low95": row.get("Low95"),
                "High95": row.get("High95"),
                "Source": "summary.nPreventedActiveTB",
                "Notes": "Python strategy-run estimate; baseline requires do-nothing parity.",
            }
        )
    if "nActiveBy20y" in summary_by_metric:
        row = summary_by_metric["nActiveBy20y"]
        metric_rows.append(
            {
                "Metric": "cumulative_intervention_active_tb_cases",
                "Median": row.get("Median"),
                "Low95": row.get("Low95"),
                "High95": row.get("High95"),
                "Source": "summary.nActiveBy20y",
            "Notes": "Legacy raw nActiveBy20y is untreated/comparator active TB; intervention active TB requires the event ledger.",
            }
        )
    return {
        "available": "partial",
        "source": "python summary without do-nothing comparator",
        "population": cfg.get("N"),
        "followHorizon": cfg.get("followHorizon"),
        "cumulative_baseline_active_tb_cases": None,
        "cumulative_intervention_active_tb_cases": None,
        "cumulative_cases_averted": _metric_median(
            summary_by_metric, "nPreventedActiveTB"
        ),
        "relative_reduction_cumulative_active_tb_cases": None,
        "metricRows": metric_rows,
        "missingFields": missing,
        "notes": "Full dynamicComparison requires do-nothing/natural-history add-on parity.",
    }


def _dynamic_from_event_ledger(results: dict[str, Any]) -> dict[str, Any] | None:
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
    comparator = wide[wide["arm"] == "comparator"]["active_tb_cases"]
    intervention = wide[wide["arm"] == "intervention"]["active_tb_cases"]
    prevented = wide[wide["arm"] == "intervention"]["active_tb_cases_prevented"]
    metric_rows = [
        _dynamic_metric_row("cumulative_baseline_active_tb_cases", comparator),
        _dynamic_metric_row("cumulative_intervention_active_tb_cases", intervention),
        _dynamic_metric_row("cumulative_cases_averted", prevented),
    ]
    rel = prevented.reset_index(drop=True) / comparator.reset_index(drop=True).replace(0, pd.NA)
    metric_rows.append(_dynamic_metric_row("relative_reduction_cumulative_active_tb_cases", rel))
    rows_by_metric = {row["Metric"]: row for row in metric_rows}
    displayed_intervention = rows_by_metric["cumulative_intervention_active_tb_cases"]["Median"]
    displayed_prevented = rows_by_metric["cumulative_cases_averted"]["Median"]
    displayed_baseline = (
        None
        if displayed_intervention is None or displayed_prevented is None
        else displayed_intervention + displayed_prevented
    )
    return {
        "available": True,
        "source": "technical.eventLedger",
        "population": results.get("interfaceConfig", {}).get("N"),
        "followHorizon": results.get("interfaceConfig", {}).get("followUpHorizonYears", results.get("interfaceConfig", {}).get("followHorizon")),
        "cumulative_baseline_active_tb_cases": displayed_baseline,
        "cumulative_intervention_active_tb_cases": displayed_intervention,
        "cumulative_cases_averted": displayed_prevented,
        "relative_reduction_cumulative_active_tb_cases": rows_by_metric["relative_reduction_cumulative_active_tb_cases"]["Median"],
        "metricRows": metric_rows,
        "missingFields": [],
        "notes": (
            "Comparator, intervention and prevented active-TB cases are derived "
            "from paired event-ledger rows. Top-level displayed values preserve "
            "comparator = intervention + prevented; metricRows retain separate "
            "distribution summaries."
        ),
    }


def _dynamic_metric_row(metric: str, values) -> dict[str, Any]:
    x = pd.Series(values, dtype="float64")
    x = x.dropna()
    if x.empty:
        median = low = high = None
    else:
        median = float(x.median())
        low = float(empirical_quantile(x.to_numpy(), 0.025))
        high = float(empirical_quantile(x.to_numpy(), 0.975))
    return {
        "Metric": metric,
        "Median": median,
        "Low95": low,
        "High95": high,
        "Source": "doNothing.derived",
        "Notes": "",
    }


def _metric_median(summary_by_metric: dict[str, dict[str, Any]], metric: str):
    row = summary_by_metric.get(metric)
    return None if row is None else row.get("Median")


def _rows(df: pd.DataFrame) -> list[dict[str, Any]]:
    if df is None:
        return []
    if isinstance(df, list):
        return df
    return df.to_dict(orient="records")


def _table_len(value) -> int | None:
    return None if value is None else int(len(value))
