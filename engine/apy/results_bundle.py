from __future__ import annotations

from typing import Any

import pandas as pd


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


def build_results_bundle(results: dict[str, Any]) -> dict[str, Any]:
    summary = results["summary"]
    key_metrics = summary[summary["Metric"].isin(KEY_METRICS)].copy()
    dynamic_comparison = _build_dynamic_comparison(results)
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
        },
        "downloads": {},
    }


def _build_dynamic_comparison(results: dict[str, Any]) -> dict[str, Any]:
    cfg = results.get("interfaceConfig", {})
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
                "Notes": "Untreated active TB by followHorizon in the strategy cohort.",
            }
        )
    return {
        "available": "partial",
        "source": "python summary without do-nothing comparator",
        "population": cfg.get("N"),
        "followHorizon": cfg.get("followHorizon"),
        "cumulative_baseline_active_tb_cases": None,
        "cumulative_intervention_active_tb_cases": _metric_median(
            summary_by_metric, "nActiveBy20y"
        ),
        "cumulative_cases_averted": _metric_median(
            summary_by_metric, "nPreventedActiveTB"
        ),
        "relative_reduction_cumulative_active_tb_cases": None,
        "metricRows": metric_rows,
        "missingFields": missing,
        "notes": "Full dynamicComparison requires do-nothing/natural-history add-on parity.",
    }


def _metric_median(summary_by_metric: dict[str, dict[str, Any]], metric: str):
    row = summary_by_metric.get(metric)
    return None if row is None else row.get("Median")


def _rows(df: pd.DataFrame) -> list[dict[str, Any]]:
    return df.to_dict(orient="records")


def _table_len(value) -> int | None:
    return None if value is None else int(len(value))
