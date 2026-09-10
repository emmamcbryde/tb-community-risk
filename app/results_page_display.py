from __future__ import annotations

from typing import Any


INTERVAL_COLUMNS = ("Median", "Low 95%", "High 95%")

KEY_METRIC_LABELS = {
    "nScreened",
    "nTestPositiveNonActive",
    "nTotalCoursesStarted",
    "nTotalCoursesCompleted",
    "nCuredInfection",
    "nPreventedActiveTB",
    "cumulative_baseline_active_tb_cases",
    "cumulative_intervention_active_tb_cases",
    "cumulative_cases_averted",
    "relative_reduction_cumulative_active_tb_cases",
}

KEY_METRIC_ORDER = (
    "nScreened",
    "nTestPositiveNonActive",
    "nTotalCoursesStarted",
    "nTotalCoursesCompleted",
    "nCuredInfection",
    "cumulative_baseline_active_tb_cases",
    "cumulative_intervention_active_tb_cases",
    "cumulative_cases_averted",
    "nPreventedActiveTB",
    "relative_reduction_cumulative_active_tb_cases",
)

FRIENDLY_METRIC_LABELS = {
    "nScreened": "People screened",
    "nTestPositiveNonActive": "Positive screening results",
    "nFalsePositiveTreated": "False-positive treatments",
    "nTotalCoursesStarted": "Preventive treatments started",
    "nTotalCoursesCompleted": "Preventive treatments completed",
    "nADRstop": "ADR-related stops",
    "nCuredInfection": "Infections effectively treated",
    "nPreventedActiveTB": "Active TB cases averted",
    "nActiveBy20y": "Active TB cases over follow-up",
    "NNS_cureInfection": "People screened per infection effectively treated",
    "NNS_preventActiveTB": "People screened per active TB case averted",
    "NNT_started_cureInfection": "Treatment starts per infection effectively treated",
    "NNT_started_preventActiveTB": "Treatment starts per active TB case averted",
    "cumulative_baseline_active_tb_cases": "Comparator active TB",
    "cumulative_intervention_active_tb_cases": "Intervention active TB",
    "cumulative_cases_averted": "Active TB cases averted",
    "relative_reduction_cumulative_active_tb_cases": "Relative reduction in active TB",
}


def results_rows_for_display(rows: list[dict[str, Any]] | None) -> list[dict[str, Any]]:
    display_rows: list[dict[str, Any]] = []
    for row in rows or []:
        metric = str(row.get("Metric", row.get("metric", "")))
        label = FRIENDLY_METRIC_LABELS.get(metric, metric)
        display_rows.append(
            {
                "Outcome": label,
                "Median": row.get("Median"),
                "Low 95%": row.get("Low95", row.get("Low 95%")),
                "High 95%": row.get("High95", row.get("High 95%")),
            }
        )
    return display_rows


def key_metric_rows_for_display(
    key_rows: list[dict[str, Any]] | None,
    dynamic_rows: list[dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    candidates = list(key_rows or []) + list(dynamic_rows or [])
    rows_by_metric = {
        str(row.get("Metric", row.get("metric", ""))): row
        for row in candidates
        if str(row.get("Metric", row.get("metric", ""))) in KEY_METRIC_LABELS
    }
    selected = []
    labels_seen = set()
    for metric in KEY_METRIC_ORDER:
        row = rows_by_metric.get(metric)
        if row is None:
            continue
        label = FRIENDLY_METRIC_LABELS.get(metric, metric)
        if label in labels_seen:
            continue
        selected.append(row)
        labels_seen.add(label)
    return results_rows_for_display(selected)


def detailed_rows_for_display(
    summary_rows: list[dict[str, Any]] | None,
    key_rows: list[dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    key_metrics = {
        str(row.get("Metric", row.get("metric", "")))
        for row in key_rows or []
    }
    details = [
        row for row in summary_rows or []
        if str(row.get("Metric", row.get("metric", ""))) not in key_metrics
    ]
    if not details:
        details = list(summary_rows or [])
    return results_rows_for_display(details)
