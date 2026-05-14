from __future__ import annotations

import math
from typing import Iterable

import pandas as pd


def compare_python_summary_to_reference(
    python_summary: pd.DataFrame,
    reference_summary: pd.DataFrame,
    metrics: list[str] | None = None,
) -> pd.DataFrame:
    python_rows = _summary_by_metric(python_summary)
    reference_rows = _summary_by_metric(reference_summary)
    metric_names = metrics or sorted(set(python_rows) | set(reference_rows))
    rows = []
    for metric in metric_names:
        py_value = _median_or_none(python_rows.get(metric))
        matlab_value = _median_or_none(reference_rows.get(metric))
        comparable = py_value is not None and matlab_value is not None
        if comparable:
            abs_diff = py_value - matlab_value
            rel_diff = math.nan if matlab_value == 0 else abs_diff / matlab_value
            notes = "Diagnostic stochastic comparison; strict equality is not expected."
        else:
            abs_diff = math.nan
            rel_diff = math.nan
            notes = "Metric missing from Python or MATLAB summary."
        rows.append(
            {
                "Metric": metric,
                "PythonMedian": py_value,
                "MatlabMedian": matlab_value,
                "AbsoluteDifference": abs_diff,
                "RelativeDifference": rel_diff,
                "Comparable": comparable,
                "Notes": notes,
            }
        )
    return pd.DataFrame(rows)


def _summary_by_metric(summary: pd.DataFrame) -> dict[str, dict]:
    if summary is None or summary.empty or "Metric" not in summary.columns:
        return {}
    return {
        str(row["Metric"]): row
        for row in summary.to_dict(orient="records")
    }


def _median_or_none(row: dict | None):
    if row is None or "Median" not in row:
        return None
    value = row["Median"]
    if pd.isna(value):
        return math.nan
    return float(value)
