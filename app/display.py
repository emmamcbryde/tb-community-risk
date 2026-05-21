from __future__ import annotations

import json
import re
from datetime import datetime
from typing import Any

import pandas as pd


def arrow_safe_dataframe(data: Any) -> pd.DataFrame:
    """Return a display-only DataFrame with Arrow-safe object columns."""
    frame = pd.DataFrame(data).copy()
    for column in frame.columns:
        if frame[column].dtype == "object":
            frame[column] = frame[column].map(display_string)
    return frame


def display_string(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True)
    return str(value)


def safe_download_stem(label: object, suffix: str) -> str:
    label_text = str(label or "scenario").strip() or "scenario"
    safe_label = re.sub(r"[^A-Za-z0-9._-]+", "_", label_text).strip("._-")
    safe_label = safe_label or "scenario"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{safe_label}_{timestamp}_{suffix}"


def economics_assumptions_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    metadata = config.get("metadata", {})
    costs = config.get("costs", {})
    test_costs = costs.get("test", {})
    regimen_costs = costs.get("regimen", {})
    return [
        {"field": "currencyCode", "value": metadata.get("currencyCode")},
        {"field": "priceYear", "value": metadata.get("priceYear")},
        {"field": "locationLabel", "value": metadata.get("locationLabel")},
        {"field": "test.IGRA", "value": test_costs.get("IGRA")},
        {"field": "test.TST", "value": test_costs.get("TST")},
        {"field": "regimen.3HP", "value": regimen_costs.get("x3HP")},
        {"field": "regimen.4R", "value": regimen_costs.get("x4R")},
        {"field": "regimen.3HR", "value": regimen_costs.get("x3HR")},
        {"field": "regimen.6H", "value": regimen_costs.get("x6H")},
        {"field": "regimen.9H", "value": regimen_costs.get("x9H")},
        {
            "field": "falsePositiveIncrementalPerPerson",
            "value": costs.get("falsePositiveIncrementalPerPerson"),
        },
        {
            "field": "activeTBDiseasePerCase",
            "value": costs.get("activeTBDiseasePerCase"),
        },
        {"field": "programSetupTotal", "value": costs.get("programSetupTotal")},
        {"field": "programRunningTotal", "value": costs.get("programRunningTotal")},
    ]


def economics_summary_rows(
    economics_results_or_rows: dict[str, Any] | list[dict[str, Any]],
) -> list[dict[str, Any]]:
    if isinstance(economics_results_or_rows, dict):
        summary_rows = economics_results_or_rows.get("summaryRows") or []
    else:
        summary_rows = economics_results_or_rows or []

    preferred_columns = [
        "label",
        "value",
        "category",
        "status",
        "includedInTotal",
        "metric",
    ]
    present_columns: list[str] = []
    for column in preferred_columns:
        if any(column in row for row in summary_rows):
            present_columns.append(column)

    extra_columns: list[str] = []
    for row in summary_rows:
        for column in row:
            if column not in present_columns and column not in extra_columns:
                extra_columns.append(column)

    ordered_columns = present_columns + extra_columns
    return [
        {
            column: row.get(column)
            for column in ordered_columns
            if column in row
        }
        for row in summary_rows
    ]


def economics_summary_csv(
    economics_results_or_rows: dict[str, Any] | list[dict[str, Any]],
) -> bytes:
    rows = economics_summary_rows(economics_results_or_rows)
    return pd.DataFrame(rows).to_csv(index=False).encode("utf-8")


def economics_assumptions_json(config: dict[str, Any]) -> str:
    return json.dumps(config, indent=2, sort_keys=True)
