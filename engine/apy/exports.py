from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


UNSUPPORTED_SUMMARY_EXPORT_MESSAGE = (
    "APY summary CSV export is unavailable because headline.summaryRows is absent."
)
UNSUPPORTED_HEADLINE_KEY_METRICS_MESSAGE = (
    "APY headline key metrics table is unavailable because headline.keyMetricsRows is absent."
)
UNSUPPORTED_HEADLINE_SUMMARY_MESSAGE = (
    "APY headline summary table is unavailable because headline.summaryRows is absent."
)
UNSUPPORTED_CHART_EXPORT_MESSAGE = (
    "APY chart series are unavailable because charts.plotDataRows is absent."
)

SUMMARY_COLUMN_ORDER = ("Metric", "Median", "Low95", "High95")


def json_export_payload(bundle: Mapping[str, Any]) -> dict[str, Any]:
    """Return a deterministic JSON-serialisable copy of an APY results bundle."""
    payload = _json_value(bundle)
    if not isinstance(payload, dict):
        raise TypeError("APY export payload requires a mapping bundle.")
    return payload


def summary_rows_for_csv(bundle: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return stable CSV-ready rows from headline.summaryRows."""
    rows = _headline_summary_rows(bundle)
    if not rows:
        return []

    columns = _stable_columns(rows, preferred=SUMMARY_COLUMN_ORDER)
    return [
        {column: _json_value(row.get(column)) for column in columns}
        for row in rows
    ]


def summary_csv_export(bundle: Mapping[str, Any]) -> dict[str, Any]:
    """Return summary CSV export status plus stable rows."""
    rows = summary_rows_for_csv(bundle)
    if not rows:
        return {
            "available": False,
            "message": UNSUPPORTED_SUMMARY_EXPORT_MESSAGE,
            "columns": [],
            "rows": [],
        }
    return {
        "available": True,
        "message": "",
        "columns": list(rows[0].keys()),
        "rows": rows,
    }


def headline_display_tables(bundle: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """Return deterministic JSON-safe display tables from headline rows."""
    return {
        "keyMetricsRows": _headline_display_table(
            _headline_key_metrics_rows(bundle),
            unavailable_message=UNSUPPORTED_HEADLINE_KEY_METRICS_MESSAGE,
        ),
        "summaryRows": _headline_display_table(
            _headline_summary_rows(bundle),
            unavailable_message=UNSUPPORTED_HEADLINE_SUMMARY_MESSAGE,
        ),
    }


def chart_numeric_series(bundle: Mapping[str, Any]) -> dict[str, Any]:
    """Return generic chart-ready numeric series from charts.plotDataRows."""
    rows = _chart_plot_data_rows(bundle)
    if not rows:
        return {
            "available": False,
            "message": UNSUPPORTED_CHART_EXPORT_MESSAGE,
            "series": [],
        }

    columns = _stable_columns(rows)
    numeric_columns = [
        column
        for column in columns
        if any(_is_number(row.get(column)) for row in rows)
    ]
    series = [
        {
            "name": column,
            "points": [
                {
                    "index": index,
                    "value": _number_or_none(row.get(column)),
                    "row": _series_row_label(row, columns),
                }
                for index, row in enumerate(rows)
            ],
        }
        for column in numeric_columns
    ]

    if not series:
        return {
            "available": False,
            "message": "APY chart series are unavailable because charts.plotDataRows contains no numeric columns.",
            "series": [],
        }
    return {"available": True, "message": "", "series": series}


def _headline_display_table(
    rows: Sequence[Mapping[str, Any]],
    *,
    unavailable_message: str,
) -> dict[str, Any]:
    if not rows:
        return {
            "available": False,
            "message": unavailable_message,
            "columns": [],
            "rows": [],
        }

    columns = _stable_columns(rows, preferred=SUMMARY_COLUMN_ORDER)
    return {
        "available": True,
        "message": "",
        "columns": columns,
        "rows": [
            {column: _json_value(row.get(column)) for column in columns}
            for row in rows
        ],
    }


def _headline_key_metrics_rows(bundle: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    headline = bundle.get("headline")
    if not isinstance(headline, Mapping):
        return []
    return _mapping_rows(headline.get("keyMetricsRows"))


def _headline_summary_rows(bundle: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    headline = bundle.get("headline")
    if not isinstance(headline, Mapping):
        return []
    return _mapping_rows(headline.get("summaryRows"))


def _chart_plot_data_rows(bundle: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    charts = bundle.get("charts")
    if not isinstance(charts, Mapping):
        return []
    return _mapping_rows(charts.get("plotDataRows"))


def _mapping_rows(value: Any) -> list[Mapping[str, Any]]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [row for row in value if isinstance(row, Mapping)]
    return []


def _stable_columns(
    rows: Sequence[Mapping[str, Any]],
    preferred: Sequence[str] = (),
) -> list[str]:
    seen = {key for row in rows for key in row}
    ordered = [column for column in preferred if column in seen]
    ordered.extend(sorted(column for column in seen if column not in ordered))
    return ordered


def _series_row_label(row: Mapping[str, Any], columns: Sequence[str]) -> dict[str, Any]:
    return {
        column: _json_value(value)
        for column in columns
        if not _is_number(value := row.get(column))
    }


def _json_value(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Mapping):
        return {str(key): _json_value(value[key]) for key in sorted(value, key=str)}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_json_value(item) for item in value]
    if hasattr(value, "item"):
        return _json_value(value.item())
    return str(value)


def _is_number(value: Any) -> bool:
    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return math.isfinite(float(value))
    if hasattr(value, "item"):
        return _is_number(value.item())
    return False


def _number_or_none(value: Any) -> float | int | None:
    if not _is_number(value):
        return None
    if isinstance(value, int):
        return value
    if hasattr(value, "item"):
        return _number_or_none(value.item())
    return float(value)
