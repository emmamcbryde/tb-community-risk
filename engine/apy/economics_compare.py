from __future__ import annotations

from typing import Any, Mapping


ECONOMICS_VALUE_COLUMNS = ("Value", "value")
ECONOMICS_KEY_FIELDS = ("Metric", "metric", "component")
METADATA_FIELDS = ("category", "status", "includedInTotal")


def coerce_number(value: Any) -> float | None:
    try:
        if value in (None, ""):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def economics_rows(economics: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not economics:
        return []
    rows = economics.get("summaryRows") or []
    return rows if isinstance(rows, list) else []


def compare_economics_rows(
    baseline_rows: list[dict[str, Any]],
    comparator_rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[str]]:
    return _compare_metric_rows(
        baseline_rows,
        comparator_rows,
        ECONOMICS_VALUE_COLUMNS,
        ECONOMICS_KEY_FIELDS,
        require_value=False,
    )


def _first_present(row: Mapping[str, Any], fields: tuple[str, ...]) -> Any:
    for field in fields:
        if field in row:
            return row.get(field)
    return None


def _index_metric_rows(
    rows: list[dict[str, Any]],
    value_columns: tuple[str, ...],
    key_fields: tuple[str, ...],
    require_value: bool,
) -> tuple[dict[str, dict[str, Any]], list[str], list[str], int]:
    values: dict[str, dict[str, Any]] = {}
    order: list[str] = []
    duplicates: list[str] = []
    skipped = 0
    for row in rows:
        metric = _first_present(row, key_fields)
        value = _first_present(row, value_columns)
        if metric in (None, "") or (require_value and value is None):
            skipped += 1
            continue
        metric_key = str(metric)
        if metric_key not in values:
            order.append(metric_key)
            values[metric_key] = {
                "value": value,
                "row": row,
            }
        else:
            duplicates.append(metric_key)
    return values, order, duplicates, skipped


def _compare_metric_rows(
    baseline_rows: list[dict[str, Any]],
    comparator_rows: list[dict[str, Any]],
    value_columns: tuple[str, ...],
    key_fields: tuple[str, ...],
    require_value: bool,
) -> tuple[list[dict[str, Any]], list[str]]:
    baseline, baseline_order, baseline_duplicates, baseline_skipped = _index_metric_rows(
        baseline_rows,
        value_columns,
        key_fields,
        require_value,
    )
    comparator, comparator_order, comparator_duplicates, comparator_skipped = _index_metric_rows(
        comparator_rows,
        value_columns,
        key_fields,
        require_value,
    )

    ordered_metrics = baseline_order + [
        metric for metric in comparator_order if metric not in baseline
    ]
    rows: list[dict[str, Any]] = []
    for metric in ordered_metrics:
        baseline_item = baseline.get(metric, {})
        comparator_item = comparator.get(metric, {})
        base_value = baseline_item.get("value")
        comp_value = comparator_item.get("value")
        base_num = coerce_number(base_value)
        comp_num = coerce_number(comp_value)
        abs_diff = None
        rel_diff = None
        if base_num is not None and comp_num is not None:
            abs_diff = comp_num - base_num
            if base_num != 0:
                rel_diff = abs_diff / base_num

        diff_row = {
            "metric": metric,
            "baseline": base_value,
            "comparator": comp_value,
            "absoluteDifference": abs_diff,
            "relativeDifference": rel_diff,
        }
        baseline_source = baseline_item.get("row")
        comparator_source = comparator_item.get("row")
        for metadata_field in METADATA_FIELDS:
            if isinstance(baseline_source, dict) and metadata_field in baseline_source:
                diff_row[metadata_field] = baseline_source.get(metadata_field)
            elif isinstance(comparator_source, dict) and metadata_field in comparator_source:
                diff_row[metadata_field] = comparator_source.get(metadata_field)
        rows.append(diff_row)

    warnings: list[str] = []
    duplicate_metrics = sorted(set(baseline_duplicates + comparator_duplicates))
    if duplicate_metrics:
        warnings.append(
            "Duplicate Metric values detected; the first row was used for: "
            + ", ".join(duplicate_metrics)
        )
    if baseline_skipped or comparator_skipped:
        value_column_label = "/".join(value_columns)
        key_label = "/".join(key_fields)
        warnings.append(
            f"Skipped rows missing {key_label} + {value_column_label}: "
            f"baseline={baseline_skipped}, comparator={comparator_skipped}."
        )
    return rows, warnings


__all__ = [
    "compare_economics_rows",
    "coerce_number",
    "economics_rows",
]
