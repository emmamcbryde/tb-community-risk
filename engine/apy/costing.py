from __future__ import annotations

from copy import deepcopy
import csv
import math
from pathlib import Path
from typing import Any


TARGET_CURRENCY = "AUD"
TARGET_PRICE_YEAR = "2025-26"
DEFAULT_INDEX_ID = "AUD_HEALTH_SYSTEM_CPI"
INFLATION_INDEX_PATH = Path(__file__).resolve().parents[2] / "data" / "inflation_indices.csv"


def load_inflation_indices(path: str | Path | None = None) -> list[dict[str, Any]]:
    source = Path(path) if path is not None else INFLATION_INDEX_PATH
    with source.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    return [_normalise_index_row(row) for row in rows]


def build_cost_item(
    *,
    cost_item_id: str,
    description: str,
    original_cost: Any,
    original_currency: str,
    original_price_year: str | int | None,
    source: str,
    resource_category: str,
    page_table_item: str = "",
    source_year_status: str = "unknown",
    target_currency: str = TARGET_CURRENCY,
    target_price_year: str = TARGET_PRICE_YEAR,
    inflation_index_id: str = DEFAULT_INDEX_ID,
    notes: str = "",
    resource_use: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "costItemId": cost_item_id,
        "description": description,
        "originalCost": original_cost,
        "originalCurrency": original_currency,
        "originalPriceYear": original_price_year,
        "sourceCitation": source,
        "pageTableItem": page_table_item,
        "sourceYearStatus": source_year_status,
        "resourceCategory": resource_category,
        "targetCurrency": target_currency,
        "targetPriceYear": target_price_year,
        "inflationIndexId": inflation_index_id,
        "indexSource": "",
        "indexVersion": "",
        "sourceYearIndexValue": None,
        "targetYearIndexValue": None,
        "inflationFactor": None,
        "convertedTargetYearCost": None,
        "conversionStatus": "not_converted",
        "warnings": [],
        "notes": notes,
        "resourceUse": deepcopy(resource_use or {}),
        "conversionApplied": False,
    }


def normalise_cost_table(
    cost_items: list[dict[str, Any]],
    index_rows: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    rows = index_rows if index_rows is not None else load_inflation_indices()
    lookup = {
        (str(row["index_id"]), str(row["price_year"])): row
        for row in rows
        if row.get("index_value") is not None
    }
    return [normalise_cost_item(item, lookup) for item in cost_items]


def normalise_cost_item(
    item: dict[str, Any],
    index_lookup: dict[tuple[str, str], dict[str, Any]],
) -> dict[str, Any]:
    out = deepcopy(item)
    warnings = list(out.get("warnings") or [])
    if out.get("conversionApplied"):
        warnings.append("Cost item already has conversionApplied=true; double inflation prevented.")
        out["conversionStatus"] = "invalid_double_inflation"
        out["warnings"] = warnings
        return out

    original_cost = _number_or_none(out.get("originalCost"))
    source_year = out.get("originalPriceYear")
    target_year = out.get("targetPriceYear")
    index_id = str(out.get("inflationIndexId") or "")

    if original_cost is None:
        warnings.append("Original cost is missing or non-numeric.")
        out["conversionStatus"] = "unresolved_original_cost"
    elif not source_year:
        warnings.append("Source price year is unknown; cost was not inflated.")
        out["conversionStatus"] = "unresolved_source_price_year"
    elif str(source_year) == str(target_year):
        out["sourceYearIndexValue"] = 1.0
        out["targetYearIndexValue"] = 1.0
        out["inflationFactor"] = 1.0
        out["convertedTargetYearCost"] = original_cost
        out["conversionStatus"] = "valid"
        out["conversionApplied"] = True
    else:
        source_index = index_lookup.get((index_id, str(source_year)))
        target_index = index_lookup.get((index_id, str(target_year)))
        if source_index is None or target_index is None:
            warnings.append(
                "Required inflation index value is missing; cost was not inflated."
            )
            out["conversionStatus"] = "unresolved_index_value"
        else:
            source_value = _number_or_none(source_index.get("index_value"))
            target_value = _number_or_none(target_index.get("index_value"))
            if source_value is None or target_value is None or source_value <= 0:
                warnings.append("Inflation index values must be positive numbers.")
                out["conversionStatus"] = "unresolved_index_value"
            else:
                factor = target_value / source_value
                converted = original_cost * factor
                out["indexSource"] = target_index.get("index_source", "")
                out["indexVersion"] = target_index.get("index_version", "")
                out["sourceYearIndexValue"] = source_value
                out["targetYearIndexValue"] = target_value
                out["inflationFactor"] = factor
                out["convertedTargetYearCost"] = converted
                out["conversionStatus"] = "valid"
                out["conversionApplied"] = True
    out["warnings"] = warnings
    if out.get("conversionStatus") == "valid":
        _assert_reconciles(out)
    return out


def valid_converted_cost(item: dict[str, Any]) -> float | None:
    if item.get("conversionStatus") != "valid":
        return None
    return _number_or_none(item.get("convertedTargetYearCost"))


def discount_value(value: Any, annual_rate: Any, years_from_present: Any) -> float | None:
    amount = _number_or_none(value)
    rate = _number_or_none(annual_rate)
    years = _number_or_none(years_from_present)
    if amount is None or rate is None or years is None:
        return None
    if rate < 0:
        raise ValueError("Discount rate must be non-negative.")
    return amount / ((1 + rate) ** years)


def unresolved_cost_warnings(cost_items: list[dict[str, Any]]) -> list[str]:
    messages: list[str] = []
    for item in cost_items:
        if item.get("conversionStatus") != "valid":
            messages.append(
                f"{item.get('costItemId')}: {item.get('conversionStatus')} - "
                + "; ".join(item.get("warnings") or [])
            )
    return messages


def _assert_reconciles(item: dict[str, Any]) -> None:
    original = _number_or_none(item.get("originalCost"))
    factor = _number_or_none(item.get("inflationFactor"))
    converted = _number_or_none(item.get("convertedTargetYearCost"))
    if original is None or factor is None or converted is None:
        raise ValueError("Converted cost is marked valid but is incomplete.")
    if not math.isclose(original * factor, converted, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("Converted cost does not reconcile with the recorded factor.")


def _normalise_index_row(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "index_id": row.get("index_id", ""),
        "index_source": row.get("index_source", ""),
        "index_version": row.get("index_version", ""),
        "price_year": row.get("price_year", ""),
        "index_value": _number_or_none(row.get("index_value")),
        "notes": row.get("notes", ""),
    }


def _number_or_none(value: Any) -> float | None:
    if value in (None, "", []):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number
