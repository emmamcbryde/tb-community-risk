from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from math import floor, isfinite
from typing import Any

import pandas as pd
import streamlit as st


GRID_SIZE = 100
GRID_COLUMNS = 10
ROUNDING_NOTE = (
    "Icons fill the nearest whole person out of 100; text shows the expected "
    "value to the nearest tenth per 100 eligible people."
)


@dataclass(frozen=True)
class OutcomeDefinition:
    outcome_id: str
    event_name: str
    label: str
    color: str
    intervention_only: bool = False


OUTCOME_DEFINITIONS: tuple[OutcomeDefinition, ...] = (
    OutcomeDefinition("active_tb_cases", "active_tb_cases", "Active TB cases", "#d97706"),
    OutcomeDefinition("tpt_started", "tpt_started_total", "Preventive treatments started", "#2563eb"),
    OutcomeDefinition("tpt_completed", "tpt_completed_total", "Preventive treatments completed", "#16a34a"),
    OutcomeDefinition("false_positive_treatments", "tpt_started_false_positive", "False-positive treatments", "#db2777"),
    OutcomeDefinition("adverse_events", "tpt_adr_stop_total", "Adverse events", "#dc2626"),
    OutcomeDefinition(
        "active_tb_cases_prevented",
        "active_tb_cases_prevented",
        "Active TB cases prevented",
        "#0d9488",
        intervention_only=True,
    ),
)


def icon_grid_value(value_per_100: Any) -> dict[str, Any]:
    value = _optional_float(value_per_100)
    if value is None:
        return {
            "available": False,
            "exactPer100": None,
            "displayPer100": "Unavailable",
            "filledIcons": 0,
            "roundingRule": ROUNDING_NOTE,
        }
    clipped = max(0.0, min(float(value), float(GRID_SIZE)))
    return {
        "available": True,
        "exactPer100": clipped,
        "displayPer100": f"{clipped:.1f} per 100",
        "filledIcons": int(max(0, min(GRID_SIZE, floor(clipped + 0.5)))),
        "roundingRule": ROUNDING_NOTE,
    }


def build_100_person_visual_data(
    event_ledger: dict[str, Any] | None,
    *,
    outcomes: Iterable[str] | None = None,
) -> list[dict[str, Any]]:
    if not isinstance(event_ledger, dict):
        return []
    totals = _records_to_frame(event_ledger.get("replicateTotals"))
    if totals.empty:
        return []
    selected = set(outcomes) if outcomes is not None else None
    values = _mean_values_by_arm_event(totals)
    denominators = {
        arm: values.get((arm, "eligible_population"))
        for arm in ("comparator", "intervention")
    }
    rows: list[dict[str, Any]] = []
    for definition in OUTCOME_DEFINITIONS:
        if selected is not None and definition.outcome_id not in selected:
            continue
        comparator_value = _per_100(
            values.get(("comparator", definition.event_name)),
            denominators.get("comparator"),
        )
        intervention_value = _per_100(
            values.get(("intervention", definition.event_name)),
            denominators.get("intervention"),
        )
        if definition.intervention_only and comparator_value is None:
            comparator_value = 0.0
        rows.append(
            {
                "outcomeId": definition.outcome_id,
                "eventName": definition.event_name,
                "label": definition.label,
                "color": definition.color,
                "interventionOnly": definition.intervention_only,
                "comparator": icon_grid_value(comparator_value),
                "intervention": icon_grid_value(intervention_value),
            }
        )
    return rows


def render_100_person_summary(rows: list[dict[str, Any]], *, title: str) -> None:
    if not rows:
        st.info("100-person visual summary is unavailable for these results.")
        return
    st.subheader(title)
    st.caption(ROUNDING_NOTE)
    legend = " ".join(
        f"<span style='display:inline-flex;align-items:center;margin-right:1rem;'>"
        f"<span style='width:0.8rem;height:0.8rem;background:{row['color']};"
        f"display:inline-block;margin-right:0.3rem;border-radius:2px;'></span>"
        f"{row['label']}</span>"
        for row in rows
    )
    st.markdown(legend, unsafe_allow_html=True)
    for row in rows:
        st.markdown(f"**{row['label']}**")
        left, right = st.columns(2)
        _render_arm_grid(left, "Without additional systematic screening", row["comparator"], row["color"])
        _render_arm_grid(right, "With targeted screening and preventive treatment", row["intervention"], row["color"])


def _render_arm_grid(container, label: str, value: dict[str, Any], color: str) -> None:
    container.markdown(label)
    container.markdown(_grid_html(value.get("filledIcons", 0), color), unsafe_allow_html=True)
    container.caption(value.get("displayPer100", "Unavailable"))


def _grid_html(filled_icons: int, color: str) -> str:
    filled = int(max(0, min(GRID_SIZE, filled_icons)))
    cells = []
    for idx in range(GRID_SIZE):
        active = idx < filled
        background = color if active else "#e5e7eb"
        border = color if active else "#d1d5db"
        cells.append(
            "<span style='width:0.72rem;height:0.72rem;"
            f"background:{background};border:1px solid {border};"
            "display:block;border-radius:2px;'></span>"
        )
    return (
        f"<div style='display:grid;grid-template-columns:repeat({GRID_COLUMNS},0.72rem);"
        "gap:0.16rem;margin:0.2rem 0 0.35rem 0;'>"
        + "".join(cells)
        + "</div>"
    )


def _records_to_frame(records: Any) -> pd.DataFrame:
    if isinstance(records, pd.DataFrame):
        return records.copy()
    if isinstance(records, list):
        return pd.DataFrame(records)
    return pd.DataFrame()


def _mean_values_by_arm_event(totals: pd.DataFrame) -> dict[tuple[str, str], float]:
    required = {"arm", "eventName", "value"}
    if not required.issubset(set(totals.columns)):
        return {}
    grouped = totals.groupby(["arm", "eventName"], dropna=False)["value"].mean()
    return {
        (str(arm), str(event)): float(value)
        for (arm, event), value in grouped.items()
        if _optional_float(value) is not None
    }


def _per_100(value: Any, denominator: Any) -> float | None:
    numerator = _optional_float(value)
    denom = _optional_float(denominator)
    if numerator is None or denom is None or denom <= 0:
        return None
    return numerator / denom * 100.0


def _optional_float(value: Any) -> float | None:
    if value in (None, "", []):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not isfinite(number):
        return None
    return number
