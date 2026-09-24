"""Health-economic helpers for the general application.

Economic recalculation reuses the completed epidemiological run; editing costs never
reruns the epidemiological model. The formulas are those of the existing economics
engine (engine/apy/event_ledger_economics.py).
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from typing import Any

from app.general.terminology import display_text


EDITABLE_COST_ITEMS = (
    ("test_igra", "IGRA screening test, per person tested"),
    ("test_tst", "TST screening test, per person tested"),
    ("regimen_3hp", "3HP preventive treatment, per course started"),
    ("regimen_4r", "4R preventive treatment, per course started"),
    ("regimen_3hr", "3HR preventive treatment, per course started"),
    ("regimen_6h", "6H preventive treatment, per course started"),
    ("regimen_9h", "9H preventive treatment, per course started"),
    ("tpt_adr_management", "Adverse-event management, per treatment stopped for adverse events"),
    ("return_for_results", "Return visit for results, per person screened"),
    ("clinical_review", "Clinical review, per preventive treatment started"),
    ("active_tb_exclusion_workup", "Active-TB exclusion work-up, per preventive treatment started"),
    ("active_tb_disease", "Active TB disease management, per case"),
    ("program_setup", "Programme setup, total"),
    ("program_running", "Programme running, total"),
)
DEMONSTRATION_COST_SOURCE = "Demonstration local-pathway working assumption"
USER_COST_SOURCE = "User-defined"


def default_economics_config() -> dict[str, Any]:
    from engine.apy.working_defaults import build_unified_working_default_preset

    return deepcopy(build_unified_working_default_preset()["economicsConfig"])


def economics_config_hash(config: dict[str, Any]) -> str:
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def source_label(citation: Any) -> str:
    text = str(citation or "")
    if text == USER_COST_SOURCE:
        return USER_COST_SOURCE
    if text.startswith("Dale KD"):
        return "Dale et al. 2022, Am J Epidemiol (Australian inputs)"
    return display_text(text, fallback=DEMONSTRATION_COST_SOURCE) or "Not recorded"


def cost_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    from engine.apy.costing import normalise_cost_table

    items = {item.get("costItemId"): item for item in normalise_cost_table(config.get("costItems") or [])}
    rows = []
    for item_id, label in EDITABLE_COST_ITEMS:
        item = items.get(item_id)
        if item is None:
            continue
        rows.append(
            {
                "id": item_id,
                "Item": label,
                "Unit cost": item.get("convertedTargetYearCost"),
                "Currency and year": f"{item.get('targetCurrency', '')} {item.get('targetPriceYear', '')}".strip(),
                "Source": source_label(item.get("sourceCitation")),
            }
        )
    return rows


def apply_cost_edits(config: dict[str, Any], edits: dict[str, float]) -> dict[str, Any]:
    """Return a copy with user-defined unit costs (target-year prices, no re-inflation)."""
    out = deepcopy(config)
    for item in out.get("costItems") or []:
        item_id = item.get("costItemId")
        if item_id not in edits:
            continue
        item["originalCost"] = float(edits[item_id])
        item["originalPriceYear"] = item.get("targetPriceYear")
        item["originalCurrency"] = item.get("targetCurrency", item.get("originalCurrency"))
        item["sourceCitation"] = USER_COST_SOURCE
        item["costRecordType"] = "source"
        item["notes"] = "User-defined value entered in the general application (target-year prices)."
    return out


def changed_cost_edits(current_rows: list[dict[str, Any]], edited_rows: list[dict[str, Any]]) -> dict[str, float]:
    before = {row["id"]: row["Unit cost"] for row in current_rows}
    edits = {}
    for row in edited_rows:
        value = row.get("Unit cost")
        if value is None or value != value:
            continue
        previous = before.get(row.get("id"))
        if previous is None or abs(float(value) - float(previous)) > 1e-9:
            edits[row["id"]] = float(value)
    return edits


def run_economics(bundle: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    from engine.apy.event_ledger_economics import run_event_ledger_health_economics

    return run_event_ledger_health_economics(bundle, config)


def icer_cloud_points(results: dict[str, Any], profile: str = "primary") -> list[dict[str, float]]:
    """Paired (DALYs averted, incremental cost) per simulated population; comparator at the origin."""
    frame = results.get("replicateResults")
    if frame is None or not hasattr(frame, "to_dict"):
        return []
    rows = frame[(frame["discountProfile"] == profile) & (frame["economicPairComplete"] == True)]  # noqa: E712
    return [
        {
            "replicate": int(row["pairedReplicateId"]),
            "dalysAverted": float(row["dalysAverted"]),
            "incrementalCost": float(row["incrementalCost"]),
        }
        for row in rows.to_dict(orient="records")
    ]
