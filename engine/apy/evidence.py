from __future__ import annotations

from copy import deepcopy
import csv
from pathlib import Path
from typing import Any

from engine.apy.config import normalise_config
from engine.apy.costing import TARGET_CURRENCY, TARGET_PRICE_YEAR, normalise_cost_table
from engine.apy.economics import default_cost_items
from engine.apy.ltbi_state import is_clinician_ready_ltbi_state


EVIDENCE_REGISTRY_VERSION = "apy_assumption_evidence_registry_v1"
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REGISTRY_PATH = REPO_ROOT / "data" / "apy_assumption_evidence_registry.csv"
REVIEWED_NUMERICAL_STATUSES = {"configured_reviewed", "model_derived_reviewed"}
REVIEWED_EXCLUSION_STATUS = "reviewed_exclusion"
NON_READY_STATUSES = {
    "",
    "unresolved",
    "provisional",
    "unreviewed_repository_input",
    "legacy_placeholder",
    "migrated_legacy_unverified",
    "unresolved_development_compatibility",
}
ESSENTIAL_COST_IDS = {
    "cost.test_igra",
    "cost.test_tst",
    "cost.regimen_3hp",
    "cost.regimen_4r",
    "cost.regimen_3hr",
    "cost.regimen_6h",
    "cost.regimen_9h",
    "cost.active_tb_disease",
}
EXPECTED_COST_BASIS = {
    "cost.test_igra": {"per_person_screened"},
    "cost.test_tst": {"per_person_screened"},
    "cost.regimen_3hp": {"per_started_course"},
    "cost.regimen_4r": {"per_started_course"},
    "cost.regimen_3hr": {"per_started_course"},
    "cost.regimen_6h": {"per_started_course"},
    "cost.regimen_9h": {"per_started_course"},
    "cost.active_tb_disease": {"per_active_tb_case"},
    "cost.false_positive_incremental": {"per_false_positive_tpt_started"},
    "cost.program_setup": {"total_once_at_program_start"},
    "cost.program_running": {"annual_during_screening_window", "total_over_screening_programme"},
    "cost.tpt_adr_management": {"per_adr_stop"},
}
REQUIRED_ASSUMPTION_IDS = {
    "epi.baseline_recent_ltbi_proportion",
    "epi.recent_to_remote_transition_rate",
    "epi.ltbi_prevalence",
    "epi.active_tb_calibration_target_2y",
    "epi.igra_sensitivity",
    "epi.igra_specificity",
    "epi.tst_sensitivity",
    "epi.tst_specificity_bcg",
    "epi.tst_specificity_no_bcg",
    "epi.tpt_initiation",
    "epi.regimen_3hp_completion",
    "epi.regimen_4r_completion",
    "epi.regimen_3hr_completion",
    "epi.regimen_6h_completion",
    "epi.regimen_9h_completion",
    "epi.regimen_3hp_adr_stop",
    "epi.regimen_4r_adr_stop",
    "epi.regimen_3hr_adr_stop",
    "epi.regimen_6h_adr_stop",
    "epi.regimen_9h_adr_stop",
    "epi.regimen_full_efficacy",
    "epi.regimen_partial_efficacy",
    "epi.active_tb_progression",
    "cost.test_igra",
    "cost.test_tst",
    "cost.regimen_3hp",
    "cost.regimen_4r",
    "cost.regimen_3hr",
    "cost.regimen_6h",
    "cost.regimen_9h",
    "cost.active_tb_disease",
    "cost.false_positive_incremental",
    "cost.program_setup",
    "cost.program_running",
    "cost.tpt_adr_management",
    "daly.active_tb_disability_weight",
    "daly.active_tb_duration",
    "daly.tb_case_fatality_risk",
    "daly.yll_per_tb_death",
    "daly.tpt_health_loss",
    "daly.adr_health_loss",
    "daly.post_tb_sequelae",
    "threshold.gdp_per_capita",
}


def load_apy_evidence_registry(path: str | Path | None = None) -> list[dict[str, Any]]:
    source = Path(path) if path is not None else DEFAULT_REGISTRY_PATH
    with source.open("r", encoding="utf-8-sig", newline="") as f:
        return [_normalise_registry_row(row) for row in csv.DictReader(f)]


def assess_apy_reference_readiness(
    config: dict[str, Any] | None = None,
    economics_config: dict[str, Any] | None = None,
    registry: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    cfg = normalise_config(config or {})
    econ = deepcopy(economics_config or {})
    rows = deepcopy(registry) if registry is not None else load_apy_evidence_registry()
    by_id = {row["assumptionId"]: row for row in rows}
    unresolved: list[str] = []
    provisional: list[str] = []
    messages: list[str] = []
    missing_required = sorted(REQUIRED_ASSUMPTION_IDS - set(by_id))
    for assumption_id in missing_required:
        _mark_unresolved(assumption_id, unresolved)
        messages.append(f"Evidence registry is missing required assumption {assumption_id}.")

    category_ready = {
        "epidemiology": _category_ready(rows, "epidemiology", unresolved, provisional),
        "cost": _costs_ready(econ, by_id, unresolved, provisional, messages),
        "daly": _dalys_ready(rows, unresolved, provisional),
        "threshold": _threshold_ready(econ, by_id, unresolved, provisional),
    }
    if not is_clinician_ready_ltbi_state(cfg):
        _mark_unresolved("epi.baseline_recent_ltbi_proportion", unresolved)
        category_ready["epidemiology"] = False
        messages.append("APY baseline recent-LTBI proportion is not clinician-ready.")

    conflicts = _evidence_conflicts(rows)
    bundling_conflicts = _bundling_conflicts(rows)
    if conflicts:
        category_ready["epidemiology"] = False
    if bundling_conflicts:
        category_ready["cost"] = False

    if missing_required:
        for key in category_ready:
            category_ready[key] = False
    overall = all(category_ready.values()) and not conflicts and not bundling_conflicts
    return {
        "contractVersion": EVIDENCE_REGISTRY_VERSION,
        "epidemiologyReady": category_ready["epidemiology"],
        "costReady": category_ready["cost"],
        "dalyReady": category_ready["daly"],
        "thresholdReady": category_ready["threshold"],
        "overallClinicianReady": overall,
        "unresolvedAssumptionIds": sorted(set(unresolved)),
        "provisionalAssumptionIds": sorted(set(provisional)),
        "evidenceConflicts": conflicts,
        "bundlingConflicts": bundling_conflicts,
        "messages": messages,
        "readinessRows": readiness_matrix_rows(rows, category_ready),
        "registryRows": rows,
    }


def readiness_matrix_rows(
    registry: list[dict[str, Any]] | None = None,
    category_ready: dict[str, bool] | None = None,
) -> list[dict[str, Any]]:
    rows = deepcopy(registry) if registry is not None else load_apy_evidence_registry()
    category_ready = category_ready or {}
    out = []
    for row in rows:
        status = row.get("reviewStatus", "")
        ready = _row_reviewed(row)
        out.append(
            {
                "assumptionId": row.get("assumptionId", ""),
                "category": row.get("category", ""),
                "categoryReady": category_ready.get(row.get("category", ""), ready),
                "ready": ready,
                "reviewStatus": status,
                "inclusionStatus": row.get("inclusionStatus", ""),
                "sourceCitation": row.get("sourceCitation", ""),
                "unresolvedReason": row.get("unresolvedReason", ""),
                "provisional": row.get("provisional", True),
            }
        )
    return out


def validate_cost_item_evidence(
    cost_item: dict[str, Any],
    registry_row: dict[str, Any],
    *,
    target_currency: str = TARGET_CURRENCY,
    target_price_year: str = TARGET_PRICE_YEAR,
) -> list[str]:
    errors: list[str] = []
    assumption_id = registry_row.get("assumptionId", "")
    status = str(registry_row.get("reviewStatus") or "")
    if status not in REVIEWED_NUMERICAL_STATUSES:
        errors.append("Cost evidence requires reviewed numerical status.")
    if registry_row.get("provisional"):
        errors.append("Cost evidence is provisional.")
    if _blank(cost_item.get("originalCost")):
        errors.append("Original cost is missing.")
    if _blank(cost_item.get("originalCurrency")):
        errors.append("Original currency is missing.")
    if _blank(cost_item.get("targetCurrency")):
        errors.append("Target currency is missing.")
    if str(cost_item.get("originalCurrency")) != str(cost_item.get("targetCurrency")):
        errors.append("Foreign-currency costs require an explicit currency-conversion record.")
    if str(cost_item.get("targetCurrency")) != str(target_currency):
        errors.append("Cost target currency does not match analysis currency.")
    if str(cost_item.get("targetPriceYear")) != str(target_price_year):
        errors.append("Cost target price year does not match analysis price year.")
    if _blank(cost_item.get("originalPriceYear")):
        errors.append("Original source price year is missing.")
    if str(cost_item.get("sourceYearStatus")) not in {"explicit", "reviewed_inferred"}:
        errors.append("Source price year must be explicit or reviewed-inferred.")
    if not str(cost_item.get("sourceCitation") or "").strip():
        errors.append("Source citation is missing.")
    if cost_item.get("conversionStatus") != "valid":
        errors.append("Target-year cost conversion is not valid.")
    if _blank(cost_item.get("convertedTargetYearCost")):
        errors.append("Converted target-year cost is missing.")
    basis = str((cost_item.get("resourceUse") or {}).get("costBasis") or "")
    if assumption_id in EXPECTED_COST_BASIS and basis not in EXPECTED_COST_BASIS[assumption_id]:
        errors.append("Cost resource-use basis is incompatible with the event mapping.")
    return errors


def _category_ready(
    rows: list[dict[str, Any]],
    category: str,
    unresolved: list[str],
    provisional: list[str],
) -> bool:
    ready = True
    for row in rows:
        if row.get("category") != category or row.get("inclusionStatus") not in {"included", "excluded", "bundled"}:
            continue
        if not _row_reviewed(row):
            ready = False
            _mark_unresolved(row["assumptionId"], unresolved)
        if row.get("provisional"):
            _mark_unresolved(row["assumptionId"], provisional)
    return ready


def _costs_ready(
    econ: dict[str, Any],
    by_id: dict[str, dict[str, Any]],
    unresolved: list[str],
    provisional: list[str],
    messages: list[str],
) -> bool:
    items = normalise_cost_table(econ.get("costItems") or default_cost_items())
    item_by_assumption = {f"cost.{item.get('costItemId')}": item for item in items}
    ready = True
    for assumption_id, row in by_id.items():
        if row.get("category") != "cost":
            continue
        if row.get("inclusionStatus") == "bundled":
            continue
        if row.get("reviewStatus") == REVIEWED_EXCLUSION_STATUS:
            if assumption_id in ESSENTIAL_COST_IDS and row.get("inclusionStatus") == "excluded":
                ready = False
                _mark_unresolved(assumption_id, unresolved)
                messages.append(f"Essential cost {assumption_id} cannot be generically excluded.")
            continue
        item = item_by_assumption.get(assumption_id, {})
        errors = validate_cost_item_evidence(item, row)
        if errors:
            ready = False
            _mark_unresolved(assumption_id, unresolved)
            messages.extend(f"{assumption_id}: {error}" for error in errors)
        if row.get("provisional"):
            _mark_unresolved(assumption_id, provisional)
    return ready


def _dalys_ready(rows: list[dict[str, Any]], unresolved: list[str], provisional: list[str]) -> bool:
    ready = True
    for row in rows:
        if row.get("category") != "daly":
            continue
        if row.get("inclusionStatus") == "excluded":
            row_ready = row.get("reviewStatus") == REVIEWED_EXCLUSION_STATUS and bool(str(row.get("notes") or "").strip())
        else:
            row_ready = row.get("reviewStatus") in REVIEWED_NUMERICAL_STATUSES
        if not row_ready:
            ready = False
            _mark_unresolved(row["assumptionId"], unresolved)
        if row.get("provisional"):
            _mark_unresolved(row["assumptionId"], provisional)
    return ready


def _threshold_ready(
    econ: dict[str, Any],
    by_id: dict[str, dict[str, Any]],
    unresolved: list[str],
    provisional: list[str],
) -> bool:
    row = by_id.get("threshold.gdp_per_capita", {})
    threshold = econ.get("threshold") or {}
    ready = (
        _row_reviewed(row)
        and not _blank(threshold.get("value"))
        and float(threshold.get("value")) > 0
        and str(threshold.get("currency") or "") == TARGET_CURRENCY
        and str(threshold.get("referenceYear") or "") == TARGET_PRICE_YEAR
        and bool(str(threshold.get("source") or "").strip())
    )
    if not ready:
        _mark_unresolved("threshold.gdp_per_capita", unresolved)
    if row.get("provisional"):
        _mark_unresolved("threshold.gdp_per_capita", provisional)
    return ready


def _evidence_conflicts(rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    conflicts = []
    for row in rows:
        text = " ".join(str(row.get(field) or "") for field in ("notes", "unresolvedReason"))
        if "conflict" in text.lower():
            conflicts.append(
                {
                    "assumptionId": row.get("assumptionId", ""),
                    "message": text,
                }
            )
    return conflicts


def _bundling_conflicts(rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    by_id = {row["assumptionId"]: row for row in rows}
    conflicts: list[dict[str, str]] = []
    for row in rows:
        assumption_id = row.get("assumptionId", "")
        dest = row.get("bundledIntoAssumptionId", "")
        if not dest:
            continue
        if dest not in by_id:
            conflicts.append({"assumptionId": assumption_id, "message": "Bundled destination does not exist."})
            continue
        if not _row_reviewed(by_id[dest]):
            conflicts.append({"assumptionId": assumption_id, "message": "Bundled destination is not reviewed and valid."})
        if dest == assumption_id or _has_bundle_cycle(assumption_id, by_id):
            conflicts.append({"assumptionId": assumption_id, "message": "Circular bundling relationship."})
        if row.get("inclusionStatus") != "bundled":
            conflicts.append({"assumptionId": assumption_id, "message": "Bundled item must have inclusionStatus=bundled."})
        if row.get("reviewStatus") not in {REVIEWED_EXCLUSION_STATUS, *REVIEWED_NUMERICAL_STATUSES}:
            conflicts.append({"assumptionId": assumption_id, "message": "Bundling requires reviewed status."})
        if str(row.get("currentValue") or "").strip():
            conflicts.append({"assumptionId": assumption_id, "message": "Component is both separately costed and bundled."})
    return conflicts


def _has_bundle_cycle(start: str, by_id: dict[str, dict[str, Any]]) -> bool:
    seen = set()
    current = start
    while current:
        if current in seen:
            return True
        seen.add(current)
        current = by_id.get(current, {}).get("bundledIntoAssumptionId", "")
    return False


def _row_reviewed(row: dict[str, Any]) -> bool:
    status = str(row.get("reviewStatus") or "")
    if row.get("provisional") or status in NON_READY_STATUSES:
        return False
    if row.get("inclusionStatus") == "excluded":
        return status == REVIEWED_EXCLUSION_STATUS and bool(str(row.get("notes") or "").strip())
    return status in REVIEWED_NUMERICAL_STATUSES and bool(str(row.get("sourceCitation") or "").strip())


def _normalise_registry_row(row: dict[str, Any]) -> dict[str, Any]:
    out = dict(row)
    out["provisional"] = str(out.get("provisional", "")).strip().lower() == "true"
    return out


def _blank(value: Any) -> bool:
    return value in (None, "", [])


def _mark_unresolved(assumption_id: str, target: list[str]) -> None:
    if assumption_id:
        target.append(assumption_id)
