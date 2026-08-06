from __future__ import annotations

from copy import deepcopy
import math
from typing import Any

import numpy as np
import pandas as pd

from engine.apy.costing import (
    TARGET_CURRENCY,
    TARGET_PRICE_YEAR,
    normalise_cost_table,
    valid_converted_cost,
)
from engine.apy.event_ledger import EVENT_LEDGER_CONTRACT_VERSION
from engine.apy.scenario import DIRECT_EFFECTS_SCOPE_STATEMENT


HEALTH_ECONOMICS_CONTRACT_VERSION = "ltbi_health_economics_results_v2"
PRIMARY_DISCOUNT_RATE = 0.03
COMPARISON_DISCOUNT_RATE = 0.0
NUMERIC_REVIEWED_STATUSES = {"configured_reviewed", "model_derived_reviewed"}
EXCLUSION_REVIEWED_STATUS = "reviewed_exclusion"
INTERPRETABLE_ICER_CLASSIFICATION = "increased cost with health gain"

EXPECTED_COST_BASES = {
    "test_igra": {"per_person_screened"},
    "test_tst": {"per_person_screened"},
    "regimen_3hp": {"per_started_course"},
    "regimen_4r": {"per_started_course"},
    "regimen_3hr": {"per_started_course"},
    "regimen_6h": {"per_started_course"},
    "regimen_9h": {"per_started_course"},
    "active_tb_disease": {"per_active_tb_case"},
    "false_positive_incremental": {"per_false_positive_tpt_started"},
    "program_setup": {"total_once_at_program_start"},
    "program_running": {"annual_during_screening_window", "total_over_screening_programme"},
    "tpt_adr_management": {"per_adr_stop"},
}

COMPONENT_SPECS = {
    "screeningTestCost": ("test", "screened"),
    "tptRegimenCost": ("regimen", "tpt_started_total"),
    "falsePositiveIncrementalCost": ("false_positive", "tpt_started_false_positive"),
    "activeTBDiseaseCost": ("active_tb", "active_tb_cases"),
    "programSetupCost": ("setup", "_setup_quantity"),
    "programRunningCost": ("running", "_running_quantity"),
    "adrManagementCost": ("adr", "tpt_adr_stop_total"),
}

PROGRAM_COMPONENTS = [
    "screeningTestCost",
    "tptRegimenCost",
    "falsePositiveIncrementalCost",
    "programSetupCost",
    "programRunningCost",
    "adrManagementCost",
]
COST_COMPONENTS = PROGRAM_COMPONENTS + ["activeTBDiseaseCost"]
DALY_COMPONENTS = ["activeTBYLD", "activeTBYLL", "tptHealthLossDALYs", "adrHealthLossDALYs"]
ID_COLS = ["discountRate", "modelType", "valueType", "replicateId", "pairedReplicateId", "replicateSeed"]


def default_daly_assumptions() -> dict[str, Any]:
    return {
        "activeTBDisabilityWeight": _assumption("Active-TB disability weight"),
        "activeTBDurationYears": _assumption("Active-TB disability duration"),
        "tbCaseFatalityRisk": _assumption("TB case-fatality risk"),
        "yllPerTBDeath": _assumption("Years of life lost per TB death"),
        "includeTPTHealthLoss": False,
        "tptHealthLossExclusionStatus": "unresolved",
        "tptHealthLossExclusionRationale": "",
        "dalyLossPerTPTStarted": _assumption("DALY loss per TPT course started"),
        "includeADRHealthLoss": False,
        "adrHealthLossExclusionStatus": "unresolved",
        "adrHealthLossExclusionRationale": "",
        "dalyLossPerADRStop": _assumption("DALY loss per ADR-related stop"),
        "includePostTBSequelae": False,
        "postTBSequelaeStatus": "unresolved",
        "postTBSequelaeExclusionRationale": "",
        "postTBSequelaeNotes": "Excluded from primary acute analysis unless reviewed assumptions are supplied.",
        "method": "scalar_expected_daly_per_active_tb_case",
        "notes": "Scalar active-TB DALY method; not an age-specific life-table calculation.",
    }


def run_event_ledger_health_economics(
    results_or_bundle: dict[str, Any],
    econ_config: dict[str, Any] | None = None,
    *,
    discount_rates: list[float] | None = None,
) -> dict[str, Any]:
    econ_config = deepcopy(econ_config or {})
    ledger = _extract_event_ledger(results_or_bundle)
    validation = ledger.get("validation") or {}
    if not ledger or validation.get("isValid") is not True:
        raise ValueError("Authoritative health economics requires a valid APY event ledger.")

    discount_rates = discount_rates or [PRIMARY_DISCOUNT_RATE, COMPARISON_DISCOUNT_RATE]
    ledger_metadata = deepcopy(ledger.get("metadata") or {})
    assumptions = _build_assumptions(econ_config, ledger_metadata)
    cost_items = normalise_cost_table(econ_config.get("costItems") or [])
    unresolved: list[dict[str, Any]] = []
    warnings: list[str] = []

    cost_lookup = _cost_lookup(cost_items, econ_config, ledger_metadata, assumptions, unresolved)
    before_daly = len(unresolved)
    daly_inputs = _resolve_daly_inputs(assumptions["daly"], unresolved)
    daly_inputs["globalIncompleteReasons"] = [
        item["field"] for item in unresolved[before_daly:] if item["field"].startswith("dalyAssumptions")
    ]
    threshold = _resolve_threshold(econ_config.get("threshold") or {}, assumptions, unresolved)
    annual_events = _complete_calendar(_events_wide(ledger["annualEvents"]), ledger_metadata, assumptions)

    annual_rows = [
        _annual_economic_row(row, ledger_metadata, cost_lookup, daly_inputs, assumptions, float(rate))
        for _, row in annual_events.iterrows()
        for rate in discount_rates
    ]
    annual = pd.DataFrame(annual_rows)
    _allocate_program_running_total(annual, cost_lookup, ledger_metadata)
    _recompute_completeness_and_totals(annual)
    replicate_results = _replicate_results(annual, threshold)
    summaries = _summaries(replicate_results, ledger_metadata, threshold)
    validation_report = _validate_economic_result(annual, replicate_results, summaries, threshold)

    provisional = bool(ledger_metadata.get("ltbiStateProvisional") or ledger_metadata.get("ltbiStateWarning"))
    if provisional:
        warnings.append("Epidemiological inputs are provisional; no clinician-ready cost-effectiveness conclusion is produced.")
    if unresolved:
        warnings.append("One or more economic inputs are unresolved; unavailable components are not replaced with zero.")
    if not bool(validation_report.get("economicallyComplete")):
        warnings.append("Economic result is structurally valid but incomplete for at least one required pair or component.")

    economically_complete = bool(validation_report["economicallyComplete"])
    threshold_complete = threshold is not None
    unresolved_daly_inputs = any(item["field"].startswith("dalyAssumptions") for item in unresolved)
    conclusion_permitted = bool(economically_complete and threshold_complete and not provisional and not unresolved_daly_inputs)
    validation_report["conclusionPermitted"] = conclusion_permitted
    result = {
        "available": True,
        "source": "event_ledger_health_economics_v2",
        "contractVersion": HEALTH_ECONOMICS_CONTRACT_VERSION,
        "metadata": {
            "economicContractVersion": HEALTH_ECONOMICS_CONTRACT_VERSION,
            "eventLedgerContractVersion": ledger_metadata.get("contractVersion"),
            "scenarioId": ledger_metadata.get("scenarioId"),
            "modelType": ledger_metadata.get("modelType"),
            "valueType": _first_value(annual, "valueType"),
            "perspective": assumptions["metadata"].get("perspective"),
            "targetCurrency": assumptions["metadata"].get("targetCurrency", TARGET_CURRENCY),
            "targetPriceYear": assumptions["metadata"].get("targetPriceYear", TARGET_PRICE_YEAR),
            "economicHorizonYears": assumptions["metadata"].get("economicHorizonYears"),
            "scopeStatement": DIRECT_EFFECTS_SCOPE_STATEMENT,
            "ltbiStateProvisional": ledger_metadata.get("ltbiStateProvisional"),
            "ltbiStateDevelopmentCompatibilityMode": (
                ledger_metadata.get("ltbiStateAssumptionStatus")
                == "unresolved_development_compatibility"
            ),
            "primaryDiscountRate": PRIMARY_DISCOUNT_RATE,
            "comparisonDiscountRate": COMPARISON_DISCOUNT_RATE,
            "isProvisional": bool(provisional or unresolved_daly_inputs or not economically_complete),
            "conclusionPermitted": conclusion_permitted,
        },
        "assumptions": assumptions,
        "costItems": cost_items,
        "annualByArm": annual,
        "replicateResults": replicate_results,
        "summaries": summaries,
        "validation": validation_report,
        "unresolvedInputs": unresolved,
        "warnings": warnings,
        "provenance": {
            "epidemiologySource": "APY event ledger",
            "costSource": "authoritative costItems with valid convertedTargetYearCost",
            "dalySource": "explicit source-aware DALY assumptions",
            "discounting": "modelYear discount factor = 1 / (1 + rate)^modelYear",
        },
    }
    _attach_legacy_compatibility_fields(result, econ_config, cost_lookup)
    result["summaryRows"] = _summary_rows(result)
    result["summaryTable"] = result["summaryRows"]
    result["status"] = {
        "isComplete": economically_complete,
        "missingInputs": [item["field"] for item in unresolved],
        "notCalculated": sorted(set(_not_calculated(result) + result["_legacyCompatibilityStatus"]["notCalculated"])),
        "messages": warnings + result["_legacyCompatibilityStatus"]["messages"],
        "partialCalculations": result["_legacyCompatibilityStatus"].get("partialCalculations", []),
        "validationReport": result["validation"],
    }
    return result


def _assumption(label: str) -> dict[str, Any]:
    return {"label": label, "value": None, "source": "", "status": "unresolved", "notes": "", "provisional": True}


def _build_assumptions(econ_config: dict[str, Any], metadata: dict[str, Any]) -> dict[str, Any]:
    assumptions = {
        "metadata": deepcopy(econ_config.get("metadata") or {}),
        "discounting": deepcopy(econ_config.get("discounting") or {}),
        "threshold": deepcopy(econ_config.get("threshold") or {}),
        "daly": default_daly_assumptions(),
        "optionalCostExclusions": deepcopy(econ_config.get("optionalCostExclusions") or {}),
    }
    if isinstance(econ_config.get("dalyAssumptions"), dict):
        _deep_update(assumptions["daly"], econ_config["dalyAssumptions"])
    follow_up = metadata.get("followUpHorizonYears") or metadata.get("followHorizon") or metadata.get("followUpHorizon")
    assumptions["metadata"].setdefault("targetCurrency", TARGET_CURRENCY)
    assumptions["metadata"].setdefault("targetPriceYear", TARGET_PRICE_YEAR)
    assumptions["metadata"].setdefault("perspective", "Australian health-system perspective")
    assumptions["metadata"].setdefault("economicHorizonYears", follow_up)
    return assumptions


def _cost_lookup(
    cost_items: list[dict[str, Any]],
    econ_config: dict[str, Any],
    metadata: dict[str, Any],
    assumptions: dict[str, Any],
    unresolved: list[dict[str, Any]],
) -> dict[str, Any]:
    test_type = str(_config_value(econ_config, metadata, "testType") or "IGRA").upper()
    regimen = str(_config_value(econ_config, metadata, "regimen") or "3HP").upper()
    ids = {
        "test": {"IGRA": "test_igra", "TST": "test_tst"}.get(test_type),
        "regimen": {
            "3HP": "regimen_3hp",
            "4R": "regimen_4r",
            "3HR": "regimen_3hr",
            "6H": "regimen_6h",
            "9H": "regimen_9h",
        }.get(regimen),
        "active_tb": "active_tb_disease",
        "false_positive": "false_positive_incremental",
        "setup": "program_setup",
        "running": "program_running",
        "adr": "tpt_adr_management",
    }
    by_id = {item.get("costItemId"): item for item in cost_items}
    out: dict[str, Any] = {"ids": ids, "items": by_id, "basisIssues": {}, "exclusions": {}}
    target_currency = assumptions["metadata"].get("targetCurrency", TARGET_CURRENCY)
    target_year = assumptions["metadata"].get("targetPriceYear", TARGET_PRICE_YEAR)
    for name, item_id in ids.items():
        if item_id is None:
            unresolved.append(_unresolved(f"costItems.{name}", "No cost item mapping is available."))
            out[name] = None
            continue
        item = by_id.get(item_id, {})
        basis = _basis(item)
        expected = EXPECTED_COST_BASES.get(item_id, set())
        if expected and basis not in expected:
            out["basisIssues"][name] = f"{item_id} cost basis '{basis}' is incompatible with event mapping."
            unresolved.append(
                _unresolved(
                    f"costItems.{item_id}.resourceUse.costBasis",
                    out["basisIssues"][name],
                )
            )
            out[name] = None
            continue
        if item_id == "test_tst" and "returnVisitForReading" not in (item.get("resourceUse") or {}):
            out["basisIssues"][name] = "TST cost must indicate whether return-reading costs are bundled."
            unresolved.append(
                _unresolved(
                    "costItems.test_tst.resourceUse.returnVisitForReading",
                    out["basisIssues"][name],
                )
            )
            out[name] = None
            continue
        if item and (
            item.get("targetCurrency") != target_currency
            or item.get("targetPriceYear") != target_year
        ):
            unresolved.append(
                _unresolved(
                    f"costItems.{item_id}.targetPriceYear",
                    f"{item_id} target currency/year does not match the analysis target.",
                )
            )
            out[name] = None
            continue
        value = valid_converted_cost(item)
        out[name] = value
        if value is None:
            exclusion = _reviewed_cost_exclusion(assumptions, item_id)
            if exclusion:
                out["exclusions"][name] = exclusion
            elif name != "adr":
                unresolved.append(_unresolved(f"costItems.{item_id}", f"{item_id} lacks valid converted target-year cost."))
    running_item = by_id.get("program_running", {})
    running_basis = _basis(running_item)
    if out.get("running") is not None and running_basis not in EXPECTED_COST_BASES["program_running"]:
        unresolved.append(_unresolved("costItems.program_running.resourceUse.costBasis", "Programme running cost requires an explicit allocation basis."))
        out["running"] = None
    if (
        out.get("running") is not None
        and running_basis == "total_over_screening_programme"
        and (running_item.get("resourceUse") or {}).get("allocationMethod") != "proportional_to_screening_volume"
    ):
        unresolved.append(_unresolved("costItems.program_running.resourceUse.allocationMethod", "Total-over-programme running cost requires an explicit allocation method."))
        out["running"] = None
    out["runningBasis"] = running_basis
    out["runningAllocation"] = (running_item.get("resourceUse") or {}).get("allocationMethod")
    return out


def _resolve_daly_inputs(daly: dict[str, Any], unresolved: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    specs = {
        "activeTBDisabilityWeight": (0.0, 1.0),
        "activeTBDurationYears": (0.0, None),
        "tbCaseFatalityRisk": (0.0, 1.0),
        "yllPerTBDeath": (0.0, None),
    }
    for field, bounds in specs.items():
        out[field] = _reviewed_numeric(daly, field, bounds, unresolved)

    include_tpt = bool(daly.get("includeTPTHealthLoss"))
    out["includeTPTHealthLoss"] = include_tpt
    if include_tpt:
        out["dalyLossPerTPTStarted"] = _reviewed_numeric(daly, "dalyLossPerTPTStarted", (0.0, None), unresolved)
        rec = daly.get("dalyLossPerTPTStarted") or {}
        if not str(rec.get("notes") or "").strip():
            unresolved.append(_unresolved("dalyAssumptions.dalyLossPerTPTStarted.notes", "Included TPT health loss requires notes."))
    else:
        out["dalyLossPerTPTStarted"] = None
        _reviewed_exclusion(
            daly,
            "tptHealthLossExclusionStatus",
            "tptHealthLossExclusionRationale",
            "dalyAssumptions.includeTPTHealthLoss",
            unresolved,
        )

    include_adr = bool(daly.get("includeADRHealthLoss"))
    out["includeADRHealthLoss"] = include_adr
    if include_adr:
        out["dalyLossPerADRStop"] = _reviewed_numeric(daly, "dalyLossPerADRStop", (0.0, None), unresolved)
        rec = daly.get("dalyLossPerADRStop") or {}
        if not str(rec.get("notes") or "").strip():
            unresolved.append(_unresolved("dalyAssumptions.dalyLossPerADRStop.notes", "Included ADR health loss requires notes."))
    else:
        out["dalyLossPerADRStop"] = None
        _reviewed_exclusion(
            daly,
            "adrHealthLossExclusionStatus",
            "adrHealthLossExclusionRationale",
            "dalyAssumptions.includeADRHealthLoss",
            unresolved,
        )

    if bool(daly.get("includePostTBSequelae")):
        unresolved.append(_unresolved("dalyAssumptions.includePostTBSequelae", "Post-TB sequelae are outside the Milestone 3 acute primary analysis."))
    else:
        _reviewed_exclusion(
            daly,
            "postTBSequelaeStatus",
            "postTBSequelaeExclusionRationale",
            "dalyAssumptions.includePostTBSequelae",
            unresolved,
        )
    return out


def _resolve_threshold(threshold: dict[str, Any], assumptions: dict[str, Any], unresolved: list[dict[str, Any]]) -> float | None:
    value = _num(threshold.get("value"))
    status = str(threshold.get("status") or "")
    target_currency = assumptions["metadata"].get("targetCurrency", TARGET_CURRENCY)
    target_year = assumptions["metadata"].get("targetPriceYear", TARGET_PRICE_YEAR)
    if value is None or value <= 0:
        unresolved.append(_unresolved("threshold.value", "Willingness-to-pay threshold must be > 0; NMB unavailable."))
        return None
    if (
        not threshold.get("currency")
        or threshold.get("referenceYear") in (None, "", [])
        or not threshold.get("source")
        or status not in NUMERIC_REVIEWED_STATUSES
    ):
        unresolved.append(_unresolved("threshold.provenance", "Threshold requires currency, year, source and reviewed numeric status."))
        return None
    if threshold.get("currency") != target_currency or str(threshold.get("referenceYear")) != str(target_year):
        unresolved.append(_unresolved("threshold.alignment", "Threshold currency/year must match the economic target currency/year."))
        return None
    return value


def _events_wide(annual_events: Any) -> pd.DataFrame:
    frame = annual_events if isinstance(annual_events, pd.DataFrame) else pd.DataFrame(annual_events)
    id_cols = [c for c in frame.columns if c not in {"eventName", "value"}]
    work = frame.copy()
    sentinel = "__MISSING__"
    for col in id_cols:
        work[col] = work[col].where(work[col].notna(), sentinel)
    wide = work.pivot_table(index=id_cols, columns="eventName", values="value", aggfunc="sum").reset_index()
    for col in id_cols:
        wide[col] = wide[col].mask(wide[col] == sentinel, np.nan)
    event_cols = [col for col in wide.columns if col not in id_cols]
    wide[event_cols] = wide[event_cols].fillna(0.0)
    return wide


def _complete_calendar(events: pd.DataFrame, metadata: dict[str, Any], assumptions: dict[str, Any]) -> pd.DataFrame:
    economic_horizon = _num(assumptions["metadata"].get("economicHorizonYears"))
    if economic_horizon is None:
        economic_horizon = _num(metadata.get("followUpHorizonYears") or metadata.get("followHorizon")) or 0
    horizon_years = list(range(max(0, int(math.ceil(economic_horizon)))))
    event_years = sorted({int(y) for y in pd.to_numeric(events.get("modelYear"), errors="coerce").dropna()})
    years = sorted(set(horizon_years + event_years + [0]))
    id_cols = [c for c in events.columns if c not in _event_columns(events) and c not in {"arm", "modelYear", "timeInterval", "withinFollowUp"}]
    pairs = events[id_cols].drop_duplicates()
    arms = ["comparator", "intervention"]
    rows = []
    for _, pair in pairs.iterrows():
        for arm in arms:
            for year in years:
                row = pair.to_dict()
                row.update(
                    {
                        "arm": arm,
                        "modelYear": year,
                        "timeInterval": f"[{year}, {year + 1})",
                        "withinFollowUp": year < economic_horizon,
                    }
                )
                rows.append(row)
    grid = pd.DataFrame(rows)
    merge_cols = id_cols + ["arm", "modelYear"]
    merged = grid.merge(events, how="left", on=merge_cols, suffixes=("", "_event"))
    for col in ["timeInterval", "withinFollowUp"]:
        if f"{col}_event" in merged:
            merged[col] = merged[f"{col}_event"].combine_first(merged[col])
            merged = merged.drop(columns=[f"{col}_event"])
    for col in _event_columns(merged):
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0.0)
    return merged


def _annual_economic_row(row: pd.Series, metadata: dict[str, Any], costs: dict[str, Any], daly: dict[str, Any], assumptions: dict[str, Any], rate: float) -> dict[str, Any]:
    year = int(row.get("modelYear", 0))
    factor = 1.0 / ((1.0 + rate) ** year)
    arm = row.get("arm")
    economic_horizon = _num(assumptions["metadata"].get("economicHorizonYears")) or 0.0
    within_economic_horizon = year < economic_horizon and bool(row.get("withinFollowUp", True))
    quantities = {
        "screened": _event(row, "screened"),
        "tpt_started_total": _event(row, "tpt_started_total"),
        "tpt_started_false_positive": _event(row, "tpt_started_false_positive"),
        "tpt_adr_stop_total": _event(row, "tpt_adr_stop_total"),
        "active_tb_cases": _event(row, "active_tb_cases"),
        "active_tb_cases_prevented": _event(row, "active_tb_cases_prevented"),
        "infection_effectively_treated_total": _event(row, "infection_effectively_treated_total"),
    }
    setup_qty = 1.0 if arm == "intervention" and year == 0 else 0.0
    running_qty = _running_quantity(row, costs, metadata) if arm == "intervention" else 0.0
    component_values = {
        "screeningTestCost": _component_cost(quantities["screened"], costs, "test"),
        "tptRegimenCost": _component_cost(quantities["tpt_started_total"], costs, "regimen"),
        "falsePositiveIncrementalCost": _component_cost(quantities["tpt_started_false_positive"], costs, "false_positive"),
        "activeTBDiseaseCost": _component_cost(quantities["active_tb_cases"], costs, "active_tb"),
        "programSetupCost": _component_cost(setup_qty, costs, "setup"),
        "programRunningCost": _component_cost(running_qty, costs, "running"),
        "adrManagementCost": _component_cost(quantities["tpt_adr_stop_total"], costs, "adr"),
    }
    missing_cost = _missing_cost_components(component_values, quantities, setup_qty, running_qty, costs)
    yld = _mul3(quantities["active_tb_cases"], daly.get("activeTBDisabilityWeight"), daly.get("activeTBDurationYears"))
    deaths = _mul(quantities["active_tb_cases"], daly.get("tbCaseFatalityRisk"))
    yll = _mul(deaths, daly.get("yllPerTBDeath"))
    tpt_loss = (
        _mul(quantities["tpt_started_total"], daly.get("dalyLossPerTPTStarted"))
        if daly.get("includeTPTHealthLoss")
        else 0.0
    )
    adr_loss = (
        _mul(quantities["tpt_adr_stop_total"], daly.get("dalyLossPerADRStop"))
        if daly.get("includeADRHealthLoss")
        else 0.0
    )
    daly_values = {
        "activeTBYLD": yld,
        "activeTBYLL": yll,
        "tptHealthLossDALYs": tpt_loss,
        "adrHealthLossDALYs": adr_loss,
    }
    missing_daly = [name for name, value in daly_values.items() if value is None]
    missing_daly.extend(daly.get("globalIncompleteReasons") or [])
    cost_complete = not missing_cost
    daly_complete = not missing_daly
    included = bool(within_economic_horizon)
    exclusion_reason = "" if included else "outside economic horizon or follow-up"
    total_cost = None if not cost_complete else float(sum(component_values.values()))
    total_dalys = None if not daly_complete else float(sum(daly_values.values()))
    return {
        "scenarioId": metadata.get("scenarioId"),
        "modelType": row.get("modelType"),
        "valueType": row.get("valueType"),
        "replicateId": row.get("replicateId"),
        "pairedReplicateId": row.get("pairedReplicateId"),
        "replicateSeed": row.get("replicateSeed"),
        "arm": arm,
        "modelYear": year,
        "timeInterval": row.get("timeInterval"),
        "withinFollowUp": bool(row.get("withinFollowUp", True)),
        "withinEconomicHorizon": bool(within_economic_horizon),
        "discountRate": rate,
        "discountFactor": factor,
        **quantities,
        **component_values,
        **daly_values,
        "costComplete": bool(cost_complete),
        "dalyComplete": bool(daly_complete),
        "missingCostComponents": ";".join(missing_cost),
        "missingDALYComponents": ";".join(missing_daly),
        "includedInEconomicAnalysis": included,
        "economicExclusionReason": exclusion_reason,
        "totalUndiscountedCost": total_cost,
        "totalDiscountedCost": None if total_cost is None else total_cost * factor,
        "totalUndiscountedDALYs": total_dalys,
        "totalDiscountedDALYs": None if total_dalys is None else total_dalys * factor,
    }


def _running_quantity(row: pd.Series, costs: dict[str, Any], metadata: dict[str, Any]) -> float:
    basis = costs.get("runningBasis")
    if basis == "annual_during_screening_window":
        window = float(metadata.get("screeningWindowYears") or metadata.get("screeningWindow") or 0)
        return 1.0 if float(row.get("modelYear", 0)) < window else 0.0
    if basis == "total_over_screening_programme" and costs.get("runningAllocation") == "proportional_to_screening_volume":
        return 0.0
    return 0.0


def _component_cost(quantity: float, costs: dict[str, Any], key: str) -> float | None:
    value = _num(costs.get(key))
    if abs(float(quantity)) <= 1e-12:
        return 0.0
    if value is None and key in costs.get("exclusions", {}):
        return 0.0
    return None if value is None else float(quantity) * value


def _missing_cost_components(
    values: dict[str, float | None],
    quantities: dict[str, float],
    setup_qty: float,
    running_qty: float,
    costs: dict[str, Any],
) -> list[str]:
    quantities_by_component = {
        "screeningTestCost": quantities["screened"],
        "tptRegimenCost": quantities["tpt_started_total"],
        "falsePositiveIncrementalCost": quantities["tpt_started_false_positive"],
        "activeTBDiseaseCost": quantities["active_tb_cases"],
        "programSetupCost": setup_qty,
        "programRunningCost": running_qty,
        "adrManagementCost": quantities["tpt_adr_stop_total"],
    }
    cost_keys = {component: key for component, (key, _) in COMPONENT_SPECS.items()}
    missing = []
    for component, value in values.items():
        if value is not None:
            continue
        if abs(float(quantities_by_component[component])) <= 1e-12:
            continue
        key = cost_keys[component]
        if key in costs.get("basisIssues", {}):
            missing.append(f"{component}: {costs['basisIssues'][key]}")
        else:
            missing.append(component)
    return missing


def _allocate_program_running_total(annual: pd.DataFrame, costs: dict[str, Any], metadata: dict[str, Any]) -> None:
    if annual.empty or costs.get("running") is None:
        return
    if costs.get("runningBasis") != "total_over_screening_programme" or costs.get("runningAllocation") != "proportional_to_screening_volume":
        return
    running_total = float(costs["running"])
    mask = annual["arm"].eq("intervention") & annual["includedInEconomicAnalysis"].astype(bool)
    screening_window = float(metadata.get("screeningWindowYears") or metadata.get("screeningWindow") or 0)
    mask &= pd.to_numeric(annual["modelYear"], errors="coerce").fillna(0) < screening_window
    for _, idx in annual[mask].groupby(ID_COLS, dropna=False).groups.items():
        screened = pd.to_numeric(annual.loc[idx, "screened"], errors="coerce").fillna(0.0)
        total_screened = float(screened.sum())
        annual.loc[idx, "programRunningCost"] = 0.0 if total_screened <= 0 else screened / total_screened * running_total
    _recompute_completeness_and_totals(annual)


def _recompute_completeness_and_totals(annual: pd.DataFrame) -> None:
    for idx, row in annual.iterrows():
        missing_cost = [col for col in COST_COMPONENTS if _num(row.get(col)) is None]
        missing_cost.extend(
            part for part in str(row.get("missingCostComponents") or "").split(";") if part and part not in missing_cost
        )
        missing_daly = [col for col in DALY_COMPONENTS if _num(row.get(col)) is None]
        missing_daly.extend(
            part for part in str(row.get("missingDALYComponents") or "").split(";") if part and part not in missing_daly
        )
        annual.at[idx, "costComplete"] = not missing_cost
        annual.at[idx, "dalyComplete"] = not missing_daly
        annual.at[idx, "missingCostComponents"] = ";".join(missing_cost)
        annual.at[idx, "missingDALYComponents"] = ";".join(missing_daly)
        cost_total = None if missing_cost else float(sum(_num(row.get(col)) or 0.0 for col in COST_COMPONENTS))
        daly_total = None if missing_daly else float(sum(_num(row.get(col)) or 0.0 for col in DALY_COMPONENTS))
        annual.at[idx, "totalUndiscountedCost"] = cost_total
        annual.at[idx, "totalDiscountedCost"] = None if cost_total is None else cost_total * float(row.get("discountFactor", 1.0))
        annual.at[idx, "totalUndiscountedDALYs"] = daly_total
        annual.at[idx, "totalDiscountedDALYs"] = None if daly_total is None else daly_total * float(row.get("discountFactor", 1.0))


def _replicate_results(annual: pd.DataFrame, threshold: float | None) -> pd.DataFrame:
    included = annual[annual["includedInEconomicAnalysis"].astype(bool)].copy()
    arm_rows = []
    group_cols = ID_COLS + ["arm"]
    for keys, group in included.groupby(group_cols, dropna=False):
        cost_complete = bool(group["costComplete"].all())
        daly_complete = bool(group["dalyComplete"].all())
        arm_rows.append(
            {
                **dict(zip(group_cols, keys)),
                "armCostComplete": cost_complete,
                "armDALYComplete": daly_complete,
                "armMissingCostComponents": _join_unique(group["missingCostComponents"]),
                "armMissingDALYComponents": _join_unique(group["missingDALYComponents"]),
                "totalCost": _strict_sum(group["totalDiscountedCost"]) if cost_complete else None,
                "totalDALYs": _strict_sum(group["totalDiscountedDALYs"]) if daly_complete else None,
                "undiscountedTotalCost": _strict_sum(group["totalUndiscountedCost"]) if cost_complete else None,
                "undiscountedTotalDALYs": _strict_sum(group["totalUndiscountedDALYs"]) if daly_complete else None,
                "activeTBCases": float(group["active_tb_cases"].sum()),
                "activeTBCasesPrevented": float(group["active_tb_cases_prevented"].sum()),
                "infectionsEffectivelyTreated": float(group["infection_effectively_treated_total"].sum()),
                **_component_totals(group),
            }
        )
    arms = pd.DataFrame(arm_rows)
    rows = []
    if arms.empty:
        return _frame(rows)
    for keys, pair in arms.groupby(ID_COLS, dropna=False):
        comp = pair[pair["arm"] == "comparator"]
        inter = pair[pair["arm"] == "intervention"]
        if comp.empty or inter.empty:
            continue
        c = comp.iloc[0]
        i = inter.iloc[0]
        cost_pair_complete = bool(c["armCostComplete"] and i["armCostComplete"])
        daly_pair_complete = bool(c["armDALYComplete"] and i["armDALYComplete"])
        economic_pair_complete = bool(cost_pair_complete and daly_pair_complete)
        inc = _subtract_if_number(i["totalCost"], c["totalCost"]) if cost_pair_complete else None
        dalys_averted = _subtract_if_number(c["totalDALYs"], i["totalDALYs"]) if daly_pair_complete else None
        classification = classify_incremental_result(inc, dalys_averted) if economic_pair_complete else "incomplete / not calculated"
        interpretable_icer = classification == INTERPRETABLE_ICER_CLASSIFICATION
        nmb = None if threshold is None or not economic_pair_complete else threshold * dalys_averted - inc
        rows.append(
            {
                **dict(zip(ID_COLS, keys)),
                "costPairComplete": cost_pair_complete,
                "dalyPairComplete": daly_pair_complete,
                "economicPairComplete": economic_pair_complete,
                "exclusionReasons": _join_unique(
                    [
                        c.get("armMissingCostComponents"),
                        i.get("armMissingCostComponents"),
                        c.get("armMissingDALYComponents"),
                        i.get("armMissingDALYComponents"),
                    ]
                ),
                "comparatorCost": c["totalCost"] if cost_pair_complete else None,
                "interventionCost": i["totalCost"] if cost_pair_complete else None,
                "incrementalCost": inc,
                "comparatorDALYs": c["totalDALYs"] if daly_pair_complete else None,
                "interventionDALYs": i["totalDALYs"] if daly_pair_complete else None,
                "dalysAverted": dalys_averted,
                "activeTBCasesPrevented": i["activeTBCasesPrevented"],
                "infectionsEffectivelyTreated": i["infectionsEffectivelyTreated"],
                "costPerActiveTBCasePrevented": _div(inc, i["activeTBCasesPrevented"]) if cost_pair_complete else None,
                "costPerInfectionEffectivelyTreated": _div(inc, i["infectionsEffectivelyTreated"]) if cost_pair_complete else None,
                "replicateICER": _div(inc, dalys_averted) if interpretable_icer else None,
                "primaryICERInterpretable": interpretable_icer,
                "netMonetaryBenefit": nmb,
                "classification": classification,
                **_paired_component_totals(c, i),
            }
        )
    return _frame(rows)


def _component_totals(group: pd.DataFrame) -> dict[str, Any]:
    out = {}
    for component in COST_COMPONENTS:
        out[f"{component}Discounted"] = _strict_sum(group[component] * group["discountFactor"])
        out[f"{component}Undiscounted"] = _strict_sum(group[component])
    out["totalProgrammeCostDiscounted"] = _strict_sum(sum(group[col] for col in PROGRAM_COMPONENTS) * group["discountFactor"])
    out["totalProgrammeCostUndiscounted"] = _strict_sum(sum(group[col] for col in PROGRAM_COMPONENTS))
    return out


def _paired_component_totals(comp: pd.Series, inter: pd.Series) -> dict[str, Any]:
    out = {}
    for side, row in {"comparator": comp, "intervention": inter}.items():
        for component in COST_COMPONENTS:
            out[f"{side}_{component}Discounted"] = row.get(f"{component}Discounted")
            out[f"{side}_{component}Undiscounted"] = row.get(f"{component}Undiscounted")
        out[f"{side}_totalProgrammeCostDiscounted"] = row.get("totalProgrammeCostDiscounted")
        out[f"{side}_totalProgrammeCostUndiscounted"] = row.get("totalProgrammeCostUndiscounted")
        out[f"{side}_activeTBDiseaseCareDiscounted"] = row.get("activeTBDiseaseCostDiscounted")
        out[f"{side}_totalArmCostDiscounted"] = row.get("totalCost")
        out[f"{side}_totalArmCostUndiscounted"] = row.get("undiscountedTotalCost")
    return out


def _summaries(reps: pd.DataFrame, metadata: dict[str, Any], threshold: float | None) -> pd.DataFrame:
    rows = []
    metrics = [
        "comparatorCost",
        "interventionCost",
        "incrementalCost",
        "comparatorDALYs",
        "interventionDALYs",
        "dalysAverted",
        "activeTBCasesPrevented",
        "netMonetaryBenefit",
    ]
    if reps.empty:
        return _frame(rows)
    for rate, group in reps.groupby("discountRate", dropna=False):
        complete = group[group["economicPairComplete"].astype(bool)]
        total_pairs = int(len(group))
        complete_pairs = int(len(complete))
        rows.append(
            {
                "discountRate": rate,
                "metric": "pairedReplicateCompleteness",
                "totalPairedReplicates": total_pairs,
                "completePairedReplicates": complete_pairs,
                "excludedPairedReplicates": total_pairs - complete_pairs,
                "exclusionReasons": _join_unique(group["exclusionReasons"]),
            }
        )
        mean_inc = float(complete["incrementalCost"].mean()) if complete_pairs else None
        mean_daly = float(complete["dalysAverted"].mean()) if complete_pairs else None
        classification = classify_incremental_result(mean_inc, mean_daly) if complete_pairs else "incomplete / not calculated"
        primary_icer = _div(mean_inc, mean_daly) if classification == INTERPRETABLE_ICER_CLASSIFICATION else None
        rows.append(
            {
                "discountRate": rate,
                "metric": "primaryICER_ratioOfMeans",
                "mean": primary_icer,
                "classification": classification,
                "n": complete_pairs,
                "totalPairedReplicates": total_pairs,
                "completePairedReplicates": complete_pairs,
                "excludedPairedReplicates": total_pairs - complete_pairs,
                "intervalLabel": "deterministic expected value" if total_pairs == 1 else "simulation distribution across replicates",
            }
        )
        if threshold is not None and metadata.get("modelType") != "expected_value":
            valid_nmb = complete[pd.to_numeric(complete["netMonetaryBenefit"], errors="coerce").notna()]
            denominator = int(len(valid_nmb))
            numerator = int((valid_nmb["netMonetaryBenefit"] > 0).sum()) if denominator else 0
            rows.append(
                {
                    "discountRate": rate,
                    "metric": "probabilityPositiveNMB_fixedParameterSimulation",
                    "mean": None if denominator == 0 else numerator / denominator,
                    "n": denominator,
                    "numerator": numerator,
                    "denominator": denominator,
                    "intervalLabel": (
                        "probability of positive NMB across finite-population "
                        "simulation replicates under fixed parameter assumptions"
                    ),
                }
            )
        for metric in metrics:
            if metric not in complete:
                continue
            values = pd.to_numeric(complete[metric], errors="coerce").dropna()
            if values.empty:
                continue
            rows.append(
                {
                    "discountRate": rate,
                    "metric": metric,
                    "n": int(values.count()),
                    "mean": float(values.mean()),
                    "sd": None if len(values) == 1 else float(values.std(ddof=1)),
                    "median": float(values.median()),
                    "p2_5": float(np.percentile(values, 2.5)),
                    "p97_5": float(np.percentile(values, 97.5)),
                    "min": float(values.min()),
                    "max": float(values.max()),
                    "intervalLabel": "deterministic expected value" if len(values) == 1 else "simulation distribution across replicates",
                }
            )
    return _frame(rows)


def classify_incremental_result(cost: Any, dalys: Any, tol: float = 1e-9) -> str:
    cost = _num(cost)
    dalys = _num(dalys)
    if cost is None or dalys is None:
        return "incomplete / not calculated"
    if abs(cost) <= tol and abs(dalys) <= tol:
        return "no material difference"
    if cost < -tol and dalys > tol:
        return "dominant"
    if cost > tol and dalys < -tol:
        return "dominated"
    if cost < -tol and dalys < -tol:
        return "cost saving with health loss"
    if abs(cost) <= tol and dalys > tol:
        return "health gain at no additional cost"
    if cost < -tol and abs(dalys) <= tol:
        return "cost saving with no material health difference"
    if cost > tol and abs(dalys) <= tol:
        return "increased cost with no material health difference"
    if cost > tol and dalys > tol:
        return INTERPRETABLE_ICER_CLASSIFICATION
    return "ICER not interpretable because health gain is zero or negative"


def _validate_economic_result(annual: pd.DataFrame, reps: pd.DataFrame, summaries: pd.DataFrame, threshold: float | None) -> dict[str, Any]:
    errors = []
    warnings = []
    if annual.empty:
        errors.append({"field": "annualByArm", "message": "Annual economic ledger is empty."})
    if reps.empty:
        errors.append({"field": "replicateResults", "message": "No paired replicate economic results were calculated."})
    if not annual.empty:
        for _, row in annual.iterrows():
            if bool(row.get("costComplete")):
                expected = sum(_num(row.get(col)) or 0.0 for col in COST_COMPONENTS)
                _check_close(errors, "annual.totalUndiscountedCost", row.get("totalUndiscountedCost"), expected)
                _check_close(errors, "annual.totalDiscountedCost", row.get("totalDiscountedCost"), expected * float(row.get("discountFactor")))
            if bool(row.get("dalyComplete")):
                expected = sum(_num(row.get(col)) or 0.0 for col in DALY_COMPONENTS)
                _check_close(errors, "annual.totalUndiscountedDALYs", row.get("totalUndiscountedDALYs"), expected)
                _check_close(errors, "annual.totalDiscountedDALYs", row.get("totalDiscountedDALYs"), expected * float(row.get("discountFactor")))
            if not bool(row.get("includedInEconomicAnalysis")) and (
                _num(row.get("totalDiscountedCost")) or _num(row.get("totalDiscountedDALYs"))
            ):
                warnings.append({"field": "annualByArm", "message": "Outside-horizon row retained for audit and excluded from replicate totals."})
    if not reps.empty:
        for _, row in reps.iterrows():
            if not bool(row.get("economicPairComplete")):
                for field in ["incrementalCost", "dalysAverted", "replicateICER", "netMonetaryBenefit"]:
                    if _num(row.get(field)) is not None:
                        errors.append({"field": field, "message": "Incomplete pair has a calculated aggregate result."})
            if bool(row.get("costPairComplete")):
                _check_close(errors, "replicate.incrementalCost", row.get("incrementalCost"), _subtract_if_number(row.get("interventionCost"), row.get("comparatorCost")))
            if bool(row.get("dalyPairComplete")):
                _check_close(errors, "replicate.dalysAverted", row.get("dalysAverted"), _subtract_if_number(row.get("comparatorDALYs"), row.get("interventionDALYs")))
            if threshold is not None and bool(row.get("economicPairComplete")):
                _check_close(errors, "replicate.netMonetaryBenefit", row.get("netMonetaryBenefit"), threshold * row["dalysAverted"] - row["incrementalCost"])
    total_pairs = int(len(reps)) if not reps.empty else 0
    complete_pairs = int(reps["economicPairComplete"].sum()) if not reps.empty and "economicPairComplete" in reps else 0
    return {
        "isValid": not errors,
        "structurallyValid": not errors,
        "economicallyComplete": bool(total_pairs and total_pairs == complete_pairs),
        "conclusionPermitted": bool(total_pairs and total_pairs == complete_pairs and threshold is not None),
        "totalPairedReplicates": total_pairs,
        "completePairedReplicates": complete_pairs,
        "excludedPairedReplicates": total_pairs - complete_pairs,
        "exclusionReasons": _join_unique(reps["exclusionReasons"]) if not reps.empty and "exclusionReasons" in reps else "",
        "errors": errors,
        "warnings": warnings,
    }


def _summary_rows(result: dict[str, Any]) -> list[dict[str, Any]]:
    rows = result["summaries"]
    return rows.to_dict(orient="records") if isinstance(rows, pd.DataFrame) else []


def _not_calculated(result: dict[str, Any]) -> list[str]:
    out = []
    reps = result.get("replicateResults")
    if isinstance(reps, pd.DataFrame) and not reps.empty:
        if not bool(reps["costPairComplete"].all()):
            out.extend(["totalCost", "incrementalCost", "ICER"])
        if not bool(reps["dalyPairComplete"].all()):
            out.extend(["totalDALYs", "dalysAverted", "ICER"])
    if any(item["field"].startswith("threshold") for item in result["unresolvedInputs"]):
        out.append("netMonetaryBenefit")
    return sorted(set(out))


def _attach_legacy_compatibility_fields(result: dict[str, Any], econ_config: dict[str, Any], costs: dict[str, Any]) -> None:
    annual = result["annualByArm"]
    primary = annual[annual["discountRate"] == PRIMARY_DISCOUNT_RATE] if isinstance(annual, pd.DataFrame) and not annual.empty else pd.DataFrame()
    comp = primary[primary["arm"] == "comparator"] if not primary.empty else pd.DataFrame()
    inter = primary[primary["arm"] == "intervention"] if not primary.empty else pd.DataFrame()

    def mean_sum(frame: pd.DataFrame, column: str) -> float | None:
        if frame.empty or column not in frame:
            return None
        frame = frame[frame["includedInEconomicAnalysis"].astype(bool)]
        grouped = frame.groupby(["replicateId", "pairedReplicateId", "replicateSeed"], dropna=False)[column].agg(
            lambda s: np.nan if s.isna().any() else s.sum()
        )
        values = pd.to_numeric(grouped, errors="coerce").dropna()
        return None if values.empty or len(values) < len(grouped) else float(values.mean())

    quantities = {
        "nScreened": mean_sum(inter, "screened"),
        "nTotalCoursesStarted": mean_sum(inter, "tpt_started_total"),
        "nFalsePositiveTreated": mean_sum(inter, "tpt_started_false_positive"),
        "nCuredInfection": mean_sum(inter, "infection_effectively_treated_total"),
        "nPreventedActiveTB": mean_sum(inter, "active_tb_cases_prevented"),
        "baselineActiveTBCases": mean_sum(comp, "active_tb_cases"),
        "interventionActiveTBCases": mean_sum(inter, "active_tb_cases"),
    }
    cost_values = {
        "testingCost": mean_sum(inter, "screeningTestCost"),
        "treatmentCost": mean_sum(inter, "tptRegimenCost"),
        "falsePositiveIncrementalCost": mean_sum(inter, "falsePositiveIncrementalCost"),
        "programSetupCost": mean_sum(inter, "programSetupCost"),
        "programRunningCost": mean_sum(inter, "programRunningCost"),
        "adrManagementCost": mean_sum(inter, "adrManagementCost"),
        "baselineTBDiseaseCost": mean_sum(comp, "activeTBDiseaseCost"),
        "interventionTBDiseaseCost": mean_sum(inter, "activeTBDiseaseCost"),
    }
    cost_values["tbDiseaseCostsAverted"] = _subtract_if_number(cost_values["baselineTBDiseaseCost"], cost_values["interventionTBDiseaseCost"])
    cost_values["totalProgramCost"] = _sum_optional([cost_values[col] for col in ["testingCost", "treatmentCost", "falsePositiveIncrementalCost", "programSetupCost", "programRunningCost", "adrManagementCost"]])
    cost_values["netCostVsBaseline"] = _subtract_if_number(cost_values["totalProgramCost"], cost_values["tbDiseaseCostsAverted"])
    legacy_status = {
        "missingInputs": [item["field"] for item in result["unresolvedInputs"]],
        "notCalculated": [field for field, value in cost_values.items() if value is None],
        "messages": [],
        "partialCalculations": [],
    }
    if costs.get("false_positive") is None:
        legacy_status["messages"].append(
            "No false-positive incremental unit cost supplied; no extra false-positive cost is included."
        )
    result["inputs"] = deepcopy(econ_config)
    result["costNormalisation"] = deepcopy(econ_config.get("costNormalisation", {}))
    result["discounting"] = deepcopy(econ_config.get("discounting", {}))
    result["healthOutcome"] = deepcopy(econ_config.get("healthOutcome", {}))
    result["threshold"] = deepcopy(econ_config.get("threshold", {}))
    result["scopeStatement"] = result["metadata"].get("scopeStatement", DIRECT_EFFECTS_SCOPE_STATEMENT)
    result["strategy"] = {
        "testType": _config_value(econ_config, result["metadata"], "testType"),
        "regimen": _config_value(econ_config, result["metadata"], "regimen"),
    }
    result["quantities"] = quantities
    result["unitCosts"] = {
        "testPerPerson": costs.get("test"),
        "treatmentPerStarted": costs.get("regimen"),
        "falsePositiveIncrementalPerPerson": costs.get("false_positive"),
        "activeTBDiseasePerCase": costs.get("active_tb"),
    }
    result["costs"] = cost_values
    result["costEffectiveness"] = {
        "costPerInfectionCured": _div(cost_values["netCostVsBaseline"], quantities["nCuredInfection"]),
        "costPerTBCasePrevented": _div(cost_values["netCostVsBaseline"], quantities["nPreventedActiveTB"]),
    }
    result["_legacyCompatibilityStatus"] = legacy_status


def _reviewed_numeric(daly: dict[str, Any], field: str, bounds: tuple[float | None, float | None], unresolved: list[dict[str, Any]]) -> float | None:
    rec = daly.get(field) or {}
    value = _num(rec.get("value") if isinstance(rec, dict) else rec)
    status = rec.get("status") if isinstance(rec, dict) else ""
    source = rec.get("source") if isinstance(rec, dict) else ""
    low, high = bounds
    if value is None or status not in NUMERIC_REVIEWED_STATUSES or not str(source or "").strip():
        unresolved.append(_unresolved(f"dalyAssumptions.{field}", f"{field} requires reviewed numeric value and source."))
        return value
    if (low is not None and value < low) or (high is not None and value > high):
        unresolved.append(_unresolved(f"dalyAssumptions.{field}", f"{field} is outside the valid range."))
    return value


def _reviewed_exclusion(daly: dict[str, Any], status_field: str, rationale_field: str, field: str, unresolved: list[dict[str, Any]]) -> None:
    if daly.get(status_field) != EXCLUSION_REVIEWED_STATUS or not str(daly.get(rationale_field) or "").strip():
        unresolved.append(_unresolved(field, f"{field} requires reviewed exclusion status and rationale."))


def _reviewed_cost_exclusion(assumptions: dict[str, Any], cost_item_id: str) -> dict[str, Any] | None:
    rec = (assumptions.get("optionalCostExclusions") or {}).get(cost_item_id) or {}
    if rec.get("status") == EXCLUSION_REVIEWED_STATUS and str(rec.get("rationale") or "").strip():
        return rec
    return None


def _basis(item: dict[str, Any]) -> str:
    return str((item.get("resourceUse") or {}).get("costBasis") or "")


def _event_columns(frame: pd.DataFrame) -> list[str]:
    non_events = {
        "scenarioId",
        "scenarioVersion",
        "populationPreset",
        "modelType",
        "backend",
        "arm",
        "replicateId",
        "pairedReplicateId",
        "replicateSeed",
        "valueType",
        "screeningWindowYears",
        "followUpHorizonYears",
        "contractVersion",
        "modelVersion",
        "modelYear",
        "timeInterval",
        "withinFollowUp",
    }
    return [col for col in frame.columns if col not in non_events]


def _strict_sum(values: Any) -> float | None:
    series = pd.to_numeric(pd.Series(values), errors="coerce")
    if series.isna().any():
        return None
    return float(series.sum())


def _frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(rows).astype(object).where(pd.notna(pd.DataFrame(rows)), None) if rows else pd.DataFrame(rows)


def _join_unique(values: Any) -> str:
    parts: list[str] = []
    if isinstance(values, pd.Series):
        iterable = values.to_list()
    elif isinstance(values, (list, tuple, set)):
        iterable = list(values)
    else:
        iterable = [values]
    for value in iterable:
        for part in str(value or "").split(";"):
            part = part.strip()
            if part and part not in parts:
                parts.append(part)
    return ";".join(parts)


def _sum_optional(values: list[Any]) -> float | None:
    total = 0.0
    for value in values:
        value = _num(value)
        if value is None:
            return None
        total += value
    return total


def _subtract_if_number(a: Any, b: Any) -> float | None:
    a = _num(a)
    b = _num(b)
    if a is None or b is None:
        return None
    return a - b


def _check_close(errors: list[dict[str, str]], field: str, observed: Any, expected: Any, tol: float = 1e-8) -> None:
    observed = _num(observed)
    expected = _num(expected)
    if observed is None or expected is None:
        return
    if abs(observed - expected) > tol:
        errors.append({"field": field, "message": f"{observed} does not reconcile with {expected}."})


def _extract_event_ledger(results: dict[str, Any]) -> dict[str, Any]:
    for candidate in (
        results.get("eventLedger"),
        results.get("technical", {}).get("eventLedger") if isinstance(results.get("technical"), dict) else None,
        results.get("bundle", {}).get("technical", {}).get("eventLedger") if isinstance(results.get("bundle"), dict) else None,
    ):
        if isinstance(candidate, dict) and candidate.get("metadata", {}).get("contractVersion") == EVENT_LEDGER_CONTRACT_VERSION:
            return candidate
    return {}


def _config_value(econ_config: dict[str, Any], metadata: dict[str, Any], field: str) -> Any:
    return econ_config.get("scenarioConfig", {}).get(field) or metadata.get(field)


def _event(row: pd.Series, name: str) -> float:
    return float(row.get(name, 0.0) or 0.0)


def _mul(a: Any, b: Any) -> float | None:
    a = _num(a)
    b = _num(b)
    if a == 0:
        return 0.0
    return None if a is None or b is None else a * b


def _mul3(a: Any, b: Any, c: Any) -> float | None:
    first = _mul(a, b)
    return None if first is None else _mul(first, c)


def _div(a: Any, b: Any) -> float | None:
    a = _num(a)
    b = _num(b)
    if a is None or b is None or abs(b) < 1e-12:
        return None
    return a / b


def _num(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str) and value == "":
        return None
    if isinstance(value, (list, tuple)) and len(value) == 0:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _unresolved(field: str, message: str) -> dict[str, Any]:
    return {"field": field, "severity": "unresolved", "message": message}


def _first_value(frame: pd.DataFrame, column: str) -> Any:
    return None if frame.empty or column not in frame else frame[column].iloc[0]


def _deep_update(target: dict[str, Any], updates: dict[str, Any]) -> None:
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            _deep_update(target[key], value)
        else:
            target[key] = value
