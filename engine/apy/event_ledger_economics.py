from __future__ import annotations

from copy import deepcopy
import math
from typing import Any

import numpy as np
import pandas as pd

from engine.apy.costing import TARGET_CURRENCY, TARGET_PRICE_YEAR, normalise_cost_table, valid_converted_cost
from engine.apy.event_ledger import EVENT_LEDGER_CONTRACT_VERSION
from engine.apy.scenario import DIRECT_EFFECTS_SCOPE_STATEMENT


HEALTH_ECONOMICS_CONTRACT_VERSION = "ltbi_health_economics_results_v2"
PRIMARY_DISCOUNT_RATE = 0.03
COMPARISON_DISCOUNT_RATE = 0.0
VALID_REVIEWED_STATUSES = {"configured_reviewed", "model_derived_reviewed", "reviewed_exclusion"}


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
        "dalyLossPerADRStop": _assumption("DALY loss per ADR-related stop"),
        "includePostTBSequelae": False,
        "postTBSequelaeStatus": "unresolved",
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
    metadata = deepcopy(ledger.get("metadata") or {})
    cost_items = normalise_cost_table(econ_config.get("costItems") or [])
    assumptions = _build_assumptions(econ_config)
    unresolved: list[dict[str, Any]] = []
    warnings: list[str] = []
    cost_lookup = _cost_lookup(cost_items, econ_config, metadata, assumptions["metadata"], unresolved)
    daly_inputs = _resolve_daly_inputs(assumptions["daly"], unresolved)
    threshold = _resolve_threshold(econ_config.get("threshold") or {}, assumptions, unresolved)
    annual_events = _events_wide(ledger["annualEvents"])
    annual_rows = []
    for _, event_row in annual_events.iterrows():
        for rate in discount_rates:
            annual_rows.append(
                _annual_economic_row(
                    event_row,
                    metadata,
                    cost_lookup,
                    daly_inputs,
                    assumptions,
                    float(rate),
                    unresolved,
                )
            )
    annual = pd.DataFrame(annual_rows)
    _allocate_program_running_total(annual, cost_lookup, metadata)
    replicate_results = _replicate_results(annual, threshold)
    summaries = _summaries(replicate_results, metadata, threshold)
    provisional = bool(metadata.get("ltbiStateProvisional") or metadata.get("ltbiStateWarning"))
    if provisional:
        warnings.append("Epidemiological inputs are provisional; no clinician-ready cost-effectiveness conclusion is produced.")
    if unresolved:
        warnings.append("One or more economic inputs are unresolved; unavailable components are not replaced with zero.")
    result = {
        "available": True,
        "source": "event_ledger_health_economics_v2",
        "contractVersion": HEALTH_ECONOMICS_CONTRACT_VERSION,
        "metadata": {
            "economicContractVersion": HEALTH_ECONOMICS_CONTRACT_VERSION,
            "eventLedgerContractVersion": metadata.get("contractVersion"),
            "scenarioId": metadata.get("scenarioId"),
            "modelType": metadata.get("modelType"),
            "valueType": _first_value(annual, "valueType"),
            "perspective": assumptions["metadata"].get("perspective"),
            "targetCurrency": assumptions["metadata"].get("targetCurrency", TARGET_CURRENCY),
            "targetPriceYear": assumptions["metadata"].get("targetPriceYear", TARGET_PRICE_YEAR),
            "scopeStatement": DIRECT_EFFECTS_SCOPE_STATEMENT,
            "ltbiStateProvisional": metadata.get("ltbiStateProvisional"),
            "ltbiStateDevelopmentCompatibilityMode": metadata.get("ltbiStateAssumptionStatus") == "unresolved_development_compatibility",
            "primaryDiscountRate": PRIMARY_DISCOUNT_RATE,
            "comparisonDiscountRate": COMPARISON_DISCOUNT_RATE,
            "isProvisional": bool(provisional or unresolved),
            "conclusionPermitted": False if provisional or unresolved else True,
        },
        "assumptions": assumptions,
        "costItems": cost_items,
        "annualByArm": annual,
        "replicateResults": replicate_results,
        "summaries": summaries,
        "validation": _validate_economic_result(annual, replicate_results),
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
        "isComplete": not _blocking_unresolved(unresolved),
        "missingInputs": [item["field"] for item in unresolved],
        "notCalculated": sorted(set(_not_calculated(result) + result["_legacyCompatibilityStatus"]["notCalculated"])),
        "messages": warnings + result["_legacyCompatibilityStatus"]["messages"],
        "partialCalculations": result["_legacyCompatibilityStatus"].get("partialCalculations", []),
        "validationReport": result["validation"],
    }
    return result


def _assumption(label: str) -> dict[str, Any]:
    return {"label": label, "value": None, "source": "", "status": "unresolved", "notes": "", "provisional": True}


def _build_assumptions(econ_config: dict[str, Any]) -> dict[str, Any]:
    assumptions = {
        "metadata": deepcopy(econ_config.get("metadata") or {}),
        "discounting": deepcopy(econ_config.get("discounting") or {}),
        "threshold": deepcopy(econ_config.get("threshold") or {}),
        "daly": default_daly_assumptions(),
    }
    if isinstance(econ_config.get("dalyAssumptions"), dict):
        _deep_update(assumptions["daly"], econ_config["dalyAssumptions"])
    assumptions["metadata"].setdefault("targetCurrency", TARGET_CURRENCY)
    assumptions["metadata"].setdefault("targetPriceYear", TARGET_PRICE_YEAR)
    assumptions["metadata"].setdefault("perspective", "Australian health-system perspective")
    return assumptions


def _cost_lookup(
    cost_items: list[dict[str, Any]],
    econ_config: dict[str, Any],
    metadata: dict[str, Any],
    analysis_metadata: dict[str, Any],
    unresolved: list[dict[str, Any]],
) -> dict[str, Any]:
    test_type = str(_config_value(econ_config, metadata, "testType") or "IGRA").upper()
    regimen = str(_config_value(econ_config, metadata, "regimen") or "3HP").upper()
    ids = {
        "test": {"IGRA": "test_igra", "TST": "test_tst"}.get(test_type),
        "regimen": {"3HP": "regimen_3hp", "4R": "regimen_4r", "3HR": "regimen_3hr", "6H": "regimen_6h", "9H": "regimen_9h"}.get(regimen),
        "active_tb": "active_tb_disease",
        "false_positive": "false_positive_incremental",
        "setup": "program_setup",
        "running": "program_running",
        "adr": "tpt_adr_management",
    }
    by_id = {item.get("costItemId"): item for item in cost_items}
    out: dict[str, Any] = {"ids": ids, "items": by_id}
    target_currency = analysis_metadata.get("targetCurrency", TARGET_CURRENCY)
    target_year = analysis_metadata.get("targetPriceYear", TARGET_PRICE_YEAR)
    for name, item_id in ids.items():
        if item_id is None:
            unresolved.append(_unresolved(f"costItems.{name}", "No cost item mapping is available."))
            out[name] = None
            continue
        item = by_id.get(item_id, {})
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
        if value is None and name != "adr":
            unresolved.append(_unresolved(f"costItems.{item_id}", f"{item_id} lacks valid converted target-year cost."))
    running_item = by_id.get("program_running", {})
    running_basis = (running_item.get("resourceUse") or {}).get("costBasis")
    if out.get("running") is not None and running_basis not in {"annual_during_screening_window", "total_over_screening_programme"}:
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
    fields = ["activeTBDisabilityWeight", "activeTBDurationYears", "tbCaseFatalityRisk", "yllPerTBDeath"]
    out = {}
    for field in fields:
        rec = daly.get(field) or {}
        value = _num(rec.get("value") if isinstance(rec, dict) else rec)
        status = rec.get("status") if isinstance(rec, dict) else ""
        source = rec.get("source") if isinstance(rec, dict) else ""
        if value is None or status not in VALID_REVIEWED_STATUSES or not str(source or "").strip():
            unresolved.append(_unresolved(f"dalyAssumptions.{field}", f"{field} requires reviewed value and source."))
        out[field] = value
    for field in ["dalyLossPerTPTStarted", "dalyLossPerADRStop"]:
        rec = daly.get(field) or {}
        out[field] = _num(rec.get("value") if isinstance(rec, dict) else rec)
    include_tpt = bool(daly.get("includeTPTHealthLoss"))
    if not include_tpt and daly.get("tptHealthLossExclusionStatus") != "reviewed_exclusion":
        unresolved.append(_unresolved("dalyAssumptions.includeTPTHealthLoss", "TPT health loss requires reviewed inclusion values or reviewed exclusion rationale."))
    out["includeTPTHealthLoss"] = include_tpt
    return out


def _resolve_threshold(threshold: dict[str, Any], assumptions: dict[str, Any], unresolved: list[dict[str, Any]]) -> float | None:
    value = _num(threshold.get("value"))
    status = str(threshold.get("status") or "")
    if value is None:
        unresolved.append(_unresolved("threshold.value", "Willingness-to-pay threshold is unresolved; NMB unavailable."))
        return None
    if not threshold.get("currency") or threshold.get("referenceYear") in (None, "", []) or not threshold.get("source") or status not in VALID_REVIEWED_STATUSES:
        unresolved.append(_unresolved("threshold.provenance", "Threshold requires currency, year, source and reviewed status."))
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


def _annual_economic_row(row: pd.Series, metadata: dict[str, Any], costs: dict[str, Any], daly: dict[str, Any], assumptions: dict[str, Any], rate: float, unresolved: list[dict[str, Any]]) -> dict[str, Any]:
    year = int(row.get("modelYear", 0))
    factor = 1.0 / ((1.0 + rate) ** year)
    arm = row.get("arm")
    screened = _event(row, "screened")
    starts = _event(row, "tpt_started_total")
    false_starts = _event(row, "tpt_started_false_positive")
    adr_stops = _event(row, "tpt_adr_stop_total")
    active_tb = _event(row, "active_tb_cases")
    setup = costs.get("setup") if arm == "intervention" and year == 0 else 0.0
    running = _running_cost(row, costs, metadata) if arm == "intervention" else 0.0
    components = {
        "screeningTestCost": _mul(screened, costs.get("test")),
        "tptRegimenCost": _mul(starts, costs.get("regimen")),
        "falsePositiveIncrementalCost": _mul(false_starts, costs.get("false_positive")),
        "activeTBDiseaseCost": _mul(active_tb, costs.get("active_tb")),
        "programSetupCost": setup,
        "programRunningCost": running,
        "adrManagementCost": _mul(adr_stops, costs.get("adr")) if costs.get("adr") is not None else 0.0,
    }
    cost_complete = all(v is not None for v in components.values())
    total_cost = None if not cost_complete else float(sum(components.values()))
    yld = _mul3(active_tb, daly.get("activeTBDisabilityWeight"), daly.get("activeTBDurationYears"))
    deaths = _mul(active_tb, daly.get("tbCaseFatalityRisk"))
    yll = _mul(deaths, daly.get("yllPerTBDeath"))
    tpt_loss = _mul(starts, daly.get("dalyLossPerTPTStarted")) if daly.get("includeTPTHealthLoss") else 0.0
    adr_loss = _mul(adr_stops, daly.get("dalyLossPerADRStop")) if daly.get("dalyLossPerADRStop") is not None else 0.0
    daly_values = [yld, yll, tpt_loss, adr_loss]
    total_dalys = None if any(v is None for v in daly_values) else float(sum(daly_values))
    out = {
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
        "discountRate": rate,
        "discountFactor": factor,
        "screened": screened,
        "tpt_started_total": starts,
        "tpt_started_false_positive": false_starts,
        "tpt_adr_stop_total": adr_stops,
        "active_tb_cases": active_tb,
        "active_tb_cases_prevented": _event(row, "active_tb_cases_prevented"),
        "infection_effectively_treated_total": _event(row, "infection_effectively_treated_total"),
        **components,
        "activeTBYLD": yld,
        "activeTBYLL": yll,
        "tptHealthLossDALYs": tpt_loss,
        "adrHealthLossDALYs": adr_loss,
        "totalUndiscountedCost": total_cost,
        "totalDiscountedCost": None if total_cost is None else total_cost * factor,
        "totalUndiscountedDALYs": total_dalys,
        "totalDiscountedDALYs": None if total_dalys is None else total_dalys * factor,
    }
    return out


def _running_cost(row: pd.Series, costs: dict[str, Any], metadata: dict[str, Any]) -> float | None:
    value = costs.get("running")
    if value is None:
        return None
    basis = costs.get("runningBasis")
    if basis == "annual_during_screening_window":
        return float(value) if float(row.get("modelYear", 0)) < float(metadata.get("screeningWindowYears") or metadata.get("screeningWindow") or 0) else 0.0
    if basis == "total_over_screening_programme" and costs.get("runningAllocation") == "proportional_to_screening_volume":
        return 0.0
    return None


def _allocate_program_running_total(annual: pd.DataFrame, costs: dict[str, Any], metadata: dict[str, Any]) -> None:
    if annual.empty or costs.get("running") is None:
        return
    if costs.get("runningBasis") != "total_over_screening_programme" or costs.get("runningAllocation") != "proportional_to_screening_volume":
        return
    running_total = float(costs["running"])
    id_cols = ["discountRate", "replicateId", "pairedReplicateId", "replicateSeed"]
    mask = annual["arm"].eq("intervention")
    screening_window = float(metadata.get("screeningWindowYears") or metadata.get("screeningWindow") or 0)
    mask &= pd.to_numeric(annual["modelYear"], errors="coerce").fillna(0) < screening_window
    for _, idx in annual[mask].groupby(id_cols, dropna=False).groups.items():
        screened = pd.to_numeric(annual.loc[idx, "screened"], errors="coerce").fillna(0.0)
        total_screened = float(screened.sum())
        annual.loc[idx, "programRunningCost"] = 0.0 if total_screened <= 0 else screened / total_screened * running_total
    _recompute_total_costs(annual)


def _recompute_total_costs(annual: pd.DataFrame) -> None:
    component_cols = [
        "screeningTestCost",
        "tptRegimenCost",
        "falsePositiveIncrementalCost",
        "activeTBDiseaseCost",
        "programSetupCost",
        "programRunningCost",
        "adrManagementCost",
    ]
    for idx, row in annual.iterrows():
        components = [_num(row.get(col)) for col in component_cols]
        total = None if any(value is None for value in components) else float(sum(components))
        annual.at[idx, "totalUndiscountedCost"] = total
        annual.at[idx, "totalDiscountedCost"] = None if total is None else total * float(row.get("discountFactor", 1.0))


def _replicate_results(annual: pd.DataFrame, threshold: float | None) -> pd.DataFrame:
    grouped = annual.groupby(["discountRate", "modelType", "valueType", "replicateId", "pairedReplicateId", "replicateSeed", "arm"], dropna=False).agg(
        totalCost=("totalDiscountedCost", lambda s: s.sum(min_count=1)),
        totalDALYs=("totalDiscountedDALYs", lambda s: s.sum(min_count=1)),
        activeTBCases=("active_tb_cases", "sum"),
        activeTBCasesPrevented=("active_tb_cases_prevented", "sum"),
        infectionsEffectivelyTreated=("infection_effectively_treated_total", "sum"),
        programCost=("screeningTestCost", lambda s: s.sum(min_count=1)),
        activeTBDiseaseCost=("activeTBDiseaseCost", lambda s: s.sum(min_count=1)),
    ).reset_index()
    rows = []
    for keys, pair in grouped.groupby(["discountRate", "modelType", "valueType", "replicateId", "pairedReplicateId", "replicateSeed"], dropna=False):
        comp = pair[pair["arm"] == "comparator"]
        inter = pair[pair["arm"] == "intervention"]
        if comp.empty or inter.empty:
            continue
        c = comp.iloc[0]
        i = inter.iloc[0]
        inc = i["totalCost"] - c["totalCost"]
        dalys_averted = c["totalDALYs"] - i["totalDALYs"]
        nmb = None if threshold is None else threshold * dalys_averted - inc
        rows.append({
            "discountRate": keys[0],
            "modelType": keys[1],
            "valueType": keys[2],
            "replicateId": keys[3],
            "pairedReplicateId": keys[4],
            "replicateSeed": keys[5],
            "comparatorCost": c["totalCost"],
            "interventionCost": i["totalCost"],
            "incrementalCost": inc,
            "comparatorDALYs": c["totalDALYs"],
            "interventionDALYs": i["totalDALYs"],
            "dalysAverted": dalys_averted,
            "activeTBCasesPrevented": i["activeTBCasesPrevented"],
            "infectionsEffectivelyTreated": i["infectionsEffectivelyTreated"],
            "costPerActiveTBCasePrevented": _div(inc, i["activeTBCasesPrevented"]),
            "costPerInfectionEffectivelyTreated": _div(inc, i["infectionsEffectivelyTreated"]),
            "replicateICER": _div(inc, dalys_averted),
            "netMonetaryBenefit": nmb,
            "classification": classify_incremental_result(inc, dalys_averted),
        })
    return pd.DataFrame(rows)


def _summaries(reps: pd.DataFrame, metadata: dict[str, Any], threshold: float | None) -> pd.DataFrame:
    rows = []
    metrics = ["comparatorCost", "interventionCost", "incrementalCost", "comparatorDALYs", "interventionDALYs", "dalysAverted", "activeTBCasesPrevented", "netMonetaryBenefit"]
    for rate, group in reps.groupby("discountRate", dropna=False):
        mean_inc = float(group["incrementalCost"].mean())
        mean_daly = float(group["dalysAverted"].mean())
        primary_icer = _div(mean_inc, mean_daly)
        rows.append(
            {
                "discountRate": rate,
                "metric": "primaryICER_ratioOfMeans",
                "mean": primary_icer,
                "n": int(len(group)),
                "intervalLabel": (
                    "simulation distribution across replicates"
                    if len(group) > 1
                    else "deterministic expected value"
                ),
            }
        )
        if threshold is not None and "netMonetaryBenefit" in group:
            rows.append(
                {
                    "discountRate": rate,
                    "metric": "probabilityPositiveNMB_fixedParameterSimulation",
                    "mean": float((group["netMonetaryBenefit"] > 0).mean()),
                    "n": int(len(group)),
                    "intervalLabel": (
                        "probability of positive NMB across finite-population "
                        "simulation replicates under fixed parameter assumptions"
                    ),
                }
            )
        for metric in metrics:
            if metric not in group:
                continue
            values = pd.to_numeric(group[metric], errors="coerce").dropna()
            if values.empty:
                continue
            rows.append({
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
                "intervalLabel": "simulation distribution across replicates" if len(values) > 1 else "deterministic expected value",
            })
    return pd.DataFrame(rows)


def classify_incremental_result(cost: float, dalys: float, tol: float = 1e-9) -> str:
    if abs(cost) <= tol and abs(dalys) <= tol:
        return "no material difference"
    if cost < -tol and dalys > tol:
        return "dominant"
    if cost > tol and dalys < -tol:
        return "dominated"
    if cost < -tol and dalys < -tol:
        return "cost saving with health loss"
    if cost > tol and dalys > tol:
        return "increased cost with health gain"
    if dalys <= tol:
        return "ICER not interpretable because health gain is zero or negative"
    return "no material difference"


def _validate_economic_result(annual: pd.DataFrame, reps: pd.DataFrame) -> dict[str, Any]:
    errors = []
    if annual.empty:
        errors.append({"field": "annualByArm", "message": "Annual economic ledger is empty."})
    if reps.empty:
        errors.append({"field": "replicateResults", "message": "No paired replicate economic results were calculated."})
    return {"isValid": not errors, "errors": errors, "warnings": []}


def _summary_rows(result: dict[str, Any]) -> list[dict[str, Any]]:
    rows = result["summaries"]
    return rows.to_dict(orient="records") if isinstance(rows, pd.DataFrame) else []


def _not_calculated(result: dict[str, Any]) -> list[str]:
    out = []
    blocking = _blocking_unresolved(result["unresolvedInputs"])
    if any(str(item["field"]).startswith("costItems") for item in blocking):
        out.extend(["totalCost", "ICER"])
    if any(str(item["field"]).startswith("dalyAssumptions") for item in blocking):
        out.extend(["totalDALYs", "ICER"])
    if any(item["field"].startswith("threshold") for item in result["unresolvedInputs"]):
        out.append("netMonetaryBenefit")
    return sorted(set(out))


def _blocking_unresolved(unresolved: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [item for item in unresolved if not str(item.get("field", "")).startswith("threshold")]


def _attach_legacy_compatibility_fields(result: dict[str, Any], econ_config: dict[str, Any], costs: dict[str, Any]) -> None:
    annual = result["annualByArm"]
    primary = annual[annual["discountRate"] == PRIMARY_DISCOUNT_RATE] if isinstance(annual, pd.DataFrame) and not annual.empty else pd.DataFrame()
    comp = primary[primary["arm"] == "comparator"] if not primary.empty else pd.DataFrame()
    inter = primary[primary["arm"] == "intervention"] if not primary.empty else pd.DataFrame()

    def mean_sum(frame: pd.DataFrame, column: str) -> float | None:
        if frame.empty or column not in frame:
            return None
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
        "baselineTBDiseaseCost": mean_sum(comp, "activeTBDiseaseCost"),
        "interventionTBDiseaseCost": mean_sum(inter, "activeTBDiseaseCost"),
    }
    cost_values["tbDiseaseCostsAverted"] = _subtract_if_number(
        cost_values["baselineTBDiseaseCost"],
        cost_values["interventionTBDiseaseCost"],
    )
    cost_values["totalProgramCost"] = _sum_optional(
        [
            cost_values["testingCost"],
            cost_values["treatmentCost"],
            cost_values["falsePositiveIncrementalCost"],
            cost_values["programSetupCost"],
            cost_values["programRunningCost"],
            mean_sum(inter, "adrManagementCost"),
        ]
    )
    cost_values["netCostVsBaseline"] = _subtract_if_number(
        cost_values["totalProgramCost"],
        cost_values["tbDiseaseCostsAverted"],
    )

    legacy_status = {
        "missingInputs": [item["field"] for item in result["unresolvedInputs"]],
        "notCalculated": [],
        "messages": [],
        "partialCalculations": [],
    }
    if costs.get("false_positive") is None:
        legacy_status["messages"].append("No false-positive incremental unit cost supplied; no extra false-positive cost is included.")
    for field, value in [
        ("testingCost", cost_values["testingCost"]),
        ("treatmentCost", cost_values["treatmentCost"]),
        ("tbDiseaseCostsAverted", cost_values["tbDiseaseCostsAverted"]),
        ("totalProgramCost", cost_values["totalProgramCost"]),
        ("netCostVsBaseline", cost_values["netCostVsBaseline"]),
    ]:
        if value is None:
            legacy_status["notCalculated"].append(field)

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
