from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any

import altair as alt
import pandas as pd
import streamlit as st

from app.display import (
    arrow_safe_dataframe,
    economics_assumptions_json,
    economics_summary_csv,
    safe_download_stem,
)
from app.health_economics_inputs import (
    assumptions_csv,
    assumptions_workbook,
    apply_assumptions_to_economics_config,
    assess_current_analysis_economic_readiness,
    conversion_audit_rows,
    fatal_validation_rows,
    mark_workspace_applied,
    mark_workspace_validated,
    new_workspace_state,
    parse_assumptions_csv,
    reconcile_workspace_state,
    rows_from_display_rows,
    update_workspace_rows,
    validate_editable_assumptions,
)
from app.health_economics_presentation import classify_icer_result
from app.state import (
    get_backend,
    has_retired_infection_history_results,
    init_session_state,
    mark_economics_changed,
    mark_economics_completed,
    record_message,
    sanitize_reference_only_state,
    sync_backend_status,
)
from engine.apy.evidence import assess_apy_reference_readiness, load_apy_evidence_registry
from engine.apy.event_ledger_economics import run_event_ledger_health_economics
from engine.apy.working_defaults import build_unified_working_default_preset


init_session_state()
st.session_state["apy_backend_name"] = "python_apy"
sanitize_reference_only_state()
backend = get_backend()


ASSUMPTION_SECTIONS = [
    "Screening and diagnostic costs",
    "Clinical pathway and preventive-treatment costs",
    "Active-TB care costs",
    "Programme and implementation costs",
    "Discounting and health outcomes",
    "Decision thresholds",
]

PROGRAMME_UNCOSTED_IDS = {
    "cost.program_setup",
    "cost.program_running",
    "cost.travel_outreach_staff_support",
}

DELIVERY_SCENARIO_DEFAULTS = {
    "includeAdditionalProgramCosts": False,
    "standaloneSetupCost": 0.0,
    "illustrativeSetupCost": 500000.0,
    "standaloneAnnualRunningCost": 0.0,
    "standaloneRunningYears": 2,
    "annualCostFirstYear": 0,
    "standaloneTravelOutreachCost": 0.0,
    "standaloneStaffSupportCost": 0.0,
    "sharedAttributionShare": 0.0,
}

POINT_COLLAPSE_TOLERANCE = {"cost": 0.01, "dalys": 1e-6}
REFERENCE_ECONOMICS_PATH = (
    Path(__file__).resolve().parents[1]
    / "app"
    / "reference_data"
    / "sa_health_report_reference_economics.json"
)

COST_COMPONENTS = [
    ("screeningTestCost", "Screening test"),
    ("returnForResultsCost", "Return/results"),
    ("clinicalReviewCost", "Clinical review"),
    ("activeTBExclusionWorkupCost", "Active-TB exclusion/work-up"),
    ("tptRegimenCost", "Preventive regimen"),
    ("adrManagementCost", "ADR management"),
    ("programSetupCost", "Programme setup"),
    ("programRunningCost", "Programme running"),
    ("travelOutreachStaffSupportCost", "Travel/outreach/staff support"),
    ("activeTBDiseaseCost", "Active-TB care"),
]

DELIVERY_COMPONENTS = {
    "screeningTestCost",
    "returnForResultsCost",
    "clinicalReviewCost",
    "activeTBExclusionWorkupCost",
    "tptRegimenCost",
    "adrManagementCost",
    "programSetupCost",
    "programRunningCost",
    "travelOutreachStaffSupportCost",
}

ECONOMIC_COST_WIDGET_SECTIONS = {
    "Screening and pathway costs": [
        "cost.test_igra",
        "cost.return_for_results",
        "cost.clinical_review",
        "cost.active_tb_exclusion_workup",
    ],
    "Preventive treatment": [
        "cost.regimen_3hp",
        "cost.tpt_adr_management",
    ],
    "Active-TB care": [
        "cost.active_tb_disease",
    ],
}

ECONOMIC_COST_STATE_LABELS = [
    "Included",
    "Excluded",
    "Bundled into another item",
    "Not locally costed",
]

ECONOMIC_COST_STATE_TO_INCLUSION = {
    "Included": "included",
    "Excluded": "excluded",
    "Bundled into another item": "bundled",
    "Not locally costed": "included",
}


def _money(value: Any, *, saving: bool = False) -> str:
    number = _number(value)
    if number is None:
        return "Unavailable"
    if saving:
        return f"AUD {abs(number):,.0f} saving"
    if number < 0:
        return f"AUD {abs(number):,.0f} saving"
    return f"AUD {number:,.0f}"


def _signed_money(value: Any) -> str:
    number = _number(value)
    if number is None:
        return "Unavailable"
    sign = "-" if number < 0 else ""
    return f"{sign}AUD {abs(number):,.0f}"


def _decimal(value: Any, digits: int = 1) -> str:
    number = _number(value)
    if number is None:
        return "Unavailable"
    return f"{number:,.{digits}f}"


def _number(value: Any) -> float | None:
    if value in (None, "", []):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _first_number(*values: Any, default: float = 0.0) -> float:
    for value in values:
        number = _number(value)
        if number is not None:
            return number
    return default


def _cost_state_label(row: dict[str, Any]) -> str:
    if row.get("reviewStatus") == "unresolved" and row.get("currentValue") in (None, "", []):
        return "Not locally costed"
    label = str(row.get("inclusionStatusLabel") or "").strip()
    if label in ECONOMIC_COST_STATE_LABELS:
        return label
    inclusion = str(row.get("inclusionStatus") or "included").strip()
    if inclusion == "excluded":
        return "Excluded"
    if inclusion == "bundled":
        return "Bundled into another item"
    return "Included"


def _bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def _selected_regimen(config: dict[str, Any] | None) -> str:
    cfg = config or {}
    return str(
        cfg.get("regimen")
        or cfg.get("preventiveRegimen")
        or (cfg.get("treatment") or {}).get("regimen")
        or "3HP"
    )


def _summary_row(econ_results: dict[str, Any] | None, metric: str, profile: str = "primary") -> dict[str, Any] | None:
    for row in (econ_results or {}).get("summaryRows") or []:
        if row.get("metric") == metric and row.get("discountProfile") == profile:
            return row
    return None


def _summary_mean(econ_results: dict[str, Any] | None, metric: str, profile: str = "primary") -> float | None:
    row = _summary_row(econ_results, metric, profile)
    return None if row is None else _number(row.get("mean"))


def icer_classification(incremental_cost: Any, dalys_averted: Any) -> dict[str, Any]:
    result = classify_icer_result(incremental_cost, dalys_averted)
    label = result["classification"]
    if label == "Dominant":
        display = "Dominant - better health and lower cost"
    elif label == "Higher cost with health gain":
        display = "Higher cost with health gain"
    elif label == "Dominated":
        display = "Dominated - higher cost without health gain"
    elif label == "Trade-off":
        display = "Trade-off - lower cost with lower health gain"
    else:
        display = "ICER not calculable"
    return {**result, "display": display}


def _signed_icer(value: Any) -> str:
    number = _number(value)
    if number is None:
        return "Not calculable"
    sign = "-" if number < 0 else ""
    return f"{sign}AUD {abs(number):,.0f} per DALY averted"


def _classification(econ_results: dict[str, Any] | None) -> str:
    row = _summary_row(econ_results, "primaryICER_ratioOfMeans")
    label = str((row or {}).get("classification") or "").strip()
    if label == "dominant":
        return (
            "Better modelled health outcomes and lower included health-system costs; "
            "provisional because local programme costs remain unresolved."
        )
    return label or "Unavailable"


def primary_economic_values(
    econ_results: dict[str, Any] | None,
    results_bundle: dict[str, Any] | None = None,
) -> dict[str, Any]:
    cost_rows = cost_category_rows(econ_results)
    incremental = _summary_mean(econ_results, "incrementalCost")
    dalys = _summary_mean(econ_results, "dalysAverted")
    active_tb_prevented = _summary_mean(econ_results, "activeTBCasesPrevented")
    gross = _gross_delivery(cost_rows)
    active_offset = next(
        (abs(row["Incremental cost"]) for row in cost_rows if row["Category"] == "Active-TB care"),
        None,
    )
    return {
        "incrementalCost": incremental,
        "dalysAverted": dalys,
        "activeTBAverted": active_tb_prevented,
        "grossDeliveryExpenditure": gross,
        "activeTBCareOffset": active_offset,
        "peopleScreened": _screened(results_bundle),
        "classification": icer_classification(incremental, dalys),
    }


def headline_panel_rows(
    *,
    econ_results: dict[str, Any] | None,
    results_bundle: dict[str, Any] | None,
) -> list[dict[str, str]]:
    values = primary_economic_values(econ_results, results_bundle)
    classification = values["classification"]
    if classification["classification"] == "Dominant":
        main = (
            f"Dominant - saves approximately {_money(values['incrementalCost'], saving=True).replace(' saving', '')} "
            f"and averts approximately {_decimal(values['dalysAverted'], 1)} DALYs"
        )
    else:
        main = classification["display"]
    return [
        {"Result": "Economic result", "Value": main},
        {"Result": "Incremental health-system cost", "Value": _signed_money(values["incrementalCost"])},
        {"Result": "DALYs averted", "Value": _decimal(values["dalysAverted"], 1)},
        {"Result": "Active TB cases averted", "Value": _decimal(values["activeTBAverted"], 1)},
        {"Result": "ICER / classification", "Value": classification["display"]},
        {"Result": "NMB and probability cost-effective", "Value": "Unavailable until a reviewed willingness-to-pay threshold is supplied"},
    ]


def _annual_primary(econ_results: dict[str, Any] | None) -> pd.DataFrame:
    annual = (econ_results or {}).get("annualByArm")
    frame = annual.copy() if isinstance(annual, pd.DataFrame) else pd.DataFrame(annual or [])
    if frame.empty:
        return frame
    if "discountProfile" in frame.columns:
        frame = frame[frame["discountProfile"].astype(str).eq("primary")].copy()
    return frame


def _screened(results_bundle: dict[str, Any] | None) -> float | None:
    totals = (((results_bundle or {}).get("technical") or {}).get("eventLedger") or {}).get("replicateTotals")
    frame = totals.copy() if isinstance(totals, pd.DataFrame) else pd.DataFrame(totals or [])
    if frame.empty:
        return None
    mask = frame["arm"].astype(str).eq("intervention") & frame["eventName"].astype(str).eq("screened")
    if not mask.any():
        return None
    return float(pd.to_numeric(frame.loc[mask, "value"], errors="coerce").mean())


def _discounted_component(frame: pd.DataFrame, component: str) -> pd.Series:
    values = pd.to_numeric(frame.get(component), errors="coerce").fillna(0.0)
    factors = pd.to_numeric(frame.get("costDiscountFactor", 1.0), errors="coerce").fillna(1.0)
    return values * factors


def cost_category_rows(econ_results: dict[str, Any] | None) -> list[dict[str, Any]]:
    frame = _annual_primary(econ_results)
    if frame.empty:
        return []
    rows = []
    for component, label in COST_COMPONENTS:
        if component not in frame.columns:
            continue
        work = frame[["replicateId", "arm"]].copy()
        work["cost"] = _discounted_component(frame, component)
        totals = work.groupby(["replicateId", "arm"], dropna=False)["cost"].sum().reset_index()
        means = totals.groupby("arm")["cost"].mean()
        comparator = float(means.get("comparator", 0.0))
        intervention = float(means.get("intervention", 0.0))
        rows.append(
            {
                "Category": label,
                "Comparator cost": comparator,
                "Intervention cost": intervention,
                "Incremental cost": intervention - comparator,
            }
        )
    return rows


def budget_impact_rows(econ_results: dict[str, Any] | None) -> list[dict[str, Any]]:
    frame = _annual_primary(econ_results)
    if frame.empty:
        return []
    rows = []
    for year, year_frame in frame.groupby("modelYear", dropna=False):
        work = year_frame[["replicateId", "arm"]].copy()
        work["delivery"] = sum(
            _discounted_component(year_frame, component)
            for component in DELIVERY_COMPONENTS
            if component in year_frame.columns
        )
        work["active_tb_care"] = _discounted_component(year_frame, "activeTBDiseaseCost")
        work["total"] = pd.to_numeric(year_frame.get("totalDiscountedCost"), errors="coerce").fillna(0.0)
        totals = work.groupby(["replicateId", "arm"], dropna=False)[["delivery", "active_tb_care", "total"]].sum().reset_index()
        means = totals.groupby("arm")[["delivery", "active_tb_care", "total"]].mean()
        comparator_total = float(means.loc["comparator", "total"]) if "comparator" in means.index else 0.0
        intervention_total = float(means.loc["intervention", "total"]) if "intervention" in means.index else 0.0
        rows.append(
            {
                "Year": int(float(year)),
                "Intervention delivery expenditure": float(means.loc["intervention", "delivery"]) if "intervention" in means.index else 0.0,
                "Comparator active-TB care": float(means.loc["comparator", "active_tb_care"]) if "comparator" in means.index else 0.0,
                "Intervention active-TB care": float(means.loc["intervention", "active_tb_care"]) if "intervention" in means.index else 0.0,
                "Annual incremental cost": intervention_total - comparator_total,
            }
        )
    cumulative = 0.0
    for row in sorted(rows, key=lambda item: item["Year"]):
        cumulative += row["Annual incremental cost"]
        row["Cumulative incremental cost"] = cumulative
    return rows


def _gross_delivery(cost_rows: list[dict[str, Any]]) -> float | None:
    if not cost_rows:
        return None
    return sum(
        float(row["Incremental cost"])
        for row in cost_rows
        if row["Category"] != "Active-TB care"
    )


def _ensure_delivery_scenario_controls() -> dict[str, Any]:
    controls = st.session_state.get("health_econ_delivery_scenarios")
    if not isinstance(controls, dict):
        controls = dict(DELIVERY_SCENARIO_DEFAULTS)
        st.session_state["health_econ_delivery_scenarios"] = controls
    for key, value in DELIVERY_SCENARIO_DEFAULTS.items():
        controls.setdefault(key, value)
    return controls


def _additional_programme_costs_enabled(controls: dict[str, Any]) -> bool:
    return _bool(controls.get("includeAdditionalProgramCosts"))


def _programme_control_values(controls: dict[str, Any]) -> dict[str, float]:
    return {
        "setup_cost": float(controls.get("standaloneSetupCost") or 0.0),
        "annual_running_cost": float(controls.get("standaloneAnnualRunningCost") or 0.0),
        "running_years": float(controls.get("standaloneRunningYears") or 0.0),
        "travel_outreach_cost": float(controls.get("standaloneTravelOutreachCost") or 0.0),
        "staff_support_cost": float(controls.get("standaloneStaffSupportCost") or 0.0),
        "attribution_share": max(0.0, min(float(controls.get("sharedAttributionShare") or 0.0), 1.0)),
    }


def _raw_programme_cost_total(values: dict[str, float]) -> float:
    return (
        values["setup_cost"]
        + values["annual_running_cost"] * values["running_years"]
        + values["travel_outreach_cost"]
        + values["staff_support_cost"]
    )


def _clone_bundle_with_running_duration(results_bundle: dict[str, Any] | None, years: int | None) -> dict[str, Any] | None:
    if not isinstance(results_bundle, dict) or years is None:
        return results_bundle
    out = deepcopy(results_bundle)
    ledger = ((out.setdefault("technical", {})).setdefault("eventLedger", {}))
    metadata = ledger.setdefault("metadata", {})
    metadata["screeningWindowYears"] = int(years)
    metadata["screeningWindow"] = int(years)
    return out


def _set_cost_item(
    economics_config: dict[str, Any],
    item_id: str,
    value: float,
    *,
    cost_basis: str | None = None,
    source: str,
    notes: str,
) -> None:
    for item in economics_config.get("costItems") or []:
        if item.get("costItemId") != item_id:
            continue
        item["originalCost"] = float(value)
        item["originalCurrency"] = "AUD"
        item["originalPriceYear"] = "2019"
        item["targetCurrency"] = "AUD"
        item["targetPriceYear"] = "2019"
        item["sourceCitation"] = source
        item["notes"] = notes
        if cost_basis:
            item.setdefault("resourceUse", {})["costBasis"] = cost_basis
        break
    for row in economics_config.get("assumptionEvidenceRegistry") or []:
        if row.get("assumptionId") != f"cost.{item_id}":
            continue
        row["currentValue"] = float(value)
        row["originalCurrency"] = "AUD"
        row["originalPriceYear"] = "2019"
        row["targetCurrency"] = "AUD"
        row["targetPriceYear"] = "2019"
        row["sourceCitation"] = source
        row["reviewStatus"] = "unreviewed_repository_input"
        row["reviewStatusLabel"] = "Unreviewed repository input"
        row["provisional"] = True
        row["inclusionStatus"] = "included"
        row["inclusionStatusLabel"] = "Included"
        row["notes"] = notes
        row["unresolvedReason"] = ""
        if cost_basis:
            row["costBasis"] = cost_basis
        break


def _delivery_scenario_config(
    base_econ: dict[str, Any],
    *,
    setup_cost: float,
    annual_running_cost: float,
    running_years: int,
    travel_outreach_cost: float,
    staff_support_cost: float,
    attribution_share: float = 1.0,
    scenario_name: str,
) -> dict[str, Any]:
    econ = deepcopy(base_econ)
    share = max(0.0, min(float(attribution_share), 1.0))
    source = f"User-defined {scenario_name} scenario assumption"
    _set_cost_item(
        econ,
        "program_setup",
        float(setup_cost) * share,
        cost_basis="total_once_at_program_start",
        source=source,
        notes="One-off programme setup or bulk implementation cost assigned to LTBI screening for this scenario.",
    )
    _set_cost_item(
        econ,
        "program_running",
        float(annual_running_cost) * share,
        cost_basis="annual_during_screening_window",
        source=source,
        notes=f"Annual programme running cost applied for {int(running_years)} year(s) in this scenario.",
    )
    _set_cost_item(
        econ,
        "travel_outreach_staff_support",
        (float(travel_outreach_cost) + float(staff_support_cost)) * share,
        cost_basis="per_person_screened",
        source=source,
        notes=(
            "Per-person-screened travel/outreach and additional staff/support cost. "
            "This does not bundle or exclude clinical pathway costs."
        ),
    )
    econ.setdefault("metadata", {})["economicScenarioName"] = scenario_name
    econ["metadata"]["programRunningDurationYears"] = int(running_years)
    econ["metadata"]["sharedProgrammeAttributionShare"] = share
    return econ


def _run_delivery_scenario(
    results_bundle: dict[str, Any],
    economics_config: dict[str, Any],
    *,
    running_years: int,
) -> dict[str, Any]:
    bundle = _clone_bundle_with_running_duration(results_bundle, running_years)
    return run_event_ledger_health_economics(bundle, economics_config)


def _effective_config_for_programme_controls(
    economics_config: dict[str, Any],
    controls: dict[str, Any],
) -> dict[str, Any]:
    if not _additional_programme_costs_enabled(controls):
        return economics_config
    values = _programme_control_values(controls)
    if _raw_programme_cost_total(values) <= 0:
        return economics_config
    return _delivery_scenario_config(
        economics_config,
        setup_cost=values["setup_cost"],
        annual_running_cost=values["annual_running_cost"],
        running_years=int(values["running_years"]),
        travel_outreach_cost=values["travel_outreach_cost"],
        staff_support_cost=values["staff_support_cost"],
        attribution_share=1.0,
        scenario_name="Current assumptions with user-defined additional programme costs",
    )


def delivery_scenario_comparison_rows(
    *,
    results_bundle: dict[str, Any] | None,
    economics_config: dict[str, Any],
    controls: dict[str, Any],
) -> list[dict[str, Any]]:
    if not isinstance(results_bundle, dict):
        return []
    values = _programme_control_values(controls)
    running_years = int(values["running_years"])
    existing = _run_delivery_scenario(results_bundle, economics_config, running_years=running_years)
    rows = [
        _delivery_scenario_row(
            "No additional programme overhead entered",
            existing,
            setup_cost=0.0,
            annual_running_cost=0.0,
            running_years=running_years,
            attribution_share=0.0,
            attributed_setup_cost=0.0,
            attributed_total_note="No additional programme overhead entered",
            point_role="current",
        )
    ]
    if _additional_programme_costs_enabled(controls) and _raw_programme_cost_total(values) > 0:
        standalone_config = _delivery_scenario_config(
            economics_config,
            setup_cost=values["setup_cost"],
            annual_running_cost=values["annual_running_cost"],
            running_years=running_years,
            travel_outreach_cost=values["travel_outreach_cost"],
            staff_support_cost=values["staff_support_cost"],
            attribution_share=1.0,
            scenario_name="Standalone programme - user-defined additional costs",
        )
        standalone = _run_delivery_scenario(results_bundle, standalone_config, running_years=running_years)
        rows.append(
            _delivery_scenario_row(
                "Standalone programme - user-defined additional costs",
                standalone,
                setup_cost=values["setup_cost"],
                annual_running_cost=values["annual_running_cost"],
                running_years=running_years,
                attribution_share=1.0,
                attributed_setup_cost=values["setup_cost"],
                attributed_total_note=f"Total additional setup assigned to LTBI screening: {_money(values['setup_cost'])}",
                point_role="scenario",
            )
        )
        if values["attribution_share"] > 0:
            shared_config = _delivery_scenario_config(
                economics_config,
                setup_cost=values["setup_cost"],
                annual_running_cost=values["annual_running_cost"],
                running_years=running_years,
                travel_outreach_cost=values["travel_outreach_cost"],
                staff_support_cost=values["staff_support_cost"],
                attribution_share=values["attribution_share"],
                scenario_name="Shared delivery - user-defined attributable share",
            )
            shared = _run_delivery_scenario(results_bundle, shared_config, running_years=running_years)
            rows.append(
                _delivery_scenario_row(
                    "Shared delivery - user-defined attributable share",
                    shared,
                    setup_cost=values["setup_cost"],
                    annual_running_cost=values["annual_running_cost"],
                    running_years=running_years,
                    attribution_share=values["attribution_share"],
                    attributed_setup_cost=values["setup_cost"] * values["attribution_share"],
                    attributed_total_note=(
                        f"Common setup {_money(values['setup_cost'])}; "
                        f"attributed share {values['attribution_share'] * 100:.0f}%; "
                        f"LTBI setup {_money(values['setup_cost'] * values['attribution_share'])}"
                    ),
                    point_role="scenario",
                )
            )
    return rows


def _delivery_scenario_row(
    label: str,
    econ_results: dict[str, Any],
    *,
    setup_cost: float,
    annual_running_cost: float,
    running_years: int,
    attribution_share: float,
    attributed_setup_cost: float,
    attributed_total_note: str,
    point_role: str,
) -> dict[str, Any]:
    values = primary_economic_values(econ_results)
    classification = values["classification"]
    arithmetic_icer = _divide(values["incrementalCost"], values["dalysAverted"])
    return {
        "Scenario": label,
        "Additional setup cost": _money(attributed_setup_cost),
        "Annual programme cost": f"{_money(annual_running_cost * attribution_share)} for {int(running_years)} year(s)",
        "Attributed share": f"{attribution_share * 100:.0f}%",
        "Incremental cost": _signed_money(values["incrementalCost"]),
        "DALYs averted": _decimal(values["dalysAverted"], 1),
        "Quadrant": _quadrant(values["dalysAverted"], values["incrementalCost"]),
        "Classification": "Dominant" if classification["classification"] == "Dominant" else classification["display"],
        "Arithmetic ICER": _signed_icer(arithmetic_icer),
        "Active TB averted": _decimal(values["activeTBAverted"], 1),
        "Programme-cost assumption": attributed_total_note,
        "_incrementalCostRaw": values["incrementalCost"],
        "_dalysAvertedRaw": values["dalysAverted"],
        "_activeTBAvertedRaw": values["activeTBAverted"],
        "_pointRole": point_role,
    }


def break_even_rows(econ_results: dict[str, Any] | None) -> list[dict[str, str]]:
    values = primary_economic_values(econ_results)
    inc = _number(values["incrementalCost"])
    if inc is None or inc >= 0:
        return []
    return [
        {
            "Planning measure": "Additional programme cost before net saving is exhausted",
            "Value": _money(abs(inc)),
            "Interpretation": "Discounted cost that could be added before the intervention ceases to be cost-saving under these assumptions.",
        }
    ]


def cost_effectiveness_plane_rows(
    *,
    results_bundle: dict[str, Any] | None,
    economics_config: dict[str, Any],
    current_economics: dict[str, Any] | None,
    controls: dict[str, Any],
) -> list[dict[str, Any]]:
    rows = [_plane_origin_row()]
    reference = _reference_report_point()
    values = _programme_control_values(controls)
    current_setup = values["setup_cost"] if _additional_programme_costs_enabled(controls) else 0.0
    current_annual = (
        f"{_money(values['annual_running_cost'])} for {int(values['running_years'])} year(s)"
        if _additional_programme_costs_enabled(controls)
        else "AUD 0"
    )
    current_share = "100%" if _additional_programme_costs_enabled(controls) and _raw_programme_cost_total(values) > 0 else "0%"
    current = _plane_row(
        "Current assumptions",
        current_economics,
        role="current",
        event_ledger_source="Current completed analysis",
        analysis_basis=_analysis_basis_label(results_bundle),
        setup_cost=_money(current_setup),
        annual_cost=current_annual,
        attribution_share=current_share,
        represents_frozen_reference=False,
    ) if current_economics else None
    if reference and current and _points_overlap(reference, current):
        collapsed = dict(current)
        collapsed["Scenario"] = "Current assumptions - same as SA Health report reference"
        collapsed["Point note"] = "Current analysis reproduces the frozen report-reference coordinates within tolerance."
        collapsed["representsFrozenReference"] = True
        rows.append(collapsed)
    else:
        if reference:
            rows.append(reference)
        if current:
            rows.append(current)
    if isinstance(results_bundle, dict):
        for scenario in delivery_scenario_comparison_rows(
            results_bundle=results_bundle,
            economics_config=economics_config,
            controls=controls,
        ):
            if scenario.get("_pointRole") != "scenario":
                continue
            rows.append(
                _plane_row_from_values(
                    label=scenario["Scenario"],
                    incremental_cost=scenario["_incrementalCostRaw"],
                    dalys_averted=scenario["_dalysAvertedRaw"],
                    active_tb_averted=scenario.get("_activeTBAvertedRaw"),
                    role="scenario",
                    event_ledger_source="Current completed analysis",
                    analysis_basis=_analysis_basis_label(results_bundle),
                    setup_cost=scenario["Additional setup cost"],
                    annual_cost=scenario["Annual programme cost"],
                    attribution_share=scenario["Attributed share"],
                    represents_frozen_reference=False,
                    point_note=scenario["Programme-cost assumption"],
                )
            )
    return _collapse_duplicate_points(rows)


def _plane_origin_row() -> dict[str, Any]:
    return {
        "Scenario": "Business as usual",
        "DALYs averted compared with business as usual": 0.0,
        "Incremental cost compared with business as usual (AUD)": 0.0,
        "Quadrant": "Origin",
        "Classification": "Business as usual",
        "Arithmetic ICER": "Not applicable",
    }


def _plane_row(
    label: str,
    econ_results: dict[str, Any] | None,
    *,
    role: str,
    event_ledger_source: str,
    analysis_basis: str,
    setup_cost: Any,
    annual_cost: Any,
    attribution_share: Any,
    represents_frozen_reference: bool,
) -> dict[str, Any]:
    if not econ_results:
        return {}
    values = primary_economic_values(econ_results)
    return _plane_row_from_values(
        label=label,
        incremental_cost=values["incrementalCost"],
        dalys_averted=values["dalysAverted"],
        active_tb_averted=values["activeTBAverted"],
        role=role,
        event_ledger_source=event_ledger_source,
        analysis_basis=analysis_basis,
        setup_cost=setup_cost,
        annual_cost=annual_cost,
        attribution_share=attribution_share,
        represents_frozen_reference=represents_frozen_reference,
        point_note="",
    )


def _plane_row_from_values(
    *,
    label: str,
    incremental_cost: Any,
    dalys_averted: Any,
    active_tb_averted: Any,
    role: str,
    event_ledger_source: str,
    analysis_basis: str,
    setup_cost: Any,
    annual_cost: Any,
    attribution_share: Any,
    represents_frozen_reference: bool,
    point_note: str,
) -> dict[str, Any]:
    inc = _number(incremental_cost)
    dalys = _number(dalys_averted)
    classification = icer_classification(inc, dalys)
    return {
        "Scenario": label,
        "DALYs averted compared with business as usual": dalys,
        "Incremental cost compared with business as usual (AUD)": inc,
        "Quadrant": _quadrant(dalys, inc),
        "Classification": "Dominant" if classification["classification"] == "Dominant" else classification["display"],
        "Arithmetic ICER": _signed_icer(_divide(inc, dalys)),
        "Active TB averted": active_tb_averted,
        "Event ledger": event_ledger_source,
        "Analysis basis": analysis_basis,
        "Additional setup cost": setup_cost,
        "Annual programme cost": annual_cost,
        "Attributed share": attribution_share,
        "Point role": role,
        "Point note": point_note,
        "representsFrozenReference": represents_frozen_reference,
    }


def _reference_report_point() -> dict[str, Any] | None:
    try:
        payload = json.loads(REFERENCE_ECONOMICS_PATH.read_text(encoding="utf-8"))
    except OSError:
        return None
    return _plane_row_from_values(
        label="SA Health report reference",
        incremental_cost=payload.get("incrementalCost"),
        dalys_averted=payload.get("dalysAverted"),
        active_tb_averted=payload.get("activeTBCasesPrevented"),
        role="reference",
        event_ledger_source="Frozen stochastic compatibility reference package",
        analysis_basis=(
            f"MATLAB-v9-compatible stochastic anchor; {payload.get('nReps')} repetitions; "
            f"seed {payload.get('seed')}"
        ),
        setup_cost="AUD 0 additional overhead entered",
        annual_cost="AUD 0",
        attribution_share="0%",
        represents_frozen_reference=True,
        point_note=str(payload.get("notes") or ""),
    )


def _analysis_basis_label(results_bundle: dict[str, Any] | None) -> str:
    ledger = ((results_bundle or {}).get("technical") or {}).get("eventLedger") or {}
    metadata = ledger.get("metadata") or {}
    model_type = metadata.get("modelType") or (results_bundle or {}).get("metadata", {}).get("modelType")
    natural_history = metadata.get("naturalHistorySemantics") or metadata.get("ltbiStateModel") or ""
    reps = metadata.get("nReps") or metadata.get("replicates")
    seed = metadata.get("seed")
    parts = [str(item) for item in [model_type, natural_history] if item]
    if reps:
        parts.append(f"{reps} repetitions")
    if seed is not None:
        parts.append(f"seed {seed}")
    return "; ".join(parts) or "Current completed analysis"


def _points_overlap(left: dict[str, Any], right: dict[str, Any]) -> bool:
    left_cost = _number(left.get("Incremental cost compared with business as usual (AUD)"))
    right_cost = _number(right.get("Incremental cost compared with business as usual (AUD)"))
    left_dalys = _number(left.get("DALYs averted compared with business as usual"))
    right_dalys = _number(right.get("DALYs averted compared with business as usual"))
    if None in {left_cost, right_cost, left_dalys, right_dalys}:
        return False
    return (
        abs(left_cost - right_cost) <= POINT_COLLAPSE_TOLERANCE["cost"]
        and abs(left_dalys - right_dalys) <= POINT_COLLAPSE_TOLERANCE["dalys"]
    )


def _collapse_duplicate_points(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    collapsed: list[dict[str, Any]] = []
    for row in rows:
        if not row:
            continue
        match = next((item for item in collapsed if _points_overlap(item, row)), None)
        if match is None:
            collapsed.append(row)
            continue
        match["Scenario"] = f"{match['Scenario']} / {row['Scenario']}"
        match["Point note"] = "Multiple scenarios share these coordinates within plotting tolerance."
        match["Point role"] = f"{match.get('Point role')}, {row.get('Point role')}"
        match["representsFrozenReference"] = bool(match.get("representsFrozenReference") or row.get("representsFrozenReference"))
    return collapsed


def _quadrant(dalys_averted: Any, incremental_cost: Any) -> str:
    x = _number(dalys_averted)
    y = _number(incremental_cost)
    if x is None or y is None:
        return "Unavailable"
    if abs(x) < 1e-12 and abs(y) < 1e-12:
        return "Origin"
    if x > 0 and y < 0:
        return "Lower right"
    if x > 0 and y > 0:
        return "Upper right"
    if x < 0 and y > 0:
        return "Upper left"
    if x < 0 and y < 0:
        return "Lower left"
    return "On axis"


def _divide(a: Any, b: Any) -> float | None:
    a_num = _number(a)
    b_num = _number(b)
    if a_num is None or b_num is None or abs(b_num) < 1e-12:
        return None
    return a_num / b_num


def _cost_effectiveness_plane_chart(rows: list[dict[str, Any]]) -> alt.Chart:
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["Chart label"] = frame.apply(
            lambda row: str(row["Scenario"])
            if str(row.get("Point role", "")) in {"current", "reference"} or "same as SA Health" in str(row.get("Scenario", ""))
            else "",
            axis=1,
        )
    x_col = "DALYs averted compared with business as usual"
    y_col = "Incremental cost compared with business as usual (AUD)"
    points = (
        alt.Chart(frame)
        .mark_point(filled=True, size=120)
        .encode(
            x=alt.X(f"{x_col}:Q", title=x_col),
            y=alt.Y(f"{y_col}:Q", title=y_col),
            color=alt.Color("Scenario:N", legend=alt.Legend(title="Scenario")),
            tooltip=[
                alt.Tooltip("Scenario:N"),
                alt.Tooltip(f"{x_col}:Q", format=".2f"),
                alt.Tooltip(f"{y_col}:Q", format=",.0f"),
                alt.Tooltip("Quadrant:N"),
                alt.Tooltip("Classification:N"),
                alt.Tooltip("Arithmetic ICER:N"),
                alt.Tooltip("Event ledger:N"),
                alt.Tooltip("Analysis basis:N"),
                alt.Tooltip("Point note:N"),
            ],
        )
    )
    labels = (
        alt.Chart(frame)
        .mark_text(align="left", dx=8, dy=-6, fontSize=12)
        .encode(
            x=f"{x_col}:Q",
            y=f"{y_col}:Q",
            text="Chart label:N",
        )
    )
    horizontal = alt.Chart(pd.DataFrame({y_col: [0]})).mark_rule(color="#6b7280").encode(y=f"{y_col}:Q")
    vertical = alt.Chart(pd.DataFrame({x_col: [0]})).mark_rule(color="#6b7280").encode(x=f"{x_col}:Q")
    return (horizontal + vertical + points + labels).properties(height=360)


def headline_rows(
    *,
    econ_results: dict[str, Any] | None,
    results_bundle: dict[str, Any] | None,
) -> list[dict[str, str]]:
    cost_rows = cost_category_rows(econ_results)
    gross = _gross_delivery(cost_rows)
    incremental = _summary_mean(econ_results, "incrementalCost")
    active_tb_prevented = _summary_mean(econ_results, "activeTBCasesPrevented")
    dalys = _summary_mean(econ_results, "dalysAverted")
    screened = _screened(results_bundle)
    active_offset = next(
        (abs(row["Incremental cost"]) for row in cost_rows if row["Category"] == "Active-TB care"),
        None,
    )
    gross_per_screened = None if gross is None or not screened else gross / screened
    net_per_case = None if incremental is None or not active_tb_prevented else incremental / active_tb_prevented
    threshold = ((st.session_state.get("economics_config") or {}).get("threshold") or {}).get("value")
    return [
        {"Result": "Gross delivery expenditure", "Value": _money(gross)},
        {"Result": "Active-TB care cost offset", "Value": _money(active_offset, saving=True)},
        {"Result": "Net incremental health-system result", "Value": _money(incremental)},
        {"Result": "Gross delivery expenditure per person screened", "Value": _money(gross_per_screened)},
        {"Result": "Net saving or cost per active TB case averted", "Value": _money(net_per_case)},
        {"Result": "Provisional DALYs averted", "Value": _decimal(dalys, 4)},
        {"Result": "Economic classification", "Value": _classification(econ_results)},
        {
            "Result": "Net monetary benefit",
            "Value": "Unavailable until a reviewed willingness-to-pay threshold is supplied"
            if threshold in (None, "", [])
            else "Available in detailed outputs",
        },
    ]


def _source_label(row: dict[str, Any]) -> str:
    source = str(row.get("sourceCitation") or "").strip()
    if source == "User-defined":
        return "User-defined"
    if not source:
        return "Unresolved or not locally costed"
    if len(source) > 90:
        return source[:87].rstrip() + "..."
    return source


def _unit_label(row: dict[str, Any]) -> str:
    unit = str(row.get("unit") or "").strip()
    basis = str(row.get("costBasis") or "").strip()
    if unit and basis:
        return f"{unit}; {basis}"
    return unit or basis


def standard_assumption_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "assumptionId": row.get("assumptionId", ""),
            "Parameter": row.get("description") or row.get("assumptionId", ""),
            "Value used by model": row.get("currentValue", ""),
            "Unit": _unit_label(row),
            "Source": _source_label(row),
        }
        for row in rows
    ]


def _assumption_section(rows: list[dict[str, Any]], section: str) -> list[dict[str, Any]]:
    if section == "Screening and diagnostic costs":
        return [row for row in rows if row.get("assumptionId") in {"cost.test_igra", "cost.test_tst"}]
    if section == "Clinical pathway and preventive-treatment costs":
        return [
            row for row in rows
            if row.get("assumptionId") in {
                "cost.regimen_3hp",
                "cost.regimen_4r",
                "cost.regimen_3hr",
                "cost.regimen_6h",
                "cost.regimen_9h",
                "cost.tpt_adr_management",
                "cost.false_positive_incremental",
                "cost.return_for_results",
                "cost.clinical_review",
                "cost.active_tb_exclusion_workup",
            }
        ]
    if section == "Active-TB care costs":
        return [row for row in rows if row.get("assumptionId") == "cost.active_tb_disease"]
    if section == "Programme and implementation costs":
        return [row for row in rows if row.get("assumptionId") in PROGRAMME_UNCOSTED_IDS]
    if section == "Discounting and health outcomes":
        return [row for row in rows if row.get("category") == "daly"]
    if section == "Decision thresholds":
        return [row for row in rows if row.get("category") == "threshold"]
    return rows


def _merge_standard_edits(
    *,
    base_rows: list[dict[str, Any]],
    edited_records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    edited_by_id = {row.get("assumptionId"): row for row in edited_records}
    out = deepcopy(base_rows)
    for row in out:
        edited = edited_by_id.get(row.get("assumptionId"))
        if edited is None:
            continue
        new_value = edited.get("Value used by model")
        old_value = row.get("currentValue", "")
        if str(new_value) != str(old_value):
            row["currentValue"] = new_value
            row["sourceCitation"] = "User-defined"
    return out


def render_standard_assumption_editor(
    *,
    rows: list[dict[str, Any]],
    working_rows: list[dict[str, Any]],
    editor_prefix: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    latest_rows = list(working_rows)
    tabs = st.tabs(ASSUMPTION_SECTIONS)
    for tab, section in zip(tabs, ASSUMPTION_SECTIONS):
        with tab:
            section_rows = _assumption_section(rows, section)
            if not section_rows:
                st.info("No assumptions in this section.")
                continue
            edited = st.data_editor(
                arrow_safe_dataframe(standard_assumption_rows(section_rows)),
                use_container_width=True,
                hide_index=True,
                num_rows="fixed",
                key=f"{editor_prefix}_{section}",
                disabled=["Parameter", "Unit", "Source"],
                column_config={"assumptionId": None},
            )
            records = edited.to_dict(orient="records") if hasattr(edited, "to_dict") else list(edited)
            latest_rows = _merge_standard_edits(base_rows=latest_rows, edited_records=records)
    updated_state = update_workspace_rows(
        st.session_state["health_econ_workspace"],
        rows_from_display_rows(latest_rows),
    )
    return updated_state["rows"], updated_state


def render_editable_cost_workspace(
    *,
    rows: list[dict[str, Any]],
    editor_prefix: str,
) -> list[dict[str, Any]]:
    updated = deepcopy(rows)
    row_by_id = {row.get("assumptionId"): row for row in updated}
    version = int(st.session_state.get("health_econ_cost_widget_version", 0))
    tabs = st.tabs(list(ECONOMIC_COST_WIDGET_SECTIONS))
    for tab, (section, assumption_ids) in zip(tabs, ECONOMIC_COST_WIDGET_SECTIONS.items()):
        with tab:
            for assumption_id in assumption_ids:
                row = row_by_id.get(assumption_id)
                if not row:
                    continue
                label = str(row.get("description") or assumption_id)
                source = str(row.get("sourceCitation") or "")
                unit = str(row.get("unit") or "AUD")
                current_state = _cost_state_label(row)
                state_key = f"{editor_prefix}_{version}_{assumption_id}_state"
                value_key = f"{editor_prefix}_{version}_{assumption_id}_value"
                state = st.selectbox(
                    f"{label} - cost state",
                    ECONOMIC_COST_STATE_LABELS,
                    index=ECONOMIC_COST_STATE_LABELS.index(current_state)
                    if current_state in ECONOMIC_COST_STATE_LABELS
                    else 0,
                    key=state_key,
                    help=(
                        "Included values enter the calculation. Bundled values are not separately costed. "
                        "Not locally costed means evidence is still needed and no additional cost is applied until a value is entered."
                    ),
                )
                disabled_value = state in {"Excluded", "Bundled into another item", "Not locally costed"}
                value = st.number_input(
                    f"{label} - Value used by model",
                    min_value=0.0,
                    value=_first_number(row.get("currentValue"), default=0.0),
                    step=10.0,
                    format="%.4f",
                    key=value_key,
                    disabled=disabled_value,
                    help=f"{unit}. Source: {source or 'Source not supplied'}.",
                )
                st.caption(f"Unit: {unit} | Source: {row.get('sourceCitation') or 'Source not supplied'}")
                old_state = _cost_state_label(row)
                old_value = _first_number(row.get("currentValue"), default=0.0)
                changed = state != old_state or (not disabled_value and abs(float(value) - old_value) > 1e-9)
                row["inclusionStatus"] = ECONOMIC_COST_STATE_TO_INCLUSION[state]
                row["inclusionStatusLabel"] = "Included" if state == "Not locally costed" else state
                if state == "Included":
                    row["currentValue"] = float(value)
                elif state in {"Excluded", "Bundled into another item"}:
                    row["currentValue"] = ""
                else:
                    row["currentValue"] = ""
                    row["reviewStatus"] = "unresolved"
                    row["reviewStatusLabel"] = "Unresolved"
                if changed:
                    row["sourceCitation"] = "User-defined"
                    if state == "Included":
                        row["reviewStatus"] = row.get("reviewStatus") or "configured_reviewed"
                        row["reviewStatusLabel"] = row.get("reviewStatusLabel") or "Reviewed numerical assumption"
                    st.caption("User-defined")
    return updated


def overridden_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "Parameter": row.get("description") or row.get("assumptionId", ""),
            "Value used by model": row.get("currentValue", ""),
            "Source": "User-defined",
        }
        for row in rows
        if str(row.get("sourceCitation") or "").strip() == "User-defined"
    ]


def cost_state_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected = []
    for row in rows:
        if row.get("category") != "cost":
            continue
        selected.append(
            {
                "Parameter": row.get("description") or row.get("assumptionId"),
                "Value used by model": row.get("currentValue", ""),
                "State": row.get("inclusionStatusLabel") or row.get("inclusionStatus"),
                "Interpretation": row.get("notes") or row.get("unresolvedReason") or "",
            }
        )
    return selected


def load_economics_config(new_config: dict[str, Any], *, mark_workspace_unsaved: bool = False) -> None:
    st.session_state["economics_config"] = new_config
    st.session_state["economics_results"] = None
    st.session_state["dirty_economics"] = False
    workspace = reconcile_workspace_state(
        st.session_state.get("health_econ_workspace"),
        new_config,
        registry=new_config.get("assumptionEvidenceRegistry"),
    )
    if mark_workspace_unsaved:
        workspace["hasUnsavedEdits"] = True
        workspace["validated"] = False
        workspace["validation"] = None
        workspace["applied"] = False
        workspace["baselineRowsHash"] = ""
    st.session_state["health_econ_workspace"] = workspace
    sync_backend_status(backend.status())


def run_authoritative_health_economics(results_bundle: dict, econ_config: dict) -> None:
    try:
        controls = _ensure_delivery_scenario_controls()
        running_years = int(float(controls.get("standaloneRunningYears") or 0))
        econ = run_event_ledger_health_economics(
            _clone_bundle_with_running_duration(results_bundle, running_years),
            econ_config,
        )
        st.session_state["economics_results"] = econ
        mark_economics_completed()
        sync_backend_status(backend.status())
        st.success("Health-economic analysis completed from the current screening outcomes.")
    except Exception as exc:
        message = f"Health-economic analysis failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)


def apply_and_recalculate(
    *,
    working_rows: list[dict[str, Any]],
    workspace_state: dict[str, Any],
    econ_config: dict[str, Any],
    config: dict[str, Any] | None,
    results_bundle: dict[str, Any] | None,
    controls: dict[str, Any],
) -> None:
    validation = validate_editable_assumptions(working_rows, econ_config, config=config or {})
    st.session_state["health_econ_workspace"] = mark_workspace_validated(workspace_state, validation)
    if not validation.get("isValidForApplication"):
        st.error("Economic assumptions contain errors that must be corrected before recalculation.")
        fatal_rows = fatal_validation_rows(validation)
        if fatal_rows:
            st.dataframe(arrow_safe_dataframe(fatal_rows), use_container_width=True, hide_index=True)
        return
    updated_config = apply_assumptions_to_economics_config(
        econ_config,
        working_rows,
        config=config or {},
    )
    applied_state = new_workspace_state(
        updated_config,
        updated_config.get("assumptionEvidenceRegistry") or working_rows,
    )
    applied_state = mark_workspace_validated(
        applied_state,
        validate_editable_assumptions(applied_state["rows"], updated_config, config=config or {}),
    )
    applied_state = mark_workspace_applied(applied_state, updated_config)
    st.session_state["economics_config"] = updated_config
    st.session_state["health_econ_workspace"] = applied_state
    mark_economics_changed()
    if results_bundle and not st.session_state.get("results_stale"):
        run_authoritative_health_economics(
            results_bundle,
            _effective_config_for_programme_controls(updated_config, controls),
        )
        st.success("Economic results recalculated from the current screening outcomes. Epidemiological results were not rerun.")
    else:
        st.warning("Run the screening analysis before recalculating health economics.")


st.title("Health Economics")
st.caption("Review economic results, optionally change cost assumptions, and export the analysis.")

config = st.session_state.get("config")
results_bundle = st.session_state.get("results_bundle")
econ_config = st.session_state.get("economics_config")
econ_results = st.session_state.get("economics_results")
scenario_label = (results_bundle or {}).get("metadata", {}).get("scenarioLabel")
can_run = bool(config and results_bundle and not st.session_state.get("results_stale"))

if not econ_config:
    econ_config = build_unified_working_default_preset()["economicsConfig"]
    st.session_state["economics_config"] = econ_config

ledger = (results_bundle or {}).get("technical", {}).get("eventLedger", {}) if isinstance(results_bundle, dict) else {}
retired_ledger = has_retired_infection_history_results(results_bundle)
workspace_state = reconcile_workspace_state(
    st.session_state.get("health_econ_workspace"),
    econ_config,
    registry=econ_config.get("assumptionEvidenceRegistry"),
)
st.session_state["health_econ_workspace"] = workspace_state
working_rows = workspace_state["rows"]
override_count = len(overridden_rows(working_rows))
controls = _ensure_delivery_scenario_controls()

if retired_ledger:
    st.warning("Previous results used an analysis pathway that is not available in this SA Health version. Please run the analysis again.")
    st.stop()

if not can_run:
    notice = st.session_state.pop("reference_only_migration_notice", "")
    if notice:
        st.warning(notice)
    st.info("Run the screening analysis first; then return here to calculate health-economic results.")
    st.stop()

if econ_results is None:
    run_authoritative_health_economics(
        results_bundle,
        _effective_config_for_programme_controls(econ_config, controls),
    )
    econ_results = st.session_state.get("economics_results")

metadata = econ_config.get("metadata") or {}
currency = metadata.get("targetCurrency") or metadata.get("currencyCode") or "AUD"
price_year = metadata.get("targetPriceYear") or metadata.get("priceYear") or "2019"
discounting = econ_config.get("discounting") or {}
primary_rate = (
    discounting.get("primaryDisplayedRate")
    or discounting.get("selectedAnnualRate")
    or ((discounting.get("profiles") or {}).get("primary") or {}).get("costRate")
    or 0.03
)
assumption_status = "SA Health working defaults" if override_count == 0 else f"{override_count} user-defined override(s)"
analysis_basis_status = _analysis_basis_label(results_bundle)
st.caption(
    f"Analysis basis: {analysis_basis_status} | Perspective: {metadata.get('perspective', 'Australian health-care system')} | "
    f"{currency} {price_year} | Primary discount rate: {float(primary_rate) * 100:.1f}% | Economic assumptions: {assumption_status}"
)
if any(row.get("assumptionId") in PROGRAMME_UNCOSTED_IDS and row.get("currentValue") in (0, 0.0, "0", "0.0") for row in working_rows):
    st.warning("Programme setup, running, travel, outreach and staff-support costs have not yet been locally costed.")

st.subheader("Economic result")
if econ_results:
    st.dataframe(
        arrow_safe_dataframe(headline_panel_rows(econ_results=econ_results, results_bundle=results_bundle)),
        use_container_width=True,
        hide_index=True,
    )
    st.caption("DALY-based results are provisional pending evidence review.")

st.subheader("Programme-delivery scenarios")
st.caption(
    "All scenarios reuse the same screening outcomes. No additional programme overhead entered means no local overhead has been "
    "entered; it does not mean implementation is costless. AUD 500,000, when selected below, is illustrative and user-changeable."
)
try:
    scenario_rows = delivery_scenario_comparison_rows(
        results_bundle=results_bundle,
        economics_config=econ_config,
        controls=controls,
    )
    st.dataframe(
        arrow_safe_dataframe([{key: value for key, value in row.items() if not key.startswith("_")} for row in scenario_rows]),
        use_container_width=True,
        hide_index=True,
    )
    st.markdown("Incremental cost-effectiveness plane")
    plane_rows = cost_effectiveness_plane_rows(
        results_bundle=results_bundle,
        economics_config=econ_config,
        current_economics=econ_results,
        controls=controls,
    )
    st.altair_chart(_cost_effectiveness_plane_chart(plane_rows), use_container_width=True)
    st.dataframe(
        arrow_safe_dataframe(
            [
                {
                    "Scenario": row.get("Scenario"),
                    "DALYs averted": _decimal(row.get("DALYs averted compared with business as usual"), 4),
                    "Incremental cost": _signed_money(row.get("Incremental cost compared with business as usual (AUD)")),
                    "Event ledger": row.get("Event ledger", ""),
                    "Analysis basis": row.get("Analysis basis", ""),
                    "Additional setup cost": row.get("Additional setup cost", ""),
                    "Annual programme cost": row.get("Annual programme cost", ""),
                    "Attributed share": row.get("Attributed share", ""),
                    "Point note": row.get("Point note", ""),
                }
                for row in plane_rows
                if row.get("Scenario") != "Business as usual"
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )
    st.caption(
        "Cost-only changes move points vertically on this plane: added programme costs move upward, "
        "greater savings move downward, and DALYs averted remain fixed."
    )
except Exception as exc:
    st.error(f"Economic scenario comparison failed: {exc}")

st.subheader("Cost breakdown and budget impact")
if econ_results:
    cost_rows = cost_category_rows(econ_results)
    values = primary_economic_values(econ_results, results_bundle)
    st.dataframe(
        arrow_safe_dataframe(
            [
                {"Measure": "Gross intervention delivery cost", "Value": _money(values["grossDeliveryExpenditure"])},
                {"Measure": "Active-TB care cost offset", "Value": _money(values["activeTBCareOffset"], saving=True)},
                {"Measure": "Net incremental health-system cost", "Value": _signed_money(values["incrementalCost"])},
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )
    if cost_rows:
        st.dataframe(arrow_safe_dataframe(cost_rows), use_container_width=True, hide_index=True)
    budget_rows = budget_impact_rows(econ_results)
    if budget_rows:
        st.line_chart(
            pd.DataFrame(budget_rows).set_index("Year")[
                ["Intervention delivery expenditure", "Comparator active-TB care", "Intervention active-TB care"]
            ],
            use_container_width=True,
        )
        st.dataframe(arrow_safe_dataframe(budget_rows), use_container_width=True, hide_index=True)
    be_rows = break_even_rows(econ_results)
    if be_rows:
        st.dataframe(arrow_safe_dataframe(be_rows), use_container_width=True, hide_index=True)
    st.caption(
        "Programme and pathway expenditure is concentrated early; active-TB care offsets occur over follow-up. "
        "The break-even amount is a planning measure, not a willingness-to-pay threshold."
    )

with st.expander("Change cost assumptions", expanded=False):
    st.caption(
        "Economic changes reuse the current screening outcomes; epidemiology is not rerun. Blank entries are not converted to zero."
    )
    st.info(
        "Cost states: numerical costs are included directly; bundled/absorbed means the activity occurs but is not costed separately; "
        "reviewed exclusion means outside the selected perspective; not locally costed means evidence is still needed; compatibility zero is a placeholder."
    )
    st.markdown("Additional programme costs")
    controls["includeAdditionalProgramCosts"] = st.toggle(
        "Include additional programme costs",
        value=_additional_programme_costs_enabled(controls),
        help="When off, no additional programme overhead is applied. This is not evidence that implementation is costless.",
    )
    st.caption(
        "AUD 500,000 is available as an illustrative user-changeable scenario assumption, not an SA Health cost estimate."
    )
    if controls["includeAdditionalProgramCosts"]:
        if st.button("Use illustrative AUD 500,000 setup scenario"):
            controls["standaloneSetupCost"] = float(controls.get("illustrativeSetupCost") or 500000.0)
            controls["sharedAttributionShare"] = 0.50
            st.rerun()
        scenario_cols = st.columns(3)
        controls["standaloneSetupCost"] = scenario_cols[0].number_input(
            "One-off setup/bulk implementation cost (AUD)",
            min_value=0.0,
            value=float(controls.get("standaloneSetupCost") or 0.0),
            step=50000.0,
            format="%.0f",
            help="AUD 2019 total one-off cost. Enter 500000 only if you want the illustrative setup scenario.",
        )
        controls["standaloneAnnualRunningCost"] = scenario_cols[1].number_input(
            "Annual programme running cost (AUD per year)",
            min_value=0.0,
            value=float(controls.get("standaloneAnnualRunningCost") or 0.0),
            step=10000.0,
            format="%.0f",
            help="AUD 2019 annual overhead applied during the selected programme years.",
        )
        controls["standaloneRunningYears"] = int(
            scenario_cols[2].number_input(
                "Number of years annual cost applies",
                min_value=0,
                max_value=20,
                value=int(float(controls.get("standaloneRunningYears") or 0)),
                step=1,
            )
        )
        scenario_cols = st.columns(3)
        controls["annualCostFirstYear"] = int(
            scenario_cols[0].number_input(
                "First annual-cost year",
                min_value=0,
                max_value=20,
                value=int(float(controls.get("annualCostFirstYear") or 0)),
                step=1,
                disabled=True,
                help="The current economics engine applies annual programme running costs from programme start during the screening window.",
            )
        )
        controls["standaloneTravelOutreachCost"] = scenario_cols[1].number_input(
            "Travel/outreach cost (AUD per person screened)",
            min_value=0.0,
            value=float(controls.get("standaloneTravelOutreachCost") or 0.0),
            step=10.0,
            format="%.0f",
            help="AUD 2019 per person screened.",
        )
        controls["standaloneStaffSupportCost"] = scenario_cols[2].number_input(
            "Additional staff/support cost (AUD per person screened)",
            min_value=0.0,
            value=float(controls.get("standaloneStaffSupportCost") or 0.0),
            step=10.0,
            format="%.0f",
            help="AUD 2019 per person screened.",
        )
        controls["sharedAttributionShare"] = st.slider(
            "Share of common programme costs attributed to LTBI screening",
            min_value=0.0,
            max_value=1.0,
            value=float(controls.get("sharedAttributionShare") or 0.0),
            step=0.05,
            format="%.2f",
            help="Scenario assumption only; no reviewed attribution percentage has been supplied.",
        )
        values = _programme_control_values(controls)
        st.caption(
            "Shared delivery attribution: "
            f"common setup {_money(values['setup_cost'])}; "
            f"share {values['attribution_share'] * 100:.0f}%; "
            f"attributed LTBI setup {_money(values['setup_cost'] * values['attribution_share'])}."
        )
        st.caption(
            "Shared programme attribution changes programme overhead only. It does not bundle clinical pathway visits into test or treatment prices."
        )
    else:
        st.caption("No additional programme overhead has been entered; this does not mean implementation is costless.")
    st.session_state["health_econ_delivery_scenarios"] = controls

    if workspace_state.get("presetConflict"):
        st.warning("The economics configuration changed while the workspace has unsaved edits.")
        conflict_cols = st.columns(3)
        conflict_cols[0].download_button(
            "Download edits before replacing",
            data=assumptions_csv(workspace_state.get("rows") or []),
            file_name=f"{safe_download_stem(scenario_label, 'unsaved_assumption_edits')}.csv",
            mime="text/csv",
        )
        if conflict_cols[1].button("Keep current working edits"):
            st.session_state["health_econ_workspace"] = reconcile_workspace_state(workspace_state, econ_config, action="keep")
            st.rerun()
        if conflict_cols[2].button("Discard and reload from current configuration"):
            st.session_state["health_econ_workspace"] = reconcile_workspace_state(workspace_state, econ_config, action="discard")
            st.rerun()

    overrides = overridden_rows(working_rows)
    if overrides:
        st.markdown("User-defined overrides")
        st.dataframe(arrow_safe_dataframe(overrides), use_container_width=True, hide_index=True)
    else:
        st.caption("No user-defined economic overrides are active.")

    st.markdown("Editable cost inputs")
    st.caption("Change the displayed cost values directly, then use Recalculate economics. Rows marked User-defined are included in exports.")
    working_rows = render_editable_cost_workspace(
        rows=working_rows,
        editor_prefix="health_econ_cost_widget",
    )
    workspace_state = update_workspace_rows(
        st.session_state["health_econ_workspace"],
        rows_from_display_rows(working_rows),
    )
    st.session_state["health_econ_workspace"] = workspace_state

    st.markdown("Discounting, health outcomes and thresholds")
    non_cost_rows = [row for row in working_rows if row.get("category") != "cost"]
    working_rows, workspace_state = render_standard_assumption_editor(
        rows=non_cost_rows,
        working_rows=working_rows,
        editor_prefix="health_econ_standard_assumption_editor",
    )
    st.session_state["health_econ_workspace"] = workspace_state

    active_validation = workspace_state.get("validation") or validate_editable_assumptions(
        working_rows,
        econ_config,
        config=config or {},
    )
    if not active_validation.get("isValidForApplication"):
        st.error("Some economic assumptions must be corrected before recalculation.")
        fatal_rows = fatal_validation_rows(active_validation)
        if fatal_rows:
            st.dataframe(arrow_safe_dataframe(fatal_rows), use_container_width=True, hide_index=True)

    uploaded_assumptions = st.file_uploader(
        "Upload edited assumptions CSV",
        type=["csv"],
        help="Uploaded assumptions are loaded into the working copy and must be recalculated before results update.",
    )
    if uploaded_assumptions is not None and st.button("Load uploaded assumptions"):
        try:
            st.session_state["health_econ_workspace"] = update_workspace_rows(
                workspace_state,
                parse_assumptions_csv(uploaded_assumptions.getvalue()),
            )
            st.success("Uploaded assumptions loaded into the working copy.")
            st.rerun()
        except Exception as exc:
            st.error(f"Could not load assumptions CSV: {exc}")

    if st.checkbox("Show evidence and technical details", value=False):
        current_readiness = assess_current_analysis_economic_readiness(config or {}, econ_config, ledger, working_rows)
        st.dataframe(
            arrow_safe_dataframe(
                [
                    {"Readiness item": "Current cost inputs", "Complete": current_readiness["currentAnalysisCostReady"]},
                    {"Readiness item": "Current DALY inputs", "Complete": current_readiness["currentAnalysisDALYReady"]},
                    {"Readiness item": "Current ICER", "Complete": current_readiness["currentAnalysisICERReady"]},
                    {"Readiness item": "Current NMB", "Complete": current_readiness["currentAnalysisNMBReady"]},
                ]
            ),
            use_container_width=True,
            hide_index=True,
        )
        audit_rows = conversion_audit_rows(working_rows)
        if audit_rows:
            st.markdown("Price-year conversion audit")
            st.dataframe(arrow_safe_dataframe(audit_rows), use_container_width=True, hide_index=True)
        st.markdown("Cost-state audit")
        st.dataframe(arrow_safe_dataframe(cost_state_rows(working_rows)), use_container_width=True, hide_index=True)
        readiness = assess_apy_reference_readiness(
            config or {},
            econ_config,
            econ_config.get("assumptionEvidenceRegistry") or load_apy_evidence_registry(),
        )
        unresolved_rows = [row for row in readiness["readinessRows"] if not row.get("ready")][:12]
        if unresolved_rows:
            st.markdown("Readiness details")
            st.dataframe(arrow_safe_dataframe(unresolved_rows), use_container_width=True, hide_index=True)

    download_cols = st.columns(2)
    download_cols[0].download_button(
        "Download edited assumptions CSV",
        data=assumptions_csv(working_rows),
        file_name=f"{safe_download_stem(scenario_label, 'edited_assumptions')}.csv",
        mime="text/csv",
    )
    download_cols[1].download_button(
        "Download edited assumptions workbook",
        data=assumptions_workbook(working_rows, active_validation),
        file_name=f"{safe_download_stem(scenario_label, 'edited_assumptions')}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

action_cols = st.columns(2)
if action_cols[0].button("Recalculate economics", type="primary", use_container_width=True):
    apply_and_recalculate(
        working_rows=working_rows,
        workspace_state=workspace_state,
        econ_config=econ_config,
        config=config or {},
        results_bundle=results_bundle,
        controls=controls,
    )
if action_cols[1].button("Restore SA Health economic defaults", use_container_width=True):
    load_economics_config(build_unified_working_default_preset()["economicsConfig"])
    st.session_state["health_econ_delivery_scenarios"] = dict(DELIVERY_SCENARIO_DEFAULTS)
    st.session_state["health_econ_cost_widget_version"] = int(st.session_state.get("health_econ_cost_widget_version", 0)) + 1
    if results_bundle and not st.session_state.get("results_stale"):
        run_authoritative_health_economics(results_bundle, st.session_state["economics_config"])
    st.rerun()

download_cols = st.columns(2)
download_cols[0].download_button(
    "Download economics assumptions JSON",
    data=economics_assumptions_json(econ_config),
    file_name=f"{safe_download_stem(scenario_label, 'economics_assumptions')}.json",
    mime="application/json",
)
if econ_results:
    download_cols[1].download_button(
        "Download economics summary CSV",
        data=economics_summary_csv(econ_results),
        file_name=f"{safe_download_stem(scenario_label, 'economics_summary')}.csv",
        mime="text/csv",
    )
else:
    download_cols[1].button("Download economics summary CSV", disabled=True)

if st.session_state.get("dirty_economics") and econ_results:
    st.warning("Economic results are stale because economics inputs changed after the last recalculation.")
if econ_results and (econ_results.get("metadata") or {}).get("isProvisional"):
    st.warning("Economic outputs are provisional and should not be interpreted as clinician-ready cost-effectiveness conclusions.")
if econ_results and econ_results.get("warnings"):
    st.warning("; ".join(str(item) for item in econ_results.get("warnings") or []))

with st.expander("Limitations and methods", expanded=False):
    st.markdown(
        """
        - Health-economic calculations use the completed event ledger from the current screening analysis; cost edits do not rerun epidemiology.
        - DALY-based results are provisional pending evidence review.
        - Programme setup, running, travel, outreach and staff-support costs require SA Health local costing.
        - Net monetary benefit and probability cost-effective remain unavailable unless a reviewed willingness-to-pay threshold is supplied.
        - Scenario comparisons vary economic assumptions only and should not be interpreted as evidence that missing local costs are truly zero.
        """
    )
st.stop()

st.subheader("Analysis status")
status_rows = [
    {
        "Item": "Screening outcomes",
        "Status": "Current and available" if can_run else "Run the screening analysis before recalculating health economics",
    },
    {
        "Item": "Economic assumptions",
        "Status": "SA Health working defaults" if override_count == 0 else f"{override_count} user-defined override(s)",
    },
    {
        "Item": "Epidemiological basis",
        "Status": EXPERIMENTAL_STATUS_LABEL if experimental_ledger else "SA Health report reference",
    },
    {"Item": "Perspective", "Status": (econ_config.get("metadata") or {}).get("perspective", "")},
    {
        "Item": "Currency and price year",
        "Status": f"{(econ_config.get('metadata') or {}).get('targetCurrency') or (econ_config.get('metadata') or {}).get('currencyCode')} {(econ_config.get('metadata') or {}).get('targetPriceYear') or (econ_config.get('metadata') or {}).get('priceYear')}",
    },
    {
        "Item": "Primary discount rate",
        "Status": str((econ_config.get("discounting") or {}).get("primaryDisplayedRate") or (econ_config.get("discounting") or {}).get("selectedAnnualRate") or ""),
    },
]
st.dataframe(arrow_safe_dataframe(status_rows), use_container_width=True, hide_index=True)
if experimental_ledger:
    st.warning(EXPERIMENTAL_ECONOMICS_GUARD_MESSAGE)
    st.info(
        "You can still download the current economic-assumption configuration, "
        "but report-facing health-economic results and scenario comparisons are disabled for this ledger."
    )
    st.download_button(
        "Download economics assumptions JSON",
        data=economics_assumptions_json(econ_config),
        file_name=f"{safe_download_stem(scenario_label, 'experimental_economics_assumptions')}.json",
        mime="application/json",
    )
    st.stop()
st.info(
    "Economic changes recalculate costs and DALYs from the current screening outcomes. "
    "Change population, testing, treatment or targeting on Set up, then run the analysis again."
)
if any(row.get("assumptionId") in PROGRAMME_UNCOSTED_IDS and row.get("currentValue") in (0, 0.0, "0", "0.0") for row in working_rows):
    st.warning("Programme setup, running, travel, outreach and staff-support costs have not yet been locally costed.")

st.subheader("Headline economic results")
if econ_results:
    st.dataframe(
        arrow_safe_dataframe(headline_rows(econ_results=econ_results, results_bundle=results_bundle)),
        use_container_width=True,
        hide_index=True,
    )
else:
    st.info("Recalculate economics to show headline results for the current screening outcomes.")

st.subheader("Cost breakdown and budget impact")
if econ_results:
    cost_rows = cost_category_rows(econ_results)
    if cost_rows:
        st.markdown("Incremental cost categories")
        st.dataframe(arrow_safe_dataframe(cost_rows), use_container_width=True, hide_index=True)
    budget_rows = budget_impact_rows(econ_results)
    if budget_rows:
        st.markdown("Annual budget impact")
        st.line_chart(
            pd.DataFrame(budget_rows).set_index("Year")[
                ["Intervention delivery expenditure", "Comparator active-TB care", "Intervention active-TB care"]
            ],
            use_container_width=True,
        )
        st.dataframe(arrow_safe_dataframe(budget_rows), use_container_width=True, hide_index=True)
    st.caption(
        "Programme and pathway expenditure is concentrated early; active-TB care offsets occur over follow-up. "
        "Exact timing remains provisional where local programme scheduling is not yet specified."
    )
else:
    st.info("Cost categories and annual budget impact will appear after economics are recalculated.")

st.subheader("Economic scenario comparison")
st.caption(
    "These scenarios reuse the same screening outcomes. Changing only economic assumptions does not rerun the epidemiological analysis."
)
if st.button("Compare economic scenarios using current screening outcomes", disabled=not can_run):
    try:
        st.session_state["economic_scenario_comparison"] = build_same_ledger_economic_scenario_comparison(
            results_bundle,
            econ_config,
        )
        st.success("Economic scenarios recalculated from the current screening outcomes.")
    except Exception as exc:
        st.error(f"Economic scenario comparison failed: {exc}")
economic_scenario_comparison = st.session_state.get("economic_scenario_comparison")
if economic_scenario_comparison:
    st.dataframe(
        arrow_safe_dataframe(economic_scenario_comparison.get("rows") or []),
        use_container_width=True,
        hide_index=True,
    )
    st.caption(
        "The setup-cost scenario is illustrative only. The bundled pathway scenario tests possible overlap "
        "between pathway, test and treatment costs."
    )

with st.expander("View or change economic assumptions", expanded=False):
    st.caption(
        "Blank entries are not converted to zero. User-edited rows are labelled User-defined and must pass validation before recalculation."
    )
    st.info(
        "Cost states: numerical costs are included directly; bundled/absorbed means the activity occurs but is not costed separately; "
        "reviewed exclusion means outside the selected perspective; not locally costed means evidence is still needed; compatibility zero is a placeholder."
    )
    preset_cols = st.columns(3)
    if preset_cols[0].button("Restore SA Health economic defaults", type="primary"):
        load_economics_config(build_unified_working_default_preset()["economicsConfig"])
        st.rerun()
    if preset_cols[1].button("Load Dale 2019 AUD working defaults"):
        load_economics_config(
            backend.economics_preset_dale2019_aud(_selected_regimen(config)),
            mark_workspace_unsaved=True,
        )
        st.rerun()
    if preset_cols[2].button("Load blank economics defaults"):
        load_economics_config(backend.default_economics_config())
        st.rerun()

    if workspace_state.get("presetConflict"):
        st.warning("The economics configuration changed while the workspace has unsaved edits.")
        conflict_cols = st.columns(3)
        conflict_cols[0].download_button(
            "Download edits before replacing",
            data=assumptions_csv(workspace_state.get("rows") or []),
            file_name=f"{safe_download_stem(scenario_label, 'unsaved_assumption_edits')}.csv",
            mime="text/csv",
        )
        if conflict_cols[1].button("Keep current working edits"):
            st.session_state["health_econ_workspace"] = reconcile_workspace_state(
                workspace_state,
                econ_config,
                action="keep",
            )
            st.rerun()
        if conflict_cols[2].button("Discard and reload from current configuration"):
            st.session_state["health_econ_workspace"] = reconcile_workspace_state(
                workspace_state,
                econ_config,
                action="discard",
            )
            st.rerun()

    overrides = overridden_rows(working_rows)
    if overrides:
        st.markdown("User-defined overrides")
        st.dataframe(arrow_safe_dataframe(overrides), use_container_width=True, hide_index=True)
    else:
        st.caption("No user-defined economic overrides are active.")

    working_rows, workspace_state = render_standard_assumption_editor(
        rows=working_rows,
        working_rows=working_rows,
        editor_prefix="health_econ_standard_assumption_editor",
    )
    st.session_state["health_econ_workspace"] = workspace_state

    recalc_disabled = not can_run
    if st.button("Recalculate economics using current screening outcomes", disabled=recalc_disabled):
        apply_and_recalculate(
            working_rows=working_rows,
            workspace_state=workspace_state,
            econ_config=econ_config,
            config=config or {},
            results_bundle=results_bundle,
            controls=controls,
        )

    active_validation = workspace_state.get("validation") or validate_editable_assumptions(
        working_rows,
        econ_config,
        config=config or {},
    )
    if active_validation.get("isValidForApplication"):
        st.success("Economic assumptions are structurally safe to apply.")
    else:
        st.error("Some economic assumptions must be corrected before recalculation.")
        fatal_rows = fatal_validation_rows(active_validation)
        if fatal_rows:
            st.dataframe(arrow_safe_dataframe(fatal_rows), use_container_width=True, hide_index=True)

    uploaded_assumptions = st.file_uploader(
        "Upload edited assumptions CSV",
        type=["csv"],
        help="Uploaded assumptions are loaded into the working copy and must be recalculated before results update.",
    )
    if uploaded_assumptions is not None and st.button("Load uploaded assumptions"):
        try:
            st.session_state["health_econ_workspace"] = update_workspace_rows(
                workspace_state,
                parse_assumptions_csv(uploaded_assumptions.getvalue()),
            )
            st.success("Uploaded assumptions loaded into the working copy.")
            st.rerun()
        except Exception as exc:
            st.error(f"Could not load assumptions CSV: {exc}")

    if st.checkbox("Show evidence and technical details", value=False):
        current_readiness = assess_current_analysis_economic_readiness(
            config or {},
            econ_config,
            ledger,
            working_rows,
        )
        readiness_flags = [
            {"Readiness item": "Current cost inputs", "Complete": current_readiness["currentAnalysisCostReady"]},
            {"Readiness item": "Current DALY inputs", "Complete": current_readiness["currentAnalysisDALYReady"]},
            {"Readiness item": "Current ICER", "Complete": current_readiness["currentAnalysisICERReady"]},
            {"Readiness item": "Current NMB", "Complete": current_readiness["currentAnalysisNMBReady"]},
            {"Readiness item": "Full strategy library", "Complete": current_readiness["fullStrategyLibraryReady"]},
            {"Readiness item": "Overall reference evidence", "Complete": current_readiness["overallReferenceEvidenceReady"]},
        ]
        st.dataframe(arrow_safe_dataframe(readiness_flags), use_container_width=True, hide_index=True)
        if current_readiness["currentBlockers"]:
            st.warning(f"{len(current_readiness['currentBlockers'])} current-analysis blocker(s) remain.")
            st.dataframe(arrow_safe_dataframe(current_readiness["currentBlockers"]), use_container_width=True, hide_index=True)
        if current_readiness["alternativeStrategyBlockers"]:
            st.markdown("Inputs needed for additional strategy comparisons")
            st.dataframe(arrow_safe_dataframe(current_readiness["alternativeStrategyBlockers"]), use_container_width=True, hide_index=True)
        audit_rows = conversion_audit_rows(working_rows)
        if audit_rows:
            st.markdown("Price-year conversion audit")
            st.dataframe(arrow_safe_dataframe(audit_rows), use_container_width=True, hide_index=True)
        st.markdown("Cost-state audit")
        st.dataframe(arrow_safe_dataframe(cost_state_rows(working_rows)), use_container_width=True, hide_index=True)
        readiness = assess_apy_reference_readiness(
            config or {},
            econ_config,
            econ_config.get("assumptionEvidenceRegistry") or load_apy_evidence_registry(),
        )
        unresolved_rows = [
            row for row in readiness["readinessRows"]
            if not row.get("ready")
        ][:12]
        if unresolved_rows:
            st.markdown("Readiness details")
            st.dataframe(arrow_safe_dataframe(unresolved_rows), use_container_width=True, hide_index=True)

    download_cols = st.columns(2)
    download_cols[0].download_button(
        "Download edited assumptions CSV",
        data=assumptions_csv(working_rows),
        file_name=f"{safe_download_stem(scenario_label, 'edited_assumptions')}.csv",
        mime="text/csv",
    )
    download_cols[1].download_button(
        "Download edited assumptions workbook",
        data=assumptions_workbook(working_rows, active_validation),
        file_name=f"{safe_download_stem(scenario_label, 'edited_assumptions')}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

st.subheader("Downloads and limitations")
download_cols = st.columns(2)
download_cols[0].download_button(
    "Download economics assumptions JSON",
    data=economics_assumptions_json(econ_config),
    file_name=f"{safe_download_stem(scenario_label, 'economics_assumptions')}.json",
    mime="application/json",
)
if econ_results:
    download_cols[1].download_button(
        "Download economics summary CSV",
        data=economics_summary_csv(econ_results),
        file_name=f"{safe_download_stem(scenario_label, 'economics_summary')}.csv",
        mime="text/csv",
    )
else:
    download_cols[1].button("Download economics summary CSV", disabled=True)

if st.session_state.get("dirty_economics") and econ_results:
    st.warning("Economic results are stale because economics inputs changed after the last recalculation.")
if econ_results and (econ_results.get("metadata") or {}).get("isProvisional"):
    st.warning("Economic outputs are provisional and should not be interpreted as clinician-ready cost-effectiveness conclusions.")
if econ_results and econ_results.get("warnings"):
    st.warning("; ".join(str(item) for item in econ_results.get("warnings") or []))
st.caption(
    "Net monetary benefit and probability cost-effective remain unavailable unless a reviewed willingness-to-pay threshold is supplied. "
    "Local setup, running, travel, outreach and staff-support costs still require SA Health review."
)
