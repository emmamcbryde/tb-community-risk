from __future__ import annotations

from copy import deepcopy
from typing import Any

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
from app.state import (
    get_backend,
    init_session_state,
    mark_economics_changed,
    mark_economics_completed,
    record_message,
    sync_backend_status,
)
from engine.apy.evidence import assess_apy_reference_readiness, load_apy_evidence_registry
from engine.apy.sa_health_reference_package import build_same_ledger_economic_scenario_comparison
from engine.apy.working_defaults import build_unified_working_default_preset


init_session_state()
st.session_state["apy_backend_name"] = "python_apy"
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


def _money(value: Any, *, saving: bool = False) -> str:
    number = _number(value)
    if number is None:
        return "Unavailable"
    if saving:
        return f"AUD {abs(number):,.0f} saving"
    if number < 0:
        return f"AUD {abs(number):,.0f} saving"
    return f"AUD {number:,.0f}"


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


def _classification(econ_results: dict[str, Any] | None) -> str:
    row = _summary_row(econ_results, "primaryICER_ratioOfMeans")
    label = str((row or {}).get("classification") or "").strip()
    if label == "dominant":
        return (
            "Better modelled health outcomes and lower included health-system costs; "
            "provisional because local programme costs remain unresolved."
        )
    return label or "Unavailable"


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
        econ = backend.run_economics(results_bundle, econ_config)
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
        run_authoritative_health_economics(results_bundle, updated_config)
        st.success(
            f"Economic results recalculated with {len(overridden_rows(applied_state['rows']))} user-defined override(s). "
            "Screening outcomes were not rerun."
        )
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
workspace_state = reconcile_workspace_state(
    st.session_state.get("health_econ_workspace"),
    econ_config,
    registry=econ_config.get("assumptionEvidenceRegistry"),
)
st.session_state["health_econ_workspace"] = workspace_state
working_rows = workspace_state["rows"]
override_count = len(overridden_rows(working_rows))

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
