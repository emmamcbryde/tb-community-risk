from __future__ import annotations

from copy import deepcopy

import streamlit as st

from app.display import (
    arrow_safe_dataframe,
    economics_assumptions_json,
    economics_summary_csv,
    safe_download_stem,
)
from app.epidemiology_inputs import (
    apply_ltbi_state_assumption_update,
    fraction_to_percent,
)
from app.state import (
    get_backend,
    init_session_state,
    mark_economics_changed,
    mark_economics_completed,
    record_message,
    sync_backend_status,
)
from app.health_economics_inputs import (
    ADVANCED_EVIDENCE_COLUMNS,
    DERIVED_CONVERSION_COLUMNS,
    INCLUSION_LABEL_TO_CODE,
    STATUS_LABEL_TO_CODE,
    assess_current_analysis_economic_readiness,
    assumptions_csv,
    assumptions_workbook,
    apply_assumptions_to_economics_config,
    conversion_audit_rows,
    editable_assumption_rows,
    fatal_validation_rows,
    group_rows,
    mark_workspace_applied,
    mark_workspace_validated,
    new_workspace_state,
    ordered_editor_rows,
    parse_assumptions_csv,
    reconcile_workspace_state,
    rows_from_display_rows,
    update_workspace_rows,
    validate_editable_assumptions,
)
from engine.apy.costing import normalise_cost_table
from engine.apy.evidence import assess_apy_reference_readiness, load_apy_evidence_registry
from engine.apy.economics import update_cost_item_original_values_from_legacy_fields
from engine.apy.ltbi_state import resolve_ltbi_state_assumptions
from engine.apy.sa_health_reference_package import build_same_ledger_economic_scenario_comparison


init_session_state()
st.session_state["apy_backend_name"] = "python_apy"
backend = get_backend()

st.title("Health Economics")
st.caption("Review economic assumptions, recalculate costs from the current screening outcomes and inspect results.")
st.info(
    "Changing economic inputs recalculates costs, DALYs and economic ratios from the current screening outcomes. "
    "It does not rerun the screening model unless you change population, testing, treatment or targeting inputs."
)

ECONOMICS_WIDGET_KEYS = [
    "econ_currency_code",
    "econ_price_year",
    "econ_location_label",
    "econ_perspective",
    "econ_target_currency",
    "econ_target_price_year",
    "econ_discount_rate",
    "econ_threshold_value",
    "econ_threshold_currency",
    "econ_threshold_year",
    "econ_threshold_source",
    "econ_test_igra",
    "econ_test_tst",
    "econ_regimen_3hp",
    "econ_regimen_4r",
    "econ_regimen_3hr",
    "econ_regimen_6h",
    "econ_regimen_9h",
    "econ_false_positive_incremental",
    "econ_active_tb_cost",
    "econ_setup_total",
    "econ_running_total",
]

ECONOMICS_WIDGET_FIELDS = {
    "econ_currency_code": ("metadata", "currencyCode"),
    "econ_price_year": ("metadata", "priceYear"),
    "econ_location_label": ("metadata", "locationLabel"),
    "econ_perspective": ("metadata", "perspective"),
    "econ_target_currency": ("metadata", "targetCurrency"),
    "econ_target_price_year": ("metadata", "targetPriceYear"),
    "econ_discount_rate": ("discounting", "selectedAnnualRate"),
    "econ_threshold_value": ("threshold", "value"),
    "econ_threshold_currency": ("threshold", "currency"),
    "econ_threshold_year": ("threshold", "referenceYear"),
    "econ_threshold_source": ("threshold", "source"),
    "econ_test_igra": ("costs", "test", "IGRA"),
    "econ_test_tst": ("costs", "test", "TST"),
    "econ_regimen_3hp": ("costs", "regimen", "x3HP"),
    "econ_regimen_4r": ("costs", "regimen", "x4R"),
    "econ_regimen_3hr": ("costs", "regimen", "x3HR"),
    "econ_regimen_6h": ("costs", "regimen", "x6H"),
    "econ_regimen_9h": ("costs", "regimen", "x9H"),
    "econ_false_positive_incremental": (
        "costs",
        "falsePositiveIncrementalPerPerson",
    ),
    "econ_active_tb_cost": ("costs", "activeTBDiseasePerCase"),
    "econ_setup_total": ("costs", "programSetupTotal"),
    "econ_running_total": ("costs", "programRunningTotal"),
}


def ensure_nested(config: dict, *path: str) -> dict:
    current = config
    for key in path:
        value = current.get(key)
        if not isinstance(value, dict):
            value = {}
            current[key] = value
        current = value
    return current


def optional_number_input(label: str, key: str) -> None:
    st.text_input(label, key=key)


def widget_value(value: object) -> str:
    if isinstance(value, list) and not value:
        return ""
    if value is None:
        return ""
    return str(value)


def sync_econ_widgets_from_config(config: dict) -> None:
    for key, path in ECONOMICS_WIDGET_FIELDS.items():
        value = nested_get(config, path)
        st.session_state[key] = widget_value(value)


def sync_econ_widgets_if_missing(config: dict) -> None:
    if any(key not in st.session_state for key in ECONOMICS_WIDGET_KEYS):
        sync_econ_widgets_from_config(config)


def parse_optional_number(key: str, integer: bool = False) -> int | float | list:
    raw = str(st.session_state.get(key, "")).strip()
    if raw == "":
        return []
    value = float(raw)
    if integer and value.is_integer():
        return int(value)
    return value


def economics_config_from_widgets(base_config: dict) -> dict:
    config = deepcopy(base_config)
    ensure_nested(config, "metadata")
    ensure_nested(config, "discounting")
    ensure_nested(config, "threshold")
    ensure_nested(config, "costs")
    ensure_nested(config, "costs", "test")
    ensure_nested(config, "costs", "regimen")
    config["metadata"]["currencyCode"] = str(
        st.session_state.get("econ_currency_code", "")
    ).strip()
    config["metadata"]["priceYear"] = str(
        st.session_state.get("econ_price_year", "")
    ).strip()
    config["metadata"]["locationLabel"] = str(
        st.session_state.get("econ_location_label", "")
    ).strip()
    config["metadata"]["perspective"] = str(
        st.session_state.get("econ_perspective", "")
    ).strip()
    config["metadata"]["targetCurrency"] = str(
        st.session_state.get("econ_target_currency", "")
    ).strip()
    config["metadata"]["targetPriceYear"] = str(
        st.session_state.get("econ_target_price_year", "")
    ).strip()
    config["discounting"]["selectedAnnualRate"] = parse_optional_number(
        "econ_discount_rate"
    )
    config["threshold"]["value"] = parse_optional_number("econ_threshold_value")
    config["threshold"]["currency"] = str(
        st.session_state.get("econ_threshold_currency", "")
    ).strip()
    config["threshold"]["referenceYear"] = parse_optional_number(
        "econ_threshold_year",
        integer=True,
    )
    config["threshold"]["source"] = str(
        st.session_state.get("econ_threshold_source", "")
    ).strip()
    config["costs"]["test"]["IGRA"] = parse_optional_number("econ_test_igra")
    config["costs"]["test"]["TST"] = parse_optional_number("econ_test_tst")
    config["costs"]["regimen"]["x3HP"] = parse_optional_number("econ_regimen_3hp")
    config["costs"]["regimen"]["x4R"] = parse_optional_number("econ_regimen_4r")
    config["costs"]["regimen"]["x3HR"] = parse_optional_number("econ_regimen_3hr")
    config["costs"]["regimen"]["x6H"] = parse_optional_number("econ_regimen_6h")
    config["costs"]["regimen"]["x9H"] = parse_optional_number("econ_regimen_9h")
    config["costs"]["falsePositiveIncrementalPerPerson"] = parse_optional_number(
        "econ_false_positive_incremental"
    )
    config["costs"]["activeTBDiseasePerCase"] = parse_optional_number(
        "econ_active_tb_cost"
    )
    config["costs"]["programSetupTotal"] = parse_optional_number("econ_setup_total")
    config["costs"]["programRunningTotal"] = parse_optional_number("econ_running_total")
    return update_cost_item_original_values_from_legacy_fields(config)


def nested_get(config: dict, path: tuple[str, ...]) -> object:
    current: object = config
    for key in path:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current


def economics_overview_rows(config: dict) -> list[dict[str, object]]:
    metadata = config.get("metadata", {})
    costs = config.get("costs", {})
    test_costs = costs.get("test", {})
    regimen_costs = costs.get("regimen", {})
    return [
        {"field": "currencyCode", "value": metadata.get("currencyCode")},
        {"field": "priceYear", "value": metadata.get("priceYear")},
        {"field": "locationLabel", "value": metadata.get("locationLabel")},
        {"field": "perspective", "value": metadata.get("perspective")},
        {"field": "targetCurrency", "value": metadata.get("targetCurrency")},
        {"field": "targetPriceYear", "value": metadata.get("targetPriceYear")},
        {"field": "discountRate", "value": config.get("discounting", {}).get("selectedAnnualRate")},
        {"field": "primaryHealthOutcome", "value": config.get("healthOutcome", {}).get("primary")},
        {"field": "threshold.value", "value": config.get("threshold", {}).get("value")},
        {"field": "threshold.currency", "value": config.get("threshold", {}).get("currency")},
        {"field": "threshold.referenceYear", "value": config.get("threshold", {}).get("referenceYear")},
        {"field": "threshold.source", "value": config.get("threshold", {}).get("source")},
        {"field": "test.IGRA", "value": test_costs.get("IGRA")},
        {"field": "test.TST", "value": test_costs.get("TST")},
        {"field": "regimen.3HP", "value": regimen_costs.get("x3HP")},
        {"field": "regimen.4R", "value": regimen_costs.get("x4R")},
        {"field": "regimen.3HR", "value": regimen_costs.get("x3HR")},
        {"field": "regimen.6H", "value": regimen_costs.get("x6H")},
        {"field": "regimen.9H", "value": regimen_costs.get("x9H")},
        {
            "field": "falsePositiveIncrementalPerPerson",
            "value": costs.get("falsePositiveIncrementalPerPerson"),
        },
        {
            "field": "activeTBDiseasePerCase",
            "value": costs.get("activeTBDiseasePerCase"),
        },
        {"field": "programSetupTotal", "value": costs.get("programSetupTotal")},
        {"field": "programRunningTotal", "value": costs.get("programRunningTotal")},
    ]


def load_economics_config(new_config: dict, *, mark_workspace_unsaved: bool = False) -> None:
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
    sync_econ_widgets_from_config(new_config)


def run_authoritative_health_economics(results_bundle: dict, econ_config: dict) -> None:
    try:
        econ = backend.run_economics(results_bundle, econ_config)
        st.session_state["economics_results"] = econ
        mark_economics_completed()
        sync_backend_status(backend.status())
        st.success("Health-economic analysis completed.")
    except Exception as exc:
        message = f"Health-economic analysis failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)


def render_assumption_editor(
    *,
    rows: list[dict],
    working_rows: list[dict],
    group_names: list[str],
    editor_prefix: str,
    show_advanced_columns: bool,
) -> tuple[list[dict], dict]:
    hidden_columns = ["reviewStatus", "inclusionStatus"]
    if not show_advanced_columns:
        hidden_columns.extend(ADVANCED_EVIDENCE_COLUMNS)
    tabs = st.tabs(group_names)
    latest_rows = list(working_rows)
    for tab, group_name in zip(tabs, group_names):
        with tab:
            group_data = group_rows(rows, group_name)
            if not group_data:
                st.info("No rows in this section.")
                continue
            display_rows = ordered_editor_rows(group_data, advanced=show_advanced_columns)
            edited = st.data_editor(
                arrow_safe_dataframe(display_rows),
                width="stretch",
                hide_index=True,
                num_rows="fixed",
                key=f"{editor_prefix}_{group_name}",
                disabled=[
                    "assumptionId",
                    "category",
                    "description",
                    *hidden_columns,
                    *DERIVED_CONVERSION_COLUMNS,
                    "convertedTargetYearCost",
                    "validationMessage",
                ],
                column_config={
                    **(
                        {
                            "assumptionId": None,
                            "category": None,
                            "sourceLocation": None,
                            "targetCurrency": None,
                            "targetPriceYear": None,
                            "sourceYearIndexValue": None,
                            "targetYearIndexValue": None,
                            "bundledIntoAssumptionId": None,
                            "doubleCountingGroup": None,
                            "unresolvedReason": None,
                        }
                        if not show_advanced_columns
                        else {}
                    ),
                    "reviewStatusLabel": st.column_config.SelectboxColumn(
                        "Review status",
                        options=list(STATUS_LABEL_TO_CODE),
                    ),
                    "inclusionStatusLabel": st.column_config.SelectboxColumn(
                        "Inclusion status",
                        options=list(INCLUSION_LABEL_TO_CODE),
                    ),
                    "provisional": st.column_config.CheckboxColumn("Provisional"),
                },
            )
            edited_records = edited.to_dict(orient="records") if hasattr(edited, "to_dict") else list(edited)
            edited_by_id = {row.get("assumptionId"): row for row in edited_records}
            for idx, row in enumerate(latest_rows):
                replacement = edited_by_id.get(row.get("assumptionId"))
                if replacement is not None:
                    latest_rows[idx] = {**row, **replacement}
    updated_state = update_workspace_rows(
        st.session_state["health_econ_workspace"],
        rows_from_display_rows(latest_rows),
    )
    return updated_state["rows"], updated_state


cols = st.columns(3)
if cols[0].button("Load economics defaults", type="primary"):
    try:
        load_economics_config(backend.default_economics_config())
        st.rerun()
    except Exception as exc:
        message = f"Could not load economics defaults: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

if cols[1].button("Load KWAB150 preset"):
    try:
        load_economics_config(backend.economics_preset_kwab150())
        st.rerun()
    except Exception as exc:
        message = f"Could not load KWAB150 preset: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

if cols[2].button("Load Dale 2019 AUD working defaults"):
    try:
        current_config = st.session_state.get("config") or {}
        selected_regimen = (
            current_config.get("regimen")
            or current_config.get("preventiveRegimen")
            or (current_config.get("treatment") or {}).get("regimen")
            or "3HP"
        )
        load_economics_config(
            backend.economics_preset_dale2019_aud(selected_regimen),
            mark_workspace_unsaved=True,
        )
        st.session_state["dale_2019_working_defaults_loaded"] = True
        st.rerun()
    except Exception as exc:
        message = f"Could not load Dale 2019 AUD working defaults: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

config = st.session_state.get("config")
results_bundle = st.session_state.get("results_bundle")
econ_config = st.session_state.get("economics_config")
scenario_label = (results_bundle or {}).get("metadata", {}).get("scenarioLabel")

if (econ_config or {}).get("metadata", {}).get("presetName") == "Dale 2019 AUD working defaults":
    st.success("Dale 2019 AUD working defaults loaded.")
    st.dataframe(
        arrow_safe_dataframe(
            [
                {"Item": "Costs", "Value": "2019 AUD"},
                {"Item": "Perspective", "Value": "Australian health-care system"},
                {"Item": "Discounting", "Value": "3% primary and 0% comparison"},
                {"Item": "Active-TB DALY", "Value": "0.333 for 0.5 years plus mortality"},
                {"Item": "Programme implementation costs", "Value": "Excluded in this working default"},
                {"Item": "Overall APY evidence status", "Value": "Still provisional"},
            ]
        ),
        width="content",
        hide_index=True,
    )
    st.warning(
        "This working default is likely optimistic where a new programme "
        "requires additional setup, staffing, engagement or delivery expenditure."
    )

if config:
    st.subheader("Recent versus remote LTBI assumption")
    ltbi_state = resolve_ltbi_state_assumptions(config)
    current_recent = ltbi_state.get("baselineRecentLTBIProportion")
    current_review_label = (
        "Reviewed model-derived estimate"
        if ltbi_state.get("baselineRecentLTBIProportionStatus") == "model_derived_reviewed"
        else "Reviewed direct evidence"
    )
    st.dataframe(
        arrow_safe_dataframe(
            [
                {
                    "Assumption": "Baseline recent-LTBI proportion",
                    "Value": current_recent,
                },
                {
                    "Assumption": "Source",
                    "Value": ltbi_state.get("baselineRecentLTBIProportionSource") or "",
                },
                {
                    "Assumption": "Provisional",
                    "Value": ltbi_state.get("provisional"),
                },
                {
                    "Assumption": "Warning",
                    "Value": ltbi_state.get("warning") or "",
                },
            ]
        ),
        width="stretch",
        hide_index=True,
    )
    with st.form("ltbi_recent_reviewed_assumption"):
        st.write(
            "Enter a reviewed recent-LTBI proportion or a documented model-derived "
            "estimate. This updates the authoritative LTBI-state assumption."
        )
        recent_percent = st.number_input(
            "Baseline recent-LTBI proportion (%)",
            min_value=0.0,
            max_value=100.0,
            value=float(fraction_to_percent(current_recent, 0.0) or 0.0),
            step=1.0,
            format="%.2f",
        )
        review_method = st.selectbox(
            "How was this value established?",
            ["Reviewed direct evidence", "Reviewed model-derived estimate"],
            index=0
            if current_review_label == "Reviewed direct evidence"
            else 1,
        )
        recent_source = st.text_input(
            "Source or derivation reference",
            value=str(ltbi_state.get("baselineRecentLTBIProportionSource") or ""),
        )
        recent_notes = st.text_area(
            "Notes",
            value=str(ltbi_state.get("notes") or ""),
        )
        save_recent = st.form_submit_button("Save recent-LTBI assumption")
    if save_recent:
        if not str(recent_source).strip():
            st.error("Provide a source or documented derivation before saving a reviewed assumption.")
        else:
            status = (
                "model_derived_reviewed"
                if review_method == "Reviewed model-derived estimate"
                else "configured_reviewed"
            )
            updated = apply_ltbi_state_assumption_update(
                config,
                baseline_recent_percent=float(recent_percent),
                transition_rate_per_year=float(
                    ltbi_state["recentToRemoteTransitionRatePerYear"]
                ),
                source=str(recent_source).strip(),
                status=status,
                notes=str(recent_notes).strip(),
            )
            st.session_state["config"] = updated
            st.session_state.pop("recent_ltbi_run_route", None)
            st.success("Recent-LTBI assumption saved.")
            st.rerun()

if not econ_config:
    st.info("Load economics defaults or the KWAB150 preset to begin.")
    st.stop()

sync_econ_widgets_if_missing(econ_config)

if st.session_state.get("dirty_economics") and st.session_state.get("economics_results"):
    st.warning("Economics results are stale because model results or economics inputs changed.")

ledger = (results_bundle or {}).get("technical", {}).get("eventLedger", {}) if isinstance(results_bundle, dict) else {}

st.subheader("Inputs required for this analysis")
st.caption(
    "Edit the assumptions that apply to the currently selected test, regimen and economic outcome. "
    "The source registry file is not changed from this page."
)

workspace_state = reconcile_workspace_state(
    st.session_state.get("health_econ_workspace"),
    econ_config,
    registry=econ_config.get("assumptionEvidenceRegistry"),
)
st.session_state["health_econ_workspace"] = workspace_state
if st.session_state.pop("health_econ_apply_message", None):
    st.success("Validated assumptions applied to the current analysis.")
if workspace_state.get("hasUnsavedEdits"):
    st.warning("Unsaved changes")
if workspace_state.get("presetConflict"):
    st.warning(
        "The economics preset or configuration changed while the assumptions "
        "workspace has unsaved edits. Choose how to handle the working copy."
    )
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
    if conflict_cols[2].button("Discard and reload from the new preset"):
        st.session_state["health_econ_workspace"] = reconcile_workspace_state(
            workspace_state,
            econ_config,
            action="discard",
        )
        st.rerun()

uploaded_assumptions = st.file_uploader(
    "Upload edited assumptions CSV",
    type=["csv"],
    help="Uploaded assumptions are loaded into the working copy and must be validated before application.",
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

if st.button("Reset working copy from current analysis"):
    st.session_state["health_econ_workspace"] = reconcile_workspace_state(
        workspace_state,
        econ_config,
        action="reset",
    )
    st.rerun()
if st.button("Discard working edits"):
    st.session_state["health_econ_workspace"] = reconcile_workspace_state(
        workspace_state,
        econ_config,
        action="discard",
    )
    st.rerun()

working_rows = workspace_state["rows"]
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
st.dataframe(arrow_safe_dataframe(readiness_flags), width="content", hide_index=True)
if st.button("Edit blocking assumptions", type="primary"):
    st.session_state["health_econ_show_blocking_editor"] = True

current_ids = set(current_readiness["currentApplicableAssumptionIds"])
current_rows = [row for row in working_rows if row.get("assumptionId") in current_ids]
show_advanced_columns = st.checkbox("Advanced evidence columns", value=False)
working_rows, workspace_state = render_assumption_editor(
    rows=current_rows,
    working_rows=working_rows,
    group_names=["Costs", "DALYs", "Threshold", "Epidemiology blockers"],
    editor_prefix="health_econ_current_assumption_editor",
    show_advanced_columns=show_advanced_columns,
)
st.session_state["health_econ_workspace"] = workspace_state
st.markdown("Price-year conversion audit")
current_rows = [row for row in working_rows if row.get("assumptionId") in current_ids]
audit_rows = conversion_audit_rows(current_rows)
if audit_rows:
    st.dataframe(arrow_safe_dataframe(audit_rows), width="stretch", hide_index=True)
else:
    st.info("No current cost rows require price-year conversion audit.")

current_readiness = assess_current_analysis_economic_readiness(
    config or {},
    econ_config,
    ledger,
    working_rows,
)
if current_readiness["currentBlockers"]:
    st.warning(f"{len(current_readiness['currentBlockers'])} current-analysis blockers require review.")
    for group_name in [
        "Current selected test and regimen",
        "Programme and delivery costs",
        "Active-TB care",
        "DALY assumptions",
        "Threshold for NMB",
    ]:
        group_blockers = [
            blocker for blocker in current_readiness["currentBlockers"]
            if blocker.get("group") == group_name
        ]
        if not group_blockers:
            continue
        with st.expander(group_name, expanded=True):
            st.dataframe(arrow_safe_dataframe(group_blockers), width="stretch", hide_index=True)
else:
    st.success("No current-analysis assumption blockers were found.")

with st.expander("All assumptions and alternative strategies", expanded=False):
    if current_readiness["alternativeStrategyBlockers"]:
        st.markdown("Inputs needed for additional strategy comparisons")
        st.dataframe(
            arrow_safe_dataframe(current_readiness["alternativeStrategyBlockers"]),
            width="stretch",
            hide_index=True,
        )
    working_rows, workspace_state = render_assumption_editor(
        rows=working_rows,
        working_rows=working_rows,
        group_names=["Costs", "DALYs", "Threshold", "Epidemiology blockers"],
        editor_prefix="health_econ_all_assumption_editor",
        show_advanced_columns=show_advanced_columns,
    )
    st.session_state["health_econ_workspace"] = workspace_state

validation_report = workspace_state.get("validation")
cols = st.columns(4)
if cols[0].button("Validate assumptions", type="primary"):
    validation_report = validate_editable_assumptions(
        working_rows,
        econ_config,
        config=config or {},
    )
    st.session_state["health_econ_workspace"] = mark_workspace_validated(
        workspace_state,
        validation_report,
    )
    st.rerun()

apply_disabled = validation_report is None or not bool(validation_report.get("isValidForApplication"))
if validation_report:
    if validation_report.get("isValidForApplication"):
        st.success("Assumptions are structurally safe to apply.")
    else:
        st.error("Assumptions contain errors that must be corrected before applying.")
        fatal_rows = fatal_validation_rows(validation_report)
        if fatal_rows:
            st.dataframe(arrow_safe_dataframe(fatal_rows), width="stretch", hide_index=True)
if apply_disabled and validation_report is not None:
    fatal_rows = fatal_validation_rows(validation_report)
    if fatal_rows:
        st.caption("Application is disabled because these rows contain fatal validation errors.")
        st.dataframe(arrow_safe_dataframe(fatal_rows), width="stretch", hide_index=True)
if cols[1].button("Apply assumptions to current analysis", disabled=apply_disabled):
    try:
        updated_config = apply_assumptions_to_economics_config(
            econ_config,
            working_rows,
            config=config or {},
        )
        st.session_state["economics_config"] = updated_config
        econ_config = updated_config
        applied_state = new_workspace_state(
            updated_config,
            updated_config.get("assumptionEvidenceRegistry") or working_rows,
        )
        applied_state = mark_workspace_validated(
            applied_state,
            validate_editable_assumptions(applied_state["rows"], updated_config, config=config or {}),
        )
        applied_state = mark_workspace_applied(applied_state, updated_config)
        applied_state["rows"] = editable_assumption_rows(
            updated_config.get("assumptionEvidenceRegistry") or working_rows,
            updated_config,
        )
        applied_state = mark_workspace_applied(applied_state, updated_config)
        st.session_state["health_econ_workspace"] = applied_state
        mark_economics_changed()
        st.session_state["health_econ_apply_message"] = True
        st.rerun()
    except Exception as exc:
        st.error(f"Assumptions were not applied: {exc}")

can_run = bool(config and results_bundle and not st.session_state.get("results_stale"))
if cols[2].button(
    "Recalculate Health Economics",
    disabled=not (can_run and bool((st.session_state.get("health_econ_workspace") or {}).get("applied"))),
):
    run_authoritative_health_economics(results_bundle, econ_config)

active_validation = validation_report or validate_editable_assumptions(
    working_rows,
    econ_config,
    config=config or {},
)
summary = active_validation["summary"]
current_readiness = assess_current_analysis_economic_readiness(
    config or {},
    econ_config,
    ledger,
    working_rows,
)
if not current_readiness["currentAnalysisICERReady"]:
    st.warning("ICER cannot yet be calculated.")
    if current_readiness["currentBlockers"]:
        st.dataframe(
            arrow_safe_dataframe(current_readiness["currentBlockers"]),
            width="stretch",
            hide_index=True,
        )
elif not current_readiness["currentAnalysisNMBReady"]:
    st.info("ICER inputs are complete, but NMB requires a reviewed threshold with matching currency and reference year.")
if current_readiness["currentAnalysisCostReady"] and not current_readiness["currentAnalysisDALYReady"]:
    st.info(
        "Cost inputs are complete but DALY inputs are incomplete. Cost-consequence outputs may be reviewed; ICER remains unavailable."
    )

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

with st.expander("Legacy/developer economic controls", expanded=False):
    st.caption(
        "The assumptions workspace above is the authoritative standard editing route. "
        "These compatibility controls are retained for local development only."
    )
    workspace_applied = bool((st.session_state.get("health_econ_workspace") or {}).get("applied"))
    if workspace_applied:
        st.info("Workspace assumptions have been applied. Legacy controls are disabled to avoid conflicting overrides.")
    metadata = ensure_nested(econ_config, "metadata")
    costs = ensure_nested(econ_config, "costs")
    ensure_nested(econ_config, "costs", "test")
    ensure_nested(econ_config, "costs", "regimen")
    ensure_nested(econ_config, "discounting")
    ensure_nested(econ_config, "threshold")

    with st.form("economics_edits"):
        st.markdown("Metadata")
        st.text_input("Currency code", key="econ_currency_code")
        optional_number_input("Price year", "econ_price_year")
        st.text_input("Location label", key="econ_location_label")
        st.text_input("Economic perspective", key="econ_perspective", disabled=True)
        st.text_input("Target currency", key="econ_target_currency")
        st.text_input("Target price year", key="econ_target_price_year")
        st.selectbox(
            "Discount rate",
            ["0.03", "0.0"],
            index=0 if str(st.session_state.get("econ_discount_rate", "0.03")) != "0.0" else 1,
            key="econ_discount_rate",
        )
        st.markdown("Primary outcome and benchmark")
        st.info("Default health outcome: DALYs averted.")
        optional_number_input("Illustrative GDP-per-capita threshold value", "econ_threshold_value")
        st.text_input("Threshold currency", key="econ_threshold_currency")
        optional_number_input("Threshold reference year", "econ_threshold_year")
        st.text_input("Threshold source", key="econ_threshold_source")
        st.caption("The GDP-per-capita benchmark is illustrative, not an official Australian funding threshold.")
        st.markdown("Test costs")
        optional_number_input("IGRA test cost", "econ_test_igra")
        optional_number_input("TST cost", "econ_test_tst")
        st.markdown("Regimen costs")
        optional_number_input("3HP regimen cost", "econ_regimen_3hp")
        optional_number_input("4R regimen cost", "econ_regimen_4r")
        optional_number_input("3HR regimen cost", "econ_regimen_3hr")
        optional_number_input("6H regimen cost", "econ_regimen_6h")
        optional_number_input("9H regimen cost", "econ_regimen_9h")
        st.markdown("Program and disease costs")
        optional_number_input("False-positive incremental cost per person", "econ_false_positive_incremental")
        optional_number_input("Active TB disease cost per case", "econ_active_tb_cost")
        optional_number_input("Program setup total", "econ_setup_total")
        optional_number_input("Program running total", "econ_running_total")
        submitted = st.form_submit_button("Apply legacy economics edits", disabled=workspace_applied)

    if submitted and not workspace_applied:
        try:
            updated_config = economics_config_from_widgets(econ_config)
            changed = updated_config != econ_config
            st.session_state["economics_config"] = updated_config
            econ_config = updated_config
            st.session_state["health_econ_workspace"] = reconcile_workspace_state(
                st.session_state.get("health_econ_workspace"),
                econ_config,
            )
            if changed:
                mark_economics_changed()
                st.success("Legacy economics edits applied.")
            else:
                st.info("No legacy economics fields changed.")
        except ValueError as exc:
            st.error(f"Invalid economics number: {exc}")

st.subheader("Current Assumptions")
warnings = []
ledger = (results_bundle or {}).get("technical", {}).get("eventLedger", {}) if isinstance(results_bundle, dict) else {}
if results_bundle:
    if ledger and ledger.get("validation", {}).get("isValid") is True:
        st.success(
            "Current screening outcomes are available for economic recalculation."
        )
    else:
        st.warning("Run the screening analysis before recalculating health economics.")
for item in normalise_cost_table(econ_config.get("costItems") or []):
    if item.get("conversionStatus") != "valid":
        warnings.append(f"{item.get('costItemId')}: {item.get('conversionStatus')}")
if econ_config.get("threshold", {}).get("value") in (None, "", []):
    warnings.append("GDP-per-capita threshold value is unresolved.")
if warnings:
    st.warning("Unresolved assumptions: " + "; ".join(warnings))
st.dataframe(
    arrow_safe_dataframe(economics_overview_rows(econ_config)),
    width="content",
    hide_index=True,
)

readiness = assess_apy_reference_readiness(
    config or {},
    econ_config,
    econ_config.get("assumptionEvidenceRegistry") or load_apy_evidence_registry(),
)
st.subheader("Evidence-readiness explanation")
status_rows = [
    {"category": "epidemiology", "ready": readiness["epidemiologyReady"]},
    {"category": "cost", "ready": readiness["costReady"]},
    {"category": "DALY", "ready": readiness["dalyReady"]},
    {"category": "threshold", "ready": readiness["thresholdReady"]},
    {"category": "overall clinician-ready", "ready": readiness["overallClinicianReady"]},
]
st.dataframe(arrow_safe_dataframe(status_rows), width="content", hide_index=True)
unresolved_rows = [
    row
    for row in readiness["readinessRows"]
    if not row.get("ready")
][:12]
if unresolved_rows:
    st.warning("APY reference evidence remains unresolved or provisional.")
    st.dataframe(
        arrow_safe_dataframe(
            [
                {
                    "category": row.get("category"),
                    "assumptionId": row.get("assumptionId"),
                    "status": row.get("reviewStatus"),
                    "source": row.get("sourceCitation"),
                    "unresolvedReason": row.get("unresolvedReason"),
                }
                for row in unresolved_rows
            ]
        ),
        width="content",
        hide_index=True,
    )
st.download_button(
    "Download economics assumptions JSON",
    data=economics_assumptions_json(econ_config),
    file_name=f"{safe_download_stem(scenario_label, 'economics_assumptions')}.json",
    mime="application/json",
)

can_run = bool(config and results_bundle and not st.session_state.get("results_stale"))
if not config:
    st.info("Load a scenario before running economics.")
elif not results_bundle:
    st.info("Run the model before running economics.")
elif st.session_state.get("results_stale"):
    st.warning("Rerun the model before running economics so the economics inputs match current results.")
else:
    semantics = str((config or {}).get("naturalHistorySemantics") or "")
    if semantics == "matlab_v9_implicit_early_late":
        st.info(
            "Epidemiological anchor: frozen APY stochastic compatibility reference. "
            "Economic edits use the current screening outcomes and do not rerun epidemiology."
        )
    else:
        st.warning(
            "This run uses the explicit recent/remote technical scenario, not the "
            "frozen SA Health epidemiological compatibility anchor."
        )

if st.button("Run health economics", type="primary", disabled=not can_run):
    run_authoritative_health_economics(results_bundle, econ_config)

st.subheader("Compare economic scenarios")
st.caption(
    "Compare alternative cost assumptions against the current screening and treatment outcomes. "
    "Changing only economic assumptions does not rerun the epidemiological analysis."
)
scenario_disabled = not can_run
if st.button("Compare economic scenarios using current results", disabled=scenario_disabled):
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
        width="stretch",
        hide_index=True,
    )
    st.caption(
        "Gross delivery expenditure ratios are before active-TB care offsets; net health-system "
        "ratios include active-TB care offsets."
    )
    st.caption(
        "The AUD 500,000 setup case is illustrative only. The bundled pathway case is a sensitivity "
        "for possible overlap between test, regimen and pathway costs."
    )

econ_results = st.session_state.get("economics_results")
if econ_results:
    metadata = econ_results.get("metadata", {})
    if metadata.get("isProvisional"):
        st.warning("Economic outputs are provisional. Do not interpret them as clinician-ready cost-effectiveness conclusions.")
    if econ_results.get("warnings"):
        st.warning("; ".join(str(item) for item in econ_results.get("warnings") or []))
    st.subheader("Health-economic Summary")
    summary_rows = econ_results.get("summaryRows") or []
    if summary_rows:
        status = econ_results.get("status", {})
        validation = econ_results.get("validation", {})
        complete = bool(status.get("isComplete"))
        st.dataframe(arrow_safe_dataframe(summary_rows), width="stretch")
        primary = [
            row for row in summary_rows
            if row.get("metric") in {"incrementalCost", "dalysAverted", "primaryICER_ratioOfMeans"}
            and float(row.get("discountRate", 0.03)) in {0.0, 0.03}
        ]
        if primary and complete:
            st.caption("Primary ICER uses mean paired incremental cost divided by mean paired DALYs averted; replicate ICERs are diagnostic only.")
            st.dataframe(arrow_safe_dataframe(primary), width="stretch")
        elif primary:
            st.markdown("Partial calculations - not a complete economic result")
            st.dataframe(arrow_safe_dataframe(primary), width="stretch")
        if validation:
            st.caption(
                "Complete pairs: "
                f"{validation.get('completePairedReplicates')} of "
                f"{validation.get('totalPairedReplicates')}."
            )
        st.markdown("Downloads")
        st.download_button(
            "Download economics summary CSV",
            data=economics_summary_csv(econ_results),
            file_name=f"{safe_download_stem(scenario_label, 'economics_summary')}.csv",
            mime="text/csv",
        )
    else:
        st.info("No health-economic summary rows were returned.")

    st.subheader("Health-economic Status")
    status = econ_results.get("status", {})
    st.json(
        {
            "last_economics_run_at": st.session_state.get("last_economics_run_at"),
            "isComplete": status.get("isComplete"),
            "missingInputs": status.get("missingInputs"),
            "notCalculated": status.get("notCalculated"),
            "messages": status.get("messages"),
            "partialCalculations": status.get("partialCalculations"),
        },
        expanded=False,
    )
else:
    st.info("Run health economics to enable the summary CSV download.")
