from __future__ import annotations

from copy import deepcopy
from typing import Any

import streamlit as st

from app.display import arrow_safe_dataframe
from app.parameter_workspace import (
    PARAMETER_GROUPS,
    apply_parameter_workspace,
    build_parameter_workspace,
    changed_parameter_count,
    parameter_summary,
    reset_all_parameters,
    reset_parameter_group,
    unified_default_session_state,
    validate_parameter_workspace,
)
from app.run_analysis_controls import TECHNICAL_DEMONSTRATION_ROUTE
from app.state import get_backend, init_session_state, record_message, sync_backend_status
from engine.apy.ltbi_state import resolve_ltbi_state_assumptions


init_session_state()


def _load_unified_defaults(*, show_workspace: bool) -> None:
    state = unified_default_session_state()
    st.session_state["apy_backend_name"] = "python_apy"
    st.session_state["config"] = state["config"]
    st.session_state["economics_config"] = state["economics_config"]
    st.session_state["parameter_workspace"] = state["parameter_workspace"]
    st.session_state["working_default_preset"] = state["working_default_preset"]
    st.session_state["parameter_workspace_visible"] = show_workspace
    st.session_state["parameter_workspace_validation"] = None
    st.session_state["validation_report"] = None
    st.session_state["results_bundle"] = None
    st.session_state["economics_results"] = None
    st.session_state["decision_analysis_results"] = None
    st.session_state["dirty_config"] = False
    st.session_state["dirty_economics"] = False
    st.session_state["results_stale"] = False
    st.session_state["last_economics_run_at"] = ""
    st.session_state["last_run_at"] = ""
    st.session_state.pop("recent_ltbi_run_route", None)
    st.session_state.pop("health_econ_workspace", None)
    sync_backend_status(get_backend().status())


def _current_workspace() -> dict[str, Any] | None:
    config = st.session_state.get("config")
    econ = st.session_state.get("economics_config")
    if not isinstance(config, dict) or not isinstance(econ, dict):
        return None
    workspace = st.session_state.get("parameter_workspace")
    if not isinstance(workspace, dict):
        workspace = build_parameter_workspace(config, econ)
        st.session_state["parameter_workspace"] = workspace
    return workspace


def _rows_from_editor(value: Any) -> list[dict[str, Any]]:
    if hasattr(value, "to_dict"):
        return value.to_dict(orient="records")
    return [dict(row) for row in value]


def _set_workspace_rows(rows: list[dict[str, Any]]) -> None:
    workspace = deepcopy(st.session_state["parameter_workspace"])
    defaults = {row["parameterId"]: row.get("defaultValue") for row in workspace.get("rows") or []}
    for row in rows:
        row["changedFromDefault"] = str(row.get("currentValue")) != str(defaults.get(row.get("parameterId")))
    workspace["rows"] = rows
    workspace["changedCount"] = changed_parameter_count(workspace)
    st.session_state["parameter_workspace"] = workspace


def _render_parameter_workspace() -> None:
    workspace = _current_workspace()
    if not workspace:
        return
    st.subheader("Review or change parameters")
    st.caption(
        "These grouped inputs update the same configuration used by the analysis. "
        "Changing a value does not change its evidence source or review status."
    )
    st.metric("Parameters changed from working defaults", changed_parameter_count(workspace))
    top_cols = st.columns([1, 1, 3])
    if top_cols[0].button("Use all defaults", use_container_width=True):
        _set_workspace_rows(reset_all_parameters(workspace["rows"]))
        st.session_state["parameter_workspace_validation"] = None
        st.rerun()
    if top_cols[1].button("Reset all changes", use_container_width=True):
        _set_workspace_rows(reset_all_parameters(workspace["rows"]))
        st.session_state["parameter_workspace_validation"] = None
        st.rerun()

    edited_rows: list[dict[str, Any]] = []
    tabs = st.tabs(PARAMETER_GROUPS)
    for tab, group in zip(tabs, PARAMETER_GROUPS):
        with tab:
            group_rows = [row for row in workspace["rows"] if row.get("group") == group]
            if st.button(f"Reset this section: {group}", key=f"reset_{group}"):
                _set_workspace_rows(reset_parameter_group(workspace["rows"], group))
                st.session_state["parameter_workspace_validation"] = None
                st.rerun()
            standard_rows = [row for row in group_rows if not row.get("advanced")]
            advanced_rows = [row for row in group_rows if row.get("advanced")]
            shown_columns = [
                "label",
                "currentValue",
                "defaultValue",
                "unit",
                "source",
                "reviewStatus",
                "provisional",
                "changedFromDefault",
                "operationalStatus",
                "notes",
            ]
            disabled_columns = [
                "parameterId",
                "group",
                "label",
                "defaultValue",
                "unit",
                "source",
                "reviewStatus",
                "provisional",
                "changedFromDefault",
                "editableType",
                "validation",
                "notes",
                "advanced",
                "selectOptions",
                "operationalStatus",
            ]
            edited = st.data_editor(
                standard_rows,
                key=f"parameter_editor_{group}",
                use_container_width=True,
                hide_index=True,
                column_order=shown_columns,
                disabled=disabled_columns,
            )
            edited_rows.extend(_rows_from_editor(edited))
            with st.expander("Advanced"):
                advanced_edited = st.data_editor(
                    advanced_rows,
                    key=f"parameter_advanced_editor_{group}",
                    use_container_width=True,
                    hide_index=True,
                    column_order=shown_columns,
                    disabled=disabled_columns,
                )
                edited_rows.extend(_rows_from_editor(advanced_edited))

    if edited_rows:
        _set_workspace_rows(edited_rows)
        workspace = st.session_state["parameter_workspace"]

    validation = st.session_state.get("parameter_workspace_validation")
    col_validate, col_apply, col_run = st.columns(3)
    if col_validate.button("Validate parameters", type="primary", use_container_width=True):
        validation = validate_parameter_workspace(workspace["rows"])
        st.session_state["parameter_workspace_validation"] = validation
        if validation["isValid"]:
            st.success("Parameters are structurally safe to apply.")
        else:
            st.error("Some parameters need correction before they can be applied.")
    elif validation:
        if validation.get("isValid"):
            st.success("Parameters are structurally safe to apply.")
        else:
            st.error("Some parameters need correction before they can be applied.")

    can_apply = bool(validation and validation.get("isValid"))
    if col_apply.button("Apply parameters", disabled=not can_apply, use_container_width=True):
        try:
            cfg, econ = apply_parameter_workspace(
                st.session_state["config"],
                st.session_state["economics_config"],
                workspace["rows"],
            )
            cfg["workingDefaultPresetId"] = workspace.get("presetId")
            cfg["workingDefaultPresetVersion"] = workspace.get("presetVersion")
            econ.setdefault("metadata", {})["workingDefaultPresetId"] = workspace.get("presetId")
            econ["metadata"]["workingDefaultPresetVersion"] = workspace.get("presetVersion")
            st.session_state["config"] = cfg
            st.session_state["economics_config"] = econ
            st.session_state["parameter_workspace"] = build_parameter_workspace(cfg, econ)
            st.session_state["dirty_config"] = True
            st.session_state["dirty_economics"] = True
            st.session_state["results_stale"] = bool(st.session_state.get("results_bundle"))
            st.success("Parameters applied to the current analysis.")
        except Exception as exc:
            record_message("error", f"Could not apply parameters: {exc}")
            st.error("Could not apply the selected parameters.")
    if not can_apply and validation and validation.get("messages"):
        st.warning("Apply is disabled until the listed parameter errors are fixed.")
        arrow_safe_dataframe(validation["messages"], use_container_width=True)

    if col_run.button("Run analysis", disabled=not isinstance(st.session_state.get("config"), dict), use_container_width=True):
        st.page_link("pages/2_Run_Model.py", label="Open Run Analysis")


st.title("LTBI Screening Decision Tool")
st.write(
    "Compare the direct benefits, harms, resource requirements and health-system "
    "consequences of targeted LTBI screening and preventive treatment."
)
st.info(
    "This tool supports planning and sequencing decisions. It does not recommend "
    "denying care to any person or group."
)

col_default, col_review, col_results = st.columns(3)
with col_default:
    if st.button("Use default parameters", type="primary", use_container_width=True):
        try:
            _load_unified_defaults(show_workspace=False)
            st.success("APY / SA Health working defaults loaded.")
        except Exception as exc:
            record_message("error", f"Could not load working defaults: {exc}")
            st.error("Could not load the working defaults.")
with col_review:
    if st.button("Review or change parameters", use_container_width=True):
        try:
            if not isinstance(st.session_state.get("config"), dict):
                _load_unified_defaults(show_workspace=True)
            else:
                st.session_state.pop("recent_ltbi_run_route", None)
                st.session_state["parameter_workspace_visible"] = True
            st.success("Parameter workspace is ready.")
        except Exception as exc:
            record_message("error", f"Could not open parameter workspace: {exc}")
            st.error("Could not open the parameter workspace.")
with col_results:
    if st.session_state.get("results_bundle"):
        st.page_link("pages/3_Results.py", label="Continue to current results")
    else:
        st.button("Continue to current results", disabled=True, use_container_width=True)

config = st.session_state.get("config")
econ = st.session_state.get("economics_config")
if isinstance(config, dict) and isinstance(econ, dict):
    st.subheader("Current working defaults")
    st.caption("These are editable working defaults. Some APY-specific evidence inputs remain provisional.")
    arrow_safe_dataframe(parameter_summary(config, econ), use_container_width=True, hide_index=True)
    ltbi_state = resolve_ltbi_state_assumptions(config)
    unresolved_recent_ltbi = ltbi_state.get("baselineRecentLTBIProportion") is None
    if unresolved_recent_ltbi:
        st.subheader("Recent versus remote LTBI assumption remains provisional")
        st.write(
            "The proportion of baseline infections that were acquired relatively "
            "recently has not yet been established for this demonstration population. "
            "A default run therefore needs an explicit provisional route before "
            "analysis."
        )
        st.caption(
            "Run provisional working defaults uses the existing 0% compatibility "
            "assumption, temporarily representing all baseline infection as remote. "
            "Outputs remain provisional and evidence-review status is unchanged."
        )
        route_cols = st.columns(2)
        if route_cols[0].button("Run provisional working defaults", use_container_width=True):
            st.session_state["recent_ltbi_run_route"] = TECHNICAL_DEMONSTRATION_ROUTE
            st.success("Provisional working-default route selected.")
            st.page_link("pages/2_Run_Model.py", label="Continue to Run Analysis")
        route_cols[1].page_link("pages/4_Economics.py", label="Review this assumption")
    action_cols = st.columns(2)
    if not unresolved_recent_ltbi or st.session_state.get("recent_ltbi_run_route") == TECHNICAL_DEMONSTRATION_ROUTE:
        st.page_link("pages/2_Run_Model.py", label="Run with these defaults")
    else:
        st.info("Choose Run provisional working defaults or review this assumption before running.")
    if action_cols[1].button("Review or change parameters", key="show_workspace_secondary", use_container_width=True):
        st.session_state["parameter_workspace_visible"] = True
        st.rerun()

if st.session_state.get("parameter_workspace_visible"):
    _render_parameter_workspace()

with st.expander("Research and Development", expanded=False):
    st.write(
        "Technical implementation details, model version identifiers and diagnostics "
        "are available in the Research and Development section."
    )
    st.page_link("pages/8_Technical_Settings.py", label="Open Technical settings")
