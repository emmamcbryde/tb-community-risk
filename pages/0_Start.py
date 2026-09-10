from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any

import streamlit as st

from app.demographic_profile import (
    restore_apy_demographic_defaults,
    risk_factor_rows,
    start_age_distribution_rows,
)
from app.display import arrow_safe_dataframe
from app.parameter_workspace import (
    PARAMETER_GROUPS,
    apply_parameter_workspace,
    build_parameter_workspace,
    changed_parameter_count,
    merge_parameter_display_edits,
    parameter_editor_rows,
    parameter_summary,
    parameter_display_rows,
    reset_all_parameters,
    reset_parameter_group,
    unified_default_session_state,
    validate_parameter_workspace,
)
from app.run_analysis_controls import TECHNICAL_DEMONSTRATION_ROUTE
from app.state import (
    get_backend,
    init_session_state,
    mark_config_changed,
    mark_economics_changed,
    mark_validation_completed,
    record_message,
    sync_backend_status,
)
from adapters.paths import scenarios_dir
from engine.apy.ltbi_state import resolve_ltbi_state_assumptions


init_session_state()


def scenario_path(filename: str) -> Path:
    base = scenarios_dir()
    base.mkdir(parents=True, exist_ok=True)
    name = filename.strip() or "streamlit_analysis.json"
    path = Path(name)
    if not path.is_absolute():
        path = base / path
    return path


def validation_rows(report: dict[str, Any]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for group in ("errors", "warnings", "infos"):
        issues = report.get(group) or []
        if isinstance(issues, dict):
            issues = [issues]
        for issue in issues:
            if isinstance(issue, dict):
                rows.append(
                    {
                        "Severity": str(issue.get("severity", group[:-1])),
                        "Field": str(issue.get("fieldLabel", issue.get("field", ""))),
                        "Message": str(issue.get("message", "")),
                    }
                )
    return rows


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
        current = row.get("currentValue")
        row["changedFromDefault"] = bool(row.get("isUserOverride")) or str(current) != str(defaults.get(row.get("parameterId")))
        if row.get("sourceObject") in {"ageDistributionBand", "config"} and str(row.get("parameterId", "")).startswith(
            ("demography.age.", "demography.risk.")
        ):
            row["valueUsedByModel"] = row.get("sourceDefaultValue") if current in (None, "") else current
            row["effectiveSource"] = "User-defined" if row.get("isUserOverride") else "Repository APY default"
    workspace["rows"] = rows
    workspace["changedCount"] = changed_parameter_count(workspace)
    st.session_state["parameter_workspace"] = workspace


def _rows_changed(before: list[dict[str, Any]], after: list[dict[str, Any]]) -> bool:
    before_by_id = {row.get("parameterId"): row for row in before}
    after_by_id = {row.get("parameterId"): row for row in after}
    if set(before_by_id) != set(after_by_id):
        return True
    for parameter_id, before_row in before_by_id.items():
        after_row = after_by_id[parameter_id]
        for key in ("currentValue", "valueUsedByModel", "effectiveSource", "isUserOverride"):
            if str(before_row.get(key)) != str(after_row.get(key)):
                return True
    return False


def _workspace_change_scope(before: list[dict[str, Any]], after: list[dict[str, Any]]) -> tuple[bool, bool]:
    before_by_id = {row.get("parameterId"): row for row in before}
    model_changed = False
    economics_changed = False
    for row in after:
        previous = before_by_id.get(row.get("parameterId"), {})
        changed = any(
            str(previous.get(key)) != str(row.get(key))
            for key in ("currentValue", "valueUsedByModel", "effectiveSource", "isUserOverride")
        )
        if not changed:
            continue
        status = row.get("operationalStatus")
        if status == "authoritative_model_input":
            model_changed = True
        elif status == "authoritative_economic_input":
            economics_changed = True
    return model_changed, economics_changed


def _editable_parameter_table(rows: list[dict[str, Any]], *, key: str) -> list[dict[str, Any]]:
    return _rows_from_editor(
        st.data_editor(
            parameter_editor_rows(rows),
            key=key,
            use_container_width=True,
            hide_index=True,
            column_order=["Parameter", "Value used by model", "Unit", "Source"],
            disabled=["Parameter", "Unit", "Source"],
            column_config={
                "Value used by model": st.column_config.TextColumn(
                    "Value used by model",
                    help="Edit this cell to create a user-defined override. Use Restore defaults to return to the source value.",
                )
            },
        )
    )


def _render_parameter_workspace() -> None:
    workspace = _current_workspace()
    if not workspace:
        return
    st.subheader("Review or change parameters")
    st.caption(
        "These grouped inputs update the same configuration used by the analysis. "
        "Changing a value does not change its evidence source or review status. "
        "Blank demographic or risk-factor override fields mean use source defaults shown above; "
        "they are not missing model inputs."
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
    override_rows = [
        {
            "Parameter": row.get("label"),
            "Value used by model": row.get("valueUsedByModel"),
            "Source": "User-defined",
        }
        for row in workspace.get("rows", [])
        if row.get("isUserOverride")
    ]
    if override_rows:
        st.warning("The parameters below include user-defined values.")
        st.dataframe(
            arrow_safe_dataframe(override_rows),
            use_container_width=True,
            hide_index=True,
        )

    visible_groups = PARAMETER_GROUPS
    tabs = st.tabs(visible_groups)
    for tab, group in zip(tabs, visible_groups):
        with tab:
            group_rows = [row for row in workspace["rows"] if row.get("group") == group]
            if st.button(f"Reset this section: {group}", key=f"reset_{group}"):
                _set_workspace_rows(reset_parameter_group(workspace["rows"], group))
                st.session_state["parameter_workspace_validation"] = None
                st.rerun()
            standard_rows = [row for row in group_rows if not row.get("advanced")]
            advanced_rows = [row for row in group_rows if row.get("advanced")]
            st.caption("Rows where Source is User-defined contain user-entered overrides.")
            read_only_rows = [row for row in standard_rows if row.get("editableType") == "read_only"]
            editable_rows = [row for row in standard_rows if row.get("editableType") != "read_only"]
            if read_only_rows:
                st.dataframe(
                    arrow_safe_dataframe(parameter_display_rows(read_only_rows)),
                    use_container_width=True,
                    hide_index=True,
                )
            edited = _editable_parameter_table(editable_rows, key=f"parameter_editor_{group}")
            edited_rows.extend(merge_parameter_display_edits(editable_rows, edited))
            edited_rows.extend(read_only_rows)
            with st.expander("Advanced"):
                advanced_read_only = [row for row in advanced_rows if row.get("editableType") == "read_only"]
                advanced_editable = [row for row in advanced_rows if row.get("editableType") != "read_only"]
                if advanced_read_only:
                    st.dataframe(
                        arrow_safe_dataframe(parameter_display_rows(advanced_read_only)),
                        use_container_width=True,
                        hide_index=True,
                    )
                advanced_edited = _editable_parameter_table(
                    advanced_editable,
                    key=f"parameter_advanced_editor_{group}",
                )
                edited_rows.extend(merge_parameter_display_edits(advanced_editable, advanced_edited))
                edited_rows.extend(advanced_read_only)

    if edited_rows:
        previous_rows = deepcopy(workspace["rows"])
        _set_workspace_rows(edited_rows)
        workspace = st.session_state["parameter_workspace"]
        if _rows_changed(previous_rows, workspace["rows"]):
            st.session_state["parameter_workspace_validation"] = None
            st.rerun()

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
            applied_workspace = build_parameter_workspace(
                st.session_state["config"],
                st.session_state["economics_config"],
            )
            model_changed, economics_changed = _workspace_change_scope(applied_workspace["rows"], workspace["rows"])
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
            if model_changed:
                mark_config_changed()
            if economics_changed and not model_changed:
                mark_economics_changed()
            st.success("Parameters applied to the current analysis.")
        except Exception as exc:
            record_message("error", f"Could not apply parameters: {exc}")
            st.error("Could not apply the selected parameters.")
    if not can_apply and validation and validation.get("messages"):
        st.warning("Apply is disabled until the listed parameter errors are fixed.")
        st.dataframe(
            arrow_safe_dataframe(validation["messages"]),
            use_container_width=True,
            hide_index=True,
        )

    if col_run.button("Run analysis", disabled=not isinstance(st.session_state.get("config"), dict), use_container_width=True):
        st.page_link("pages/2_Run_Model.py", label="Open Run Analysis")


def _render_age_risk_summary(config: dict[str, Any]) -> None:
    with st.expander("Age distribution and risk factors", expanded=False):
        st.caption(
            "Repository APY demographic default; external provenance not independently reviewed in this workflow. "
            "These are the resolved values currently used by the model."
        )
        st.dataframe(
            arrow_safe_dataframe(start_age_distribution_rows(config)),
            use_container_width=True,
            hide_index=True,
        )
        st.dataframe(
            arrow_safe_dataframe(risk_factor_rows(config)),
            use_container_width=True,
            hide_index=True,
        )
        if st.button("Restore APY demographic defaults"):
            restored = restore_apy_demographic_defaults(config)
            if restored != config:
                st.session_state["config"] = restored
                st.session_state.pop("recent_ltbi_run_route", None)
                st.session_state["parameter_workspace"] = build_parameter_workspace(
                    restored,
                    st.session_state.get("economics_config") or {},
                )
                mark_config_changed()
                st.success("APY demographic defaults restored. Rerun epidemiology before using previous results.")
                st.rerun()
            st.info("APY demographic defaults are already loaded.")


def _render_configuration_status() -> None:
    st.subheader("Validation and next action")
    if st.session_state.get("dirty_config"):
        st.warning("Configuration changes require a new analysis run before previous results are used.")
    elif st.session_state.get("results_bundle"):
        st.success("Current analysis results match the saved configuration.")
    else:
        st.info("No analysis has been run in this session.")
    if st.session_state.get("results_stale"):
        st.warning("Existing results are stale because setup inputs changed after the last run.")
    if st.session_state.get("dirty_economics") and not st.session_state.get("results_stale"):
        st.info("Economic inputs changed. Health economics can be recalculated without rerunning screening outcomes.")
    changed = changed_parameter_count(st.session_state.get("parameter_workspace") or {})
    if changed:
        st.info(f"{changed} parameters differ from the SA Health working reference.")
    else:
        st.success("No parameter overrides are currently applied.")

    cols = st.columns(2)
    if cols[0].button("Validate setup", use_container_width=True):
        try:
            report = get_backend().validate_config(st.session_state["config"])
            st.session_state["validation_report"] = report
            mark_validation_completed()
            if report.get("isValid") is True:
                st.success("Setup is valid.")
            elif report.get("isValid") is False:
                st.error("Setup has validation errors.")
            else:
                st.info("Validation completed.")
        except Exception as exc:
            record_message("error", f"Validation failed: {exc}")
            st.error("Validation failed.")
    cols[1].page_link("pages/2_Run_Model.py", label="Proceed to Run Analysis")

    report = st.session_state.get("validation_report")
    if report:
        rows = validation_rows(report)
        if rows:
            st.dataframe(arrow_safe_dataframe(rows), use_container_width=True, hide_index=True)


def _render_analysis_file_controls() -> None:
    with st.expander("Save or load setup", expanded=False):
        st.caption("Use these controls to save or reload the current setup and economic assumptions.")
        file_name = st.text_input("Analysis file", value="streamlit_analysis.json")
        cols = st.columns(2)
        if cols[0].button("Save setup"):
            try:
                path = scenario_path(file_name)
                st.session_state["save_info"] = get_backend().save_scenario(
                    st.session_state["config"],
                    str(path),
                    st.session_state.get("economics_config"),
                )
                sync_backend_status(get_backend().status())
                st.success(f"Saved setup to {path}")
            except Exception as exc:
                record_message("error", f"Save failed: {exc}")
                st.error("Save failed.")

        if cols[1].button("Load setup"):
            try:
                path = scenario_path(file_name)
                payload = json.loads(path.read_text(encoding="utf-8"))
                loaded_config = payload.get("config")
                if not isinstance(loaded_config, dict):
                    raise ValueError("Analysis JSON does not contain a configuration object.")
                loaded_econ = payload.get("economics")
                if not isinstance(loaded_econ, dict):
                    loaded_econ = st.session_state.get("economics_config") or {}
                st.session_state["config"] = loaded_config
                st.session_state["economics_config"] = loaded_econ
                st.session_state["parameter_workspace"] = build_parameter_workspace(loaded_config, loaded_econ)
                st.session_state["parameter_workspace_visible"] = True
                st.session_state.pop("recent_ltbi_run_route", None)
                st.session_state["validation_report"] = None
                st.session_state["load_info"] = {
                    "filename": str(path),
                    "contractVersion": payload.get("contractVersion", ""),
                    "scenarioLabel": payload.get("scenarioLabel", ""),
                }
                st.session_state["results_bundle"] = None
                st.session_state["results_stale"] = False
                st.session_state["dirty_config"] = True
                sync_backend_status(get_backend().status())
                st.success(f"Loaded setup from {path}")
                st.rerun()
            except Exception as exc:
                record_message("error", f"Load failed: {exc}")
                st.error("Load failed.")


st.title("Set up")
st.write(
    "Define the population, screening strategy and model assumptions for the "
    "LTBI Screening Decision Tool. Changes to population, testing, treatment "
    "or targeting inputs require a new analysis run."
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
    st.dataframe(
        arrow_safe_dataframe(parameter_summary(config, econ)),
        use_container_width=True,
        hide_index=True,
    )
    _render_age_risk_summary(config)
    st.subheader("Strategy controls")
    st.caption(
        "Use Review or change parameters below for the single editable setup table. "
        "The same values are used by Run Analysis; there is no second strategy editor."
    )
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
        route_cols[1].page_link("pages/6_Evidence_Assumptions.py", label="Review this assumption")
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

if isinstance(st.session_state.get("config"), dict):
    _render_configuration_status()
    _render_analysis_file_controls()

with st.expander("Research and Development", expanded=False):
    st.write(
        "Technical implementation details, model version identifiers and diagnostics "
        "are available in the Research and Development section."
    )
    st.page_link("pages/8_Technical_Settings.py", label="Open Technical settings")
