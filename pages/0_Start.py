from __future__ import annotations

from copy import deepcopy
from typing import Any

import streamlit as st

from app.demographic_profile import (
    start_age_distribution_rows,
    start_risk_factor_rows,
)
from app.display import arrow_safe_dataframe
from app.parameter_workspace import (
    MODEL_METHOD_LABELS,
    PARAMETER_GROUPS,
    SIMULATION_MODE_LABELS,
    apply_parameter_workspace,
    build_parameter_workspace,
    changed_analysis_settings_count,
    changed_parameter_count,
    merge_parameter_display_edits,
    parameter_editor_rows,
    parameter_display_rows,
    unified_default_session_state,
    validate_parameter_workspace,
)
from app.state import (
    get_backend,
    init_session_state,
    mark_config_changed,
    mark_economics_changed,
    record_message,
    sanitize_reference_only_state,
    sync_backend_status,
)

def _page_link(path: str, *, label: str) -> None:
    try:
        st.page_link(path, label=label)
    except Exception:
        st.button(label, disabled=True)


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


def _restore_unified_defaults() -> None:
    previous_config = deepcopy(st.session_state.get("config"))
    previous_econ = deepcopy(st.session_state.get("economics_config"))
    state = unified_default_session_state()
    st.session_state["apy_backend_name"] = "python_apy"
    st.session_state["config"] = state["config"]
    st.session_state["economics_config"] = state["economics_config"]
    st.session_state["parameter_workspace"] = state["parameter_workspace"]
    st.session_state["working_default_preset"] = state["working_default_preset"]
    st.session_state["parameter_workspace_visible"] = True
    st.session_state["parameter_workspace_validation"] = None
    st.session_state.pop("recent_ltbi_run_route", None)
    st.session_state.pop("health_econ_workspace", None)
    st.session_state["setup_analysis_method"] = _analysis_mode_label(state["config"])
    st.session_state["setup_restore_defaults_message"] = "APY / SA Health defaults restored."
    _bump_parameter_editor_version()
    if previous_config != state["config"]:
        mark_config_changed()
    if previous_econ != state["economics_config"] and previous_config == state["config"]:
        mark_economics_changed()
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
    workspace["analysisSettingsChangedCount"] = changed_analysis_settings_count(workspace)
    st.session_state["parameter_workspace"] = workspace


def _bump_parameter_editor_version(group: str | None = None) -> None:
    st.session_state["parameter_editor_version"] = int(st.session_state.get("parameter_editor_version", 0)) + 1
    if group:
        st.session_state[f"parameter_editor_version_{group}"] = int(
            st.session_state.get(f"parameter_editor_version_{group}", 0)
        ) + 1


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


def _sync_workspace_after_direct_config_change() -> None:
    st.session_state["parameter_workspace"] = build_parameter_workspace(
        st.session_state["config"],
        st.session_state.get("economics_config") or {},
    )
    st.session_state["parameter_workspace_validation"] = None
    mark_config_changed()


def _analysis_mode_label(config: dict[str, Any]) -> str:
    return MODEL_METHOD_LABELS.get(
        str(config.get("analysisMethod") or "expected_value"),
        "Quick deterministic preview - single expected-value calculation",
    )


def _simulation_mode_label(config: dict[str, Any]) -> str:
    return SIMULATION_MODE_LABELS.get(str(config.get("simulationMode") or "standard"), "SA Health reference: 2,000 repetitions")


def _apply_simulation_mode(config: dict[str, Any], mode_label: str) -> None:
    code = next((key for key, label in SIMULATION_MODE_LABELS.items() if label == mode_label), "custom")
    config["simulationMode"] = code
    config["simulationModeLabel"] = mode_label
    if code == "quick_preview":
        config["nReps"] = 100
    elif code == "intermediate":
        config["nReps"] = 500
    elif code == "standard":
        config["nReps"] = 2000


def _run_classification(config: dict[str, Any]) -> str:
    if str(config.get("analysisMethod")) == "expected_value":
        return "Quick deterministic preview"
    reps = int(float(config.get("nReps") or 0))
    seed = int(float(config.get("seed") or 0))
    if reps == 2000 and seed == 1:
        return "SA Health reference run"
    if reps < 2000:
        return "Exploratory preview run"
    return "Custom stochastic run"


def _render_epidemiological_basis() -> None:
    st.subheader("Epidemiological basis")
    st.success("Future TB outcomes use the assumptions applied in the SA Health report.")
    st.caption(
        "Incidence, infection-history, calibration and progression assumptions are fixed "
        "for this SA Health version and documented in Evidence & Assumptions."
    )


def _render_analysis_settings_controls(config: dict[str, Any]) -> None:
    st.subheader("Analysis settings")
    before = deepcopy(config)
    method_options = list(MODEL_METHOD_LABELS.values())
    widget_key = "setup_analysis_method"
    configured_label = _analysis_mode_label(config)
    if st.session_state.get(widget_key) not in method_options:
        st.session_state[widget_key] = configured_label
    method_label = st.radio(
        "Analysis type",
        method_options,
        index=None,
        horizontal=True,
        key=widget_key,
    )
    if method_label not in method_options:
        method_label = configured_label
    config["analysisMethod"] = next(
        (code for code, label in MODEL_METHOD_LABELS.items() if label == method_label),
        "expected_value",
    )
    config["analysisMethodLabel"] = method_label
    if config["analysisMethod"] == "agent_based":
        mode_options = list(SIMULATION_MODE_LABELS.values())
        mode_label = st.selectbox(
            "Repetitions",
            mode_options,
            index=mode_options.index(_simulation_mode_label(config)) if _simulation_mode_label(config) in mode_options else mode_options.index("Custom"),
        )
        _apply_simulation_mode(config, mode_label)
        if config.get("simulationMode") == "custom":
            config["nReps"] = st.number_input(
                "Custom repetitions",
                min_value=1,
                max_value=5000,
                value=int(float(config.get("nReps") or 2000)),
                step=1,
                help="Use a positive whole number. Very large runs are not suitable for interactive review.",
            )
        seed_value = st.number_input(
            "Random seed",
            min_value=0,
            max_value=2_147_483_647,
            value=int(float(config.get("seed") or 1)),
            step=1,
        )
        config["seed"] = int(seed_value)
        if int(config.get("nReps") or 0) < 2000:
            st.warning("Preview analyses run faster but do not reproduce the SA Health reference.")
        st.info(f"{_run_classification(config)}: {int(config.get('nReps')):,} repetitions, seed {int(config.get('seed'))}.")
    else:
        st.info("Quick deterministic preview: no stochastic repetitions or random seed are used.")
    if config != before:
        st.session_state["config"] = config
        _sync_workspace_after_direct_config_change()


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
    if not rows:
        return []
    version = int(st.session_state.get("parameter_editor_version", 0))
    group_version = int(st.session_state.get(f"parameter_editor_version_{key.rsplit('_', 1)[-1]}", 0))
    return _rows_from_editor(
        st.data_editor(
            parameter_editor_rows(rows),
            key=f"{key}_{version}_{group_version}",
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
    st.subheader("Parameters")
    st.caption(
        "Edit the Value used by model column to create a user-defined override. "
        "Risk-factor prevalence values are shown as percentages; entering 43% is applied as 0.43."
    )
    scientific_changes = changed_parameter_count(workspace)
    analysis_changes = changed_analysis_settings_count(workspace)
    st.metric("Scientific parameter overrides", scientific_changes)
    if analysis_changes:
        st.info("Analysis settings changed. Rerun analysis before interpreting previous results.")

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

    visible_by_group: dict[str, list[dict[str, Any]]] = {}
    for group in PARAMETER_GROUPS:
        group_rows = [row for row in workspace["rows"] if row.get("group") == group]
        visible_rows = [
            row
            for row in group_rows
            if not str(row.get("parameterId") or "").startswith("demography.age.")
            and row.get("operationalStatus") != "descriptive_metadata"
            and not row.get("advanced")
        ]
        if visible_rows:
            visible_by_group[group] = visible_rows
    visible_groups = list(visible_by_group)
    tabs = st.tabs(visible_groups)
    for tab, group in zip(tabs, visible_groups):
        with tab:
            visible_rows = visible_by_group[group]
            st.caption("Rows where Source is User-defined contain user-entered overrides.")
            edited = _editable_parameter_table(visible_rows, key=f"parameter_editor_{group}")
            edited_rows.extend(merge_parameter_display_edits(visible_rows, edited))

    if edited_rows:
        edited_by_id = {row.get("parameterId"): row for row in edited_rows}
        complete_rows = [edited_by_id.get(row.get("parameterId"), row) for row in workspace["rows"]]
        previous_rows = deepcopy(workspace["rows"])
        _set_workspace_rows(complete_rows)
        workspace = st.session_state["parameter_workspace"]
        if _rows_changed(previous_rows, workspace["rows"]):
            validation = validate_parameter_workspace(workspace["rows"])
            st.session_state["parameter_workspace_validation"] = validation
            if validation["isValid"]:
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
                    st.rerun()
                except Exception as exc:
                    record_message("error", f"Could not apply parameters: {exc}")
                    st.error("Could not apply the selected parameters.")
            else:
                st.rerun()

    validation = st.session_state.get("parameter_workspace_validation")
    if validation and not validation.get("isValid"):
        st.error("Some parameters need correction before they can be used.")
        st.dataframe(arrow_safe_dataframe(validation["messages"]), use_container_width=True, hide_index=True)


def _render_age_risk_summary(config: dict[str, Any]) -> None:
    st.caption(
        "Repository APY demographic default; external provenance not independently reviewed in this workflow. "
        "These are the resolved values currently used by the model."
    )
    age_rows = [
        {
            "Age group": row["Age group"],
            "Proportion used by model": row["Current proportion used by model"],
        }
        for row in start_age_distribution_rows(config)
    ]
    risk_rows = start_risk_factor_rows(config)
    with st.expander("View age distribution", expanded=False):
        st.dataframe(
            arrow_safe_dataframe(age_rows),
            use_container_width=True,
            hide_index=True,
        )
    with st.expander("View risk factors", expanded=False):
        st.dataframe(
            arrow_safe_dataframe(risk_rows),
            use_container_width=True,
            hide_index=True,
        )


init_session_state()
if not isinstance(st.session_state.get("config"), dict) or not isinstance(st.session_state.get("economics_config"), dict):
    _load_unified_defaults(show_workspace=True)
sanitize_reference_only_state()

st.title("Set up")
restore_message = st.session_state.pop("setup_restore_defaults_message", "")
if restore_message:
    st.success(restore_message)
st.write(
    "Define the population, screening strategy and model assumptions for the "
    "LTBI Screening Decision Tool. Changes to population, testing, treatment "
    "or targeting inputs require a new analysis run."
)

config = st.session_state.get("config")
econ = st.session_state.get("economics_config")
if isinstance(config, dict) and isinstance(econ, dict):
    notice = st.session_state.pop("reference_only_migration_notice", "")
    if notice:
        st.warning(notice)
    _render_epidemiological_basis()
    _render_analysis_settings_controls(config)
    _render_parameter_workspace()
    _render_age_risk_summary(st.session_state["config"])

st.button("Restore APY defaults", use_container_width=True, on_click=_restore_unified_defaults)

if isinstance(st.session_state.get("config"), dict):
    effective_config = st.session_state["config"]
    if str(effective_config.get("analysisMethod")) == "expected_value":
        st.success("Run Analysis will use: Quick deterministic preview - single expected-value calculation.")
    else:
        st.success(
            "Run Analysis will use: SA Health report analysis - "
            f"{int(effective_config.get('nReps') or 0):,} runs, seed {int(effective_config.get('seed') or 1)}."
        )
    _page_link("pages/2_Run_Model.py", label="Open Run Analysis")
