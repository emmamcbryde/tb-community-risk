from __future__ import annotations

from copy import deepcopy
from typing import Any

import streamlit as st

from app.demographic_profile import (
    restore_apy_demographic_defaults,
    risk_factor_rows,
    start_age_distribution_rows,
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
    parameter_summary,
    parameter_display_rows,
    reset_all_parameters,
    reset_parameter_group,
    unified_default_session_state,
    validate_parameter_workspace,
)
from app.state import (
    get_backend,
    init_session_state,
    mark_config_changed,
    mark_economics_changed,
    record_message,
    sync_backend_status,
)
from engine.apy.config import normalise_config
from engine.apy.data import load_parameters_from_config
from engine.apy.infection_history import (
    EXPERIMENTAL_STATUS_LABEL,
    TRAJECTORY_PRESETS,
    calibrate_infection_history,
    configure_compatibility_reference_assumptions,
    configure_infection_history_assumptions,
    has_explicit_infection_history,
    infection_history_readiness_status,
    is_experimental_infection_history_results,
)


init_session_state()


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
        "Expected outcomes — single deterministic run",
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
        return "Expected-outcomes run"
    reps = int(float(config.get("nReps") or 0))
    seed = int(float(config.get("seed") or 0))
    if reps == 2000 and seed == 1:
        return "SA Health reference run"
    if reps < 2000:
        return "Exploratory preview run"
    return "Custom stochastic run"


def _format_percent(value: Any) -> str:
    try:
        return f"{float(value) * 100:.2f}%"
    except (TypeError, ValueError):
        return ""


def _decimal(value: Any, digits: int = 2) -> str:
    try:
        return f"{float(value):,.{digits}f}"
    except (TypeError, ValueError):
        return "Unavailable"


def _target_ltbi_prevalence(config: dict[str, Any]) -> float:
    value = config.get("ltbiPrevalence")
    return 47 / 624 if value in (None, "", []) else float(value)


def _render_infection_history_controls(config: dict[str, Any]) -> None:
    current_explicit = has_explicit_infection_history(config)
    st.subheader("Epidemiological basis")
    if current_explicit:
        st.warning(
            f"{EXPERIMENTAL_STATUS_LABEL}. Restore the SA Health report reference "
            "before reproducing report-facing epidemiology or economics."
        )
    else:
        st.success("SA Health report reference: validated software-compatible stochastic anchor.")
        st.caption(
            "The reference preserves the validated software-compatible epidemiological "
            "anchor. Its implicit early phase is not a measured recent-infection fraction."
        )

    with st.expander("Experimental infection-history analysis", expanded=current_explicit):
        st.warning(
            "This exploratory analysis reconstructs the timing of past infection from "
            "prevalence, age pattern and an assumed historical trajectory. It is not "
            "the SA Health report reference. Its current TB progression calibration is "
            "not validated for decision-making or health-economic reporting."
        )
        if not current_explicit:
            if st.button("Enable experimental infection-history analysis"):
                st.session_state["config"] = configure_infection_history_assumptions(config, "steady")
                st.session_state.pop("recent_ltbi_run_route", None)
                _sync_workspace_after_direct_config_change()
                st.rerun()
            return

        if st.button("Restore SA Health report reference"):
            st.session_state["config"] = configure_compatibility_reference_assumptions(config)
            st.session_state.pop("recent_ltbi_run_route", None)
            _sync_workspace_after_direct_config_change()
            st.rerun()

        labels = [preset["label"] for preset in TRAJECTORY_PRESETS.values()]
        label_to_key = {preset["label"]: key for key, preset in TRAJECTORY_PRESETS.items()}
        nested = config.get("ltbiStateAssumptions") or {}
        current_key = str(nested.get("infectionPressureTrajectory") or "steady").lower()
        current_label = TRAJECTORY_PRESETS.get(current_key, TRAJECTORY_PRESETS["steady"])["label"]
        trajectory_label = st.selectbox(
            "Historical TB infection pressure",
            labels,
            index=labels.index(current_label) if current_label in labels else labels.index("Steady"),
            help=(
                "Rising means infection pressure has increased toward the present; "
                "falling means it has declined toward the present."
            ),
        )
        trajectory = label_to_key[trajectory_label]
        configured = configure_infection_history_assumptions(config, trajectory)
        if configured != config:
            st.session_state["config"] = configured
            st.session_state.pop("recent_ltbi_run_route", None)
            _sync_workspace_after_direct_config_change()
            st.rerun()

        try:
            cfg = normalise_config(configured)
            pars = load_parameters_from_config(cfg)
            nested = cfg.get("ltbiStateAssumptions") or {}
            derived = calibrate_infection_history(
                pars,
                target_prevalence=_target_ltbi_prevalence(cfg),
                target_age_or=float(cfg.get("targetAgeOR") or 7.54),
                trajectory=trajectory,
                trend_rate_per_year=nested.get("infectionPressureTrendRatePerYear"),
                recent_window_years=nested.get("recentDefinitionYears") or 2.0,
            )
            diagnostic_rows = [
                {"Quantity": "Status", "Value": EXPERIMENTAL_STATUS_LABEL},
                {"Quantity": "Assumed trajectory", "Value": TRAJECTORY_PRESETS[trajectory]["label"]},
                {"Quantity": "Baseline LTBI prevalence", "Value": _format_percent(derived["expectedPrevalence"])},
                {
                    "Quantity": "Infection acquired within preceding two years",
                    "Value": _format_percent(derived["recentFraction"]),
                },
                {
                    "Quantity": "Infection acquired earlier than two years",
                    "Value": _format_percent(derived["remoteFraction"]),
                },
            ]
            results_bundle = st.session_state.get("results_bundle")
            if is_experimental_infection_history_results(results_bundle):
                totals = (((results_bundle or {}).get("technical") or {}).get("eventLedger") or {}).get("replicateTotals")
                try:
                    import pandas as pd

                    frame = totals if isinstance(totals, pd.DataFrame) else pd.DataFrame(totals or [])
                    comparator = _event_value(frame, "comparator", "active_tb_cases")
                    averted = _event_value(frame, "intervention", "active_tb_cases_prevented")
                    diagnostic_rows.extend(
                        [
                            {"Quantity": "Comparator active TB in current run", "Value": _decimal(comparator, 2)},
                            {"Quantity": "Active TB averted in current run", "Value": _decimal(averted, 2)},
                            {
                                "Quantity": "Difference from SA Health reference comparator",
                                "Value": _decimal(None if comparator is None else comparator - 37.305, 2),
                            },
                            {
                                "Quantity": "Difference from SA Health reference averted",
                                "Value": _decimal(None if averted is None else averted - 11.9665, 2),
                            },
                        ]
                    )
                except Exception:
                    pass
            else:
                diagnostic_rows.extend(
                    [
                        {"Quantity": "Comparator active TB in current run", "Value": "Run analysis to calculate"},
                        {"Quantity": "Active TB averted in current run", "Value": "Run analysis to calculate"},
                    ]
                )
            st.dataframe(
                arrow_safe_dataframe(diagnostic_rows),
                use_container_width=True,
                hide_index=True,
            )
            st.markdown("Recent fraction by age group")
            st.dataframe(
                arrow_safe_dataframe(
                    [
                        {
                            "Age group": row["ageGroup"],
                            "LTBI prevalence": _format_percent(row["ltbiPrevalence"]),
                            "Infection acquired within preceding two years": _format_percent(
                                row["recentFractionAmongPrevalent"]
                            ),
                        }
                        for row in derived["ageRows"]
                    ]
                ),
                use_container_width=True,
                hide_index=True,
            )
            st.caption(
                "The two-year infection-timing fraction is model-derived from the selected "
                "transmission-history assumption. It is not directly observed and is not "
                "equivalent to the early higher-progression-risk state, which currently has "
                "a five-year mean residence time."
            )
            readiness = infection_history_readiness_status()
            st.markdown("Experimental readiness guard")
            st.dataframe(arrow_safe_dataframe(readiness["items"]), use_container_width=True, hide_index=True)
        except Exception as exc:
            st.error(f"Could not derive recent-versus-remote LTBI from this setup: {exc}")


def _event_value(frame, arm: str, event: str) -> float | None:
    if frame is None or getattr(frame, "empty", True):
        return None
    subset = frame[
        frame["arm"].astype(str).eq(arm)
        & frame["eventName"].astype(str).eq(event)
    ]
    if subset.empty:
        return None
    return float(subset["value"].astype(float).mean())


def _render_analysis_settings_controls(config: dict[str, Any]) -> None:
    st.subheader("Analysis settings")
    st.caption("Choose faster preview runs for exploration. The SA Health reference uses 2,000 repetitions and seed 1.")
    before = deepcopy(config)
    method_options = list(MODEL_METHOD_LABELS.values())
    method_label = st.radio(
        "Analysis type",
        method_options,
        index=method_options.index(_analysis_mode_label(config)) if _analysis_mode_label(config) in method_options else 0,
        horizontal=True,
    )
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
        st.info("Expected outcomes run: no stochastic repetitions are used.")
    if config != before:
        st.session_state["config"] = config
        _sync_workspace_after_direct_config_change()
        st.rerun()


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
    st.subheader("Review or change parameters")
    st.caption(
        "These grouped inputs update the same configuration used by the analysis. "
        "Changing a value does not change its evidence source or review status. "
        "Blank demographic or risk-factor override fields mean use source defaults shown above; "
        "they are not missing model inputs."
    )
    scientific_changes = changed_parameter_count(workspace)
    analysis_changes = changed_analysis_settings_count(workspace)
    st.metric("Scientific parameter overrides", scientific_changes)
    if analysis_changes:
        st.info("Analysis settings changed. Rerun analysis before interpreting previous results.")
    top_cols = st.columns([1, 1, 3])
    if top_cols[0].button("Use all defaults", use_container_width=True):
        _set_workspace_rows(reset_all_parameters(workspace["rows"]))
        _bump_parameter_editor_version()
        st.session_state["parameter_workspace_validation"] = None
        st.success("Defaults restored.")
        st.stop()
    if top_cols[1].button("Reset all changes", use_container_width=True):
        _set_workspace_rows(reset_all_parameters(workspace["rows"]))
        _bump_parameter_editor_version()
        st.session_state["parameter_workspace_validation"] = None
        st.success("Defaults restored.")
        st.stop()

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
                _bump_parameter_editor_version(group)
                st.session_state["parameter_workspace_validation"] = None
                st.success(f"{group} defaults restored.")
                st.stop()
            standard_rows = [row for row in group_rows if not row.get("advanced")]
            advanced_rows = [
                row
                for row in group_rows
                if row.get("advanced")
                and not str(row.get("parameterId") or "").startswith("demography.age.")
            ]
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
    col_validate, col_apply = st.columns(2)
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
    risk_rows = [
        {
            "Risk factor": row["Risk factor"],
            "Proportion used by model": _format_percent(row.get("Prevalence")),
        }
        for row in risk_factor_rows(config)
    ]
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
        _page_link("pages/3_Results.py", label="Continue to current results")
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
    _render_infection_history_controls(config)
    _render_analysis_settings_controls(config)

if st.session_state.get("parameter_workspace_visible"):
    _render_parameter_workspace()

if isinstance(st.session_state.get("config"), dict):
    _page_link("pages/2_Run_Model.py", label="Open Run Analysis")
