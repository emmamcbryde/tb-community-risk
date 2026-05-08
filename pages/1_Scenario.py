from __future__ import annotations

import json
from pathlib import Path

import streamlit as st

from app.display import arrow_safe_dataframe, display_string
from app.state import (
    get_backend,
    init_session_state,
    mark_config_changed,
    mark_validation_completed,
    record_message,
    sync_backend_status,
)
from adapters.paths import scenarios_dir


init_session_state()
backend = get_backend()

st.title("Scenario")
st.caption("Initial APY v9 scenario shell backed by MATLAB.")


def display_backend_status() -> None:
    status = backend.status()
    sync_backend_status(status)
    cols = st.columns(3)
    cols[0].metric("Backend", status.get("name", "unknown"))
    cols[1].metric("MATLAB started", "yes" if status.get("started") else "no")
    cols[2].metric("Adapter", "error" if status.get("error") else "ready")
    st.caption(f"ABM path: {status.get('abm_path', '')}")
    if status.get("error"):
        st.error(status["error"])


def config_overview_rows(config: dict) -> list[dict[str, object]]:
    fields = [
        "scenarioLabel",
        "configVersion",
        "modelVersion",
        "N",
        "nReps",
        "seed",
        "screenCoverage",
        "screenWindow",
        "followHorizon",
        "screeningStrategy",
        "testType",
        "regimen",
        "outputDir",
    ]
    return [
        {"field": field, "value": config.get(field)}
        for field in fields
        if field in config
    ]


def validation_rows(report: dict) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for group in ("errors", "warnings", "infos"):
        issues = report.get(group) or []
        if isinstance(issues, dict):
            issues = [issues]
        for issue in issues:
            if isinstance(issue, dict):
                rows.append(
                    {
                        "severity": display_string(issue.get("severity", group[:-1])),
                        "field": display_string(issue.get("field", "")),
                        "code": display_string(issue.get("code", "")),
                        "message": display_string(issue.get("message", "")),
                    }
                )
    return rows


def choice_index(value: object, options: list[str], fallback: int = 0) -> int:
    try:
        return options.index(str(value))
    except ValueError:
        return fallback


def reset_run_state() -> None:
    st.session_state["validation_report"] = None
    st.session_state["results_bundle"] = None
    st.session_state["results_stale"] = False
    st.session_state["dirty_config"] = False


def scenario_path(filename: str) -> Path:
    base = scenarios_dir()
    base.mkdir(parents=True, exist_ok=True)
    name = filename.strip() or "streamlit_scenario.json"
    path = Path(name)
    if not path.is_absolute():
        path = base / path
    return path


st.subheader("Backend Status")
display_backend_status()

if st.button("Load MATLAB defaults", type="primary"):
    try:
        st.session_state["config"] = backend.default_config()
        reset_run_state()
        sync_backend_status(backend.status())
        st.success("Default scenario loaded.")
    except Exception as exc:
        message = f"Could not load defaults: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

config = st.session_state.get("config")
if config:
    if st.session_state.get("dirty_config"):
        st.warning("Config has unvalidated changes.")
    if st.session_state.get("results_stale"):
        st.warning("Existing results are stale because config changed after the last run.")

    st.subheader("Edit Scenario")
    test_options = ["IGRA", "TST"]
    regimen_options = ["3HP", "4R", "3HR", "6H", "9H"]
    strategy_options = ["random", "ltbi", "cure", "prevent"]

    with st.form("scenario_edits"):
        scenario_label = st.text_input(
            "Scenario label",
            value=str(config.get("scenarioLabel", "")),
        )
        screen_coverage = st.number_input(
            "Screen coverage",
            min_value=0.0,
            max_value=1.0,
            value=float(config.get("screenCoverage", 0.0)),
            step=0.05,
        )
        test_type = st.selectbox(
            "Test type",
            test_options,
            index=choice_index(config.get("testType"), test_options),
        )
        regimen = st.selectbox(
            "Regimen",
            regimen_options,
            index=choice_index(config.get("regimen"), regimen_options),
        )
        screening_strategy = st.selectbox(
            "Screening strategy",
            strategy_options,
            index=choice_index(config.get("screeningStrategy"), strategy_options, 3),
        )
        n_reps = st.number_input(
            "Simulation replicates",
            min_value=1,
            value=int(config.get("nReps", 1)),
            step=100,
        )
        submitted = st.form_submit_button("Apply edits")

    if submitted:
        updates = {
            "scenarioLabel": scenario_label,
            "screenCoverage": float(screen_coverage),
            "testType": test_type,
            "regimen": regimen,
            "screeningStrategy": screening_strategy,
            "nReps": int(n_reps),
        }
        changed = any(config.get(key) != value for key, value in updates.items())
        if changed:
            config.update(updates)
            st.session_state["config"] = config
            mark_config_changed()
            st.success("Scenario edits applied.")
        else:
            st.info("No scenario fields changed.")

    st.subheader("Current Config")
    st.dataframe(
        arrow_safe_dataframe(config_overview_rows(config)),
        width="content",
        hide_index=True,
    )

    with st.expander("Full config"):
        st.json(config, expanded=False)

    if st.button("Validate current config"):
        try:
            st.session_state["validation_report"] = backend.validate_config(config)
            sync_backend_status(backend.status())
            mark_validation_completed()
            st.success("Validation completed.")
        except Exception as exc:
            message = f"Validation failed: {exc}"
            sync_backend_status(backend.status())
            record_message("error", message)
            st.error(message)

    report = st.session_state.get("validation_report")
    if report:
        st.subheader("Validation Report")
        valid = report.get("isValid")
        if valid is True:
            st.success("Config is valid.")
        elif valid is False:
            st.error("Config has validation errors.")

        rows = validation_rows(report)
        if rows:
            st.dataframe(
                arrow_safe_dataframe(rows),
                width="stretch",
                hide_index=True,
            )
        else:
            st.json(report, expanded=False)

    st.subheader("Scenario Files")
    file_name = st.text_input("Scenario file", value="streamlit_scenario.json")
    cols = st.columns(2)
    if cols[0].button("Save scenario"):
        try:
            path = scenario_path(file_name)
            st.session_state["save_info"] = backend.save_scenario(
                config,
                str(path),
                st.session_state.get("economics_config"),
            )
            sync_backend_status(backend.status())
            st.success(f"Saved scenario to {path}")
        except Exception as exc:
            message = f"Save failed: {exc}"
            sync_backend_status(backend.status())
            record_message("error", message)
            st.error(message)

    if cols[1].button("Load scenario"):
        try:
            path = scenario_path(file_name)
            payload = json.loads(path.read_text(encoding="utf-8"))
            loaded_config = payload.get("config")
            if not isinstance(loaded_config, dict):
                raise ValueError("Scenario JSON does not contain a config object.")
            st.session_state["config"] = loaded_config
            st.session_state["economics_config"] = payload.get("economics")
            st.session_state["validation_report"] = None
            st.session_state["load_info"] = {
                "filename": str(path),
                "contractVersion": payload.get("contractVersion", ""),
                "scenarioLabel": payload.get("scenarioLabel", ""),
            }
            st.session_state["results_bundle"] = None
            st.session_state["results_stale"] = False
            st.session_state["dirty_config"] = True
            sync_backend_status(backend.status())
            st.success(f"Loaded scenario from {path}")
            st.rerun()
        except Exception as exc:
            message = f"Load failed: {exc}"
            sync_backend_status(backend.status())
            record_message("error", message)
            st.error(message)
else:
    st.info("Load MATLAB defaults to initialize the scenario state.")
