from __future__ import annotations

import streamlit as st

from app.state import (
    get_backend,
    init_session_state,
    mark_run_completed,
    mark_validation_completed,
    record_message,
    sync_backend_status,
)


init_session_state()
backend = get_backend()

st.title("Run Model")
st.caption("Runs the frozen MATLAB v9 backend and stores the portable results bundle.")

config = st.session_state.get("config")
if not config:
    st.info("Load a scenario first.")
    st.stop()

status = backend.status()
sync_backend_status(status)
st.caption(
    f"Backend: {status.get('name', 'unknown')} | MATLAB started: "
    f"{'yes' if status.get('started') else 'no'}"
)
if status.get("error"):
    st.error(status["error"])

if st.session_state.get("dirty_config"):
    st.warning("Current config has changed since the last validation or run.")
if st.session_state.get("results_stale"):
    st.warning("Stored results are stale because inputs changed after the last run.")
elif st.session_state.get("results_bundle"):
    st.success("Stored results match the current config state.")

if st.button("Validate scenario"):
    try:
        st.session_state["validation_report"] = backend.validate_config(config)
        mark_validation_completed()
        sync_backend_status(backend.status())
        st.success("Validation completed.")
    except Exception as exc:
        message = f"Validation failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

report = st.session_state.get("validation_report")
if report:
    st.subheader("Validation Report")
    if report.get("isValid") is True:
        st.success("Config is valid.")
    elif report.get("isValid") is False:
        st.error("Config has validation errors.")
    st.json(report, expanded=False)

run_label = "Run MATLAB model"
if st.session_state.get("dirty_config"):
    run_label = "Validate and run MATLAB model"

if st.button(run_label, type="primary"):
    try:
        bundle = backend.run_scenario_bundle(
            config,
            validation_report=st.session_state.get("validation_report"),
        )
        st.session_state["results_bundle"] = bundle
        validation = bundle.get("validation", {})
        if isinstance(validation, dict) and isinstance(validation.get("report"), dict):
            st.session_state["validation_report"] = validation["report"]
        st.session_state["economics_results"] = None
        st.session_state["dirty_economics"] = True
        sync_backend_status(backend.status())
        mark_run_completed()
        st.success("Model run completed.")
    except Exception as exc:
        message = f"Model run failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

if st.session_state.get("results_bundle"):
    metadata = st.session_state["results_bundle"].get("metadata", {})
    st.subheader("Latest Run")
    st.json(
        {
            "last_run_at": st.session_state.get("last_run_at"),
            "results_stale": st.session_state.get("results_stale"),
            "scenarioLabel": metadata.get("scenarioLabel"),
            "modelVersion": metadata.get("modelVersion"),
            "contractVersion": metadata.get("contractVersion"),
        },
        expanded=False,
    )
