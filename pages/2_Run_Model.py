from __future__ import annotations

import streamlit as st

from app.run_analysis_controls import (
    TECHNICAL_DEMONSTRATION_ROUTE,
    prepare_run_config_for_recent_ltbi_route,
)
from app.run_progress import StreamlitProgressDisplay, finalising_status, initialising_status
from app.state import (
    get_backend,
    init_session_state,
    mark_run_completed,
    mark_validation_completed,
    record_message,
    sync_backend_status,
)
from engine.apy.ltbi_state import resolve_ltbi_state_assumptions


init_session_state()
st.session_state["apy_backend_name"] = "python_apy"
backend = get_backend()

st.title("Run Analysis")
st.caption("Validate inputs and run the LTBI screening analysis.")

config = st.session_state.get("config")
if not config:
    st.info("Define a strategy first.")
    st.stop()

status = backend.status()
sync_backend_status(status)
if status.get("error"):
    st.error(status["error"])

if st.session_state.get("dirty_config"):
    st.warning("Current inputs have changed since the last validation or run.")
if st.session_state.get("results_stale"):
    st.warning("Stored results are stale because inputs changed after the last run.")
elif st.session_state.get("results_bundle"):
    st.success("Stored results match the current inputs.")

ltbi_dev_compatibility_requested = False
ltbi_state = resolve_ltbi_state_assumptions(config)
unresolved_ltbi_state = ltbi_state.get("baselineRecentLTBIProportion") is None
if ltbi_state.get("warning"):
    st.warning(str(ltbi_state["warning"]))
if unresolved_ltbi_state:
    st.subheader("Recent versus remote LTBI assumption")
    st.write(
        "The proportion of baseline infections that were acquired relatively "
        "recently has not yet been established for this demonstration population. "
        "This affects the estimated risk of progression to active TB."
    )
    decision_cols = st.columns(2)
    if decision_cols[0].button("Run a technical demonstration"):
        st.session_state["recent_ltbi_run_route"] = TECHNICAL_DEMONSTRATION_ROUTE
        st.info(
            "Technical demonstration selected. The analysis will use the existing "
            "0% compatibility placeholder, temporarily representing all baseline "
            "infection as remote. Outputs will be provisional and not reference, "
            "reviewed or clinician-ready results."
        )
    decision_cols[1].page_link(
        "pages/4_Economics.py",
        label="Review or enter the assumption",
    )
    if st.session_state.get("recent_ltbi_run_route") == TECHNICAL_DEMONSTRATION_ROUTE:
        ltbi_dev_compatibility_requested = True
        st.warning(
            "Technical demonstration mode is selected. Every output from this run "
            "will remain provisional, and evidence-review status will not be promoted."
        )

if st.button("Validate inputs"):
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
    st.subheader("Input Checks")
    if report.get("isValid") is True:
        st.success("Inputs are valid.")
    elif report.get("isValid") is False:
        st.error("Inputs have validation errors.")
    st.json(report, expanded=False)

run_label = "Run analysis"
if st.session_state.get("dirty_config"):
    run_label = "Validate and run analysis"

if st.button(run_label, type="primary"):
    try:
        progress = StreamlitProgressDisplay()
        progress.update(initialising_status())
        ltbi_state = resolve_ltbi_state_assumptions(config)
        if ltbi_state.get("baselineRecentLTBIProportion") is None:
            if not ltbi_dev_compatibility_requested:
                st.info(
                    "Choose Run a technical demonstration, or review and enter the "
                    "recent-versus-remote LTBI assumption, before running the analysis."
                )
                st.stop()
        run_config = prepare_run_config_for_recent_ltbi_route(
            config,
            selected_route=st.session_state.get("recent_ltbi_run_route"),
        )
        if run_config != config:
            st.session_state["config"] = run_config
        bundle = backend.run_scenario_bundle(
            run_config,
            validation_report=st.session_state.get("validation_report"),
            progress_callback=progress.callback,
        )
        progress.update(finalising_status())
        st.session_state["results_bundle"] = bundle
        validation = bundle.get("validation", {})
        if isinstance(validation, dict) and isinstance(validation.get("report"), dict):
            st.session_state["validation_report"] = validation["report"]
        st.session_state["economics_results"] = None
        st.session_state["dirty_economics"] = True
        st.session_state["economics_config"] = None
        sync_backend_status(backend.status())
        mark_run_completed()
        st.success("Analysis completed.")
    except Exception as exc:
        message = f"Analysis failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

if st.session_state.get("results_bundle"):
    metadata = st.session_state["results_bundle"].get("metadata", {})
    st.subheader("Latest Analysis")
    st.json(
        {
            "last_run_at": st.session_state.get("last_run_at"),
            "results_stale": st.session_state.get("results_stale"),
            "analysisLabel": metadata.get("scenarioLabel"),
        },
        expanded=False,
    )
    with st.expander("Technical information", expanded=False):
        st.json(
            {
                "modelVersion": metadata.get("modelVersion"),
                "contractVersion": metadata.get("contractVersion"),
            },
            expanded=False,
        )
