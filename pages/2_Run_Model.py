from __future__ import annotations

import streamlit as st

from app.display import arrow_safe_dataframe
from app.parameter_workspace import MODEL_METHOD_LABELS
from app.run_analysis_controls import (
    TECHNICAL_DEMONSTRATION_ROUTE,
    prepare_run_config_for_recent_ltbi_route,
)
from app.run_progress import StreamlitProgressDisplay, finalising_status, initialising_status
from app.state import (
    get_backend,
    init_session_state,
    mark_run_completed,
    record_message,
    sync_backend_status,
)
from engine.apy.ltbi_state import resolve_ltbi_state_assumptions


init_session_state()
st.session_state["apy_backend_name"] = "python_apy"
backend = get_backend()

st.title("Run Analysis")

config = st.session_state.get("config")
if not config:
    st.info("Set up the analysis before running it.")
    st.stop()

status = backend.status()
sync_backend_status(status)
if status.get("error"):
    st.error(status["error"])

method_label = MODEL_METHOD_LABELS.get(str(config.get("analysisMethod") or "expected_value"), "Expected outcomes")
is_stochastic = str(config.get("analysisMethod")) == "agent_based"
reps = int(float(config.get("nReps") or 0))
seed = int(float(config.get("seed") or 1))
if is_stochastic and reps == 2000 and seed == 1:
    run_type = "SA Health reference"
elif is_stochastic and reps < 2000:
    run_type = "Stochastic preview"
elif is_stochastic:
    run_type = "Modified stochastic run"
else:
    run_type = "Deterministic exploratory run"

st.subheader("Current run")
summary_rows = [
    {
        "Setting": "Analysis type",
        "Value": "Stochastic individual-based analysis" if is_stochastic else "Deterministic expected-value analysis",
    },
    {"Setting": "Run type", "Value": run_type},
    {"Setting": "Screening test", "Value": config.get("testType")},
    {"Setting": "Preventive treatment", "Value": config.get("regimen")},
    {"Setting": "Coverage", "Value": config.get("screenCoverage")},
]
if is_stochastic:
    summary_rows.insert(2, {"Setting": "Repetitions", "Value": f"{reps:,}"})
    summary_rows.insert(3, {"Setting": "Random seed", "Value": seed})
st.dataframe(arrow_safe_dataframe(summary_rows), use_container_width=True, hide_index=True)
if is_stochastic and reps < 2000:
    st.warning("Preview analyses are useful for checking setup but do not reproduce the SA Health reference.")

if st.session_state.get("results_stale") or st.session_state.get("dirty_config"):
    st.warning("Inputs have changed. Run the analysis again before interpreting results.")
elif st.session_state.get("results_bundle"):
    st.success("Current results correspond to this configuration.")
else:
    st.info("No results have been generated for this setup.")

ltbi_dev_compatibility_requested = False
ltbi_state = resolve_ltbi_state_assumptions(config)
unresolved_ltbi_state = ltbi_state.get("baselineRecentLTBIProportion") is None
if unresolved_ltbi_state:
    st.subheader("Recent versus remote LTBI assumption")
    st.warning("Choose the provisional working route on Set up, or review this assumption, before running.")
    decision_cols = st.columns(2)
    if decision_cols[0].button("Use provisional working route"):
        st.session_state["recent_ltbi_run_route"] = TECHNICAL_DEMONSTRATION_ROUTE
        st.rerun()
    decision_cols[1].page_link(
        "pages/6_Evidence_Assumptions.py",
        label="Review or enter the assumption",
    )
    if st.session_state.get("recent_ltbi_run_route") == TECHNICAL_DEMONSTRATION_ROUTE:
        ltbi_dev_compatibility_requested = True
        st.caption("Provisional route selected. Detailed caveats are in Evidence & Assumptions.")

run_label = "Run analysis"

if st.button(run_label, type="primary"):
    try:
        progress = StreamlitProgressDisplay()
        progress.update(initialising_status())
        ltbi_state = resolve_ltbi_state_assumptions(config)
        if ltbi_state.get("baselineRecentLTBIProportion") is None:
            if not ltbi_dev_compatibility_requested:
                st.info(
                    "Choose the provisional working route, or review and enter the "
                    "recent-versus-remote LTBI assumption, before running the analysis."
                )
                st.stop()
        run_config = prepare_run_config_for_recent_ltbi_route(
            config,
            selected_route=st.session_state.get("recent_ltbi_run_route"),
        )
        if run_config != config:
            st.session_state["config"] = run_config
        validation_report = backend.validate_config(run_config)
        if not validation_report.get("isValid"):
            st.session_state["validation_report"] = validation_report
            blocking = []
            for issue in validation_report.get("errors") or []:
                if isinstance(issue, dict):
                    blocking.append(issue.get("message") or str(issue))
                else:
                    blocking.append(str(issue))
            st.error("This setup is not valid. Return to Set up before running.")
            if blocking:
                st.write(blocking[0])
            st.page_link("pages/0_Start.py", label="Return to Set up")
            st.stop()
        bundle = backend.run_scenario_bundle(
            run_config,
            validation_report=validation_report,
            progress_callback=progress.callback,
        )
        progress.update(finalising_status())
        st.session_state["results_bundle"] = bundle
        st.session_state["validation_report"] = validation_report
        st.session_state["economics_results"] = None
        st.session_state["dirty_economics"] = True
        st.session_state["economics_config"] = None
        sync_backend_status(backend.status())
        mark_run_completed()
        st.success("Analysis completed.")
        st.page_link("pages/3_Results.py", label="Open Results")
    except Exception as exc:
        message = f"Analysis failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

if st.session_state.get("results_bundle"):
    st.page_link("pages/3_Results.py", label="Open Results")
