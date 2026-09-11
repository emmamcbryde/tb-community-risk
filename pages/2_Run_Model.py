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
    mark_validation_completed,
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
run_type = "SA Health reference" if is_stochastic and reps == 2000 and seed == 1 else ("Exploratory preview" if is_stochastic and reps < 2000 else method_label)

st.subheader("Current run")
summary_rows = [
    {"Setting": "Analysis type", "Value": method_label},
    {"Setting": "Repetitions", "Value": f"{reps:,}" if is_stochastic else "Not used"},
    {"Setting": "Random seed", "Value": seed if is_stochastic else "Not used"},
    {"Setting": "Run type", "Value": run_type},
    {"Setting": "Screening test", "Value": config.get("testType")},
    {"Setting": "Preventive treatment", "Value": config.get("regimen")},
    {"Setting": "Coverage", "Value": config.get("screenCoverage")},
]
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
    if decision_cols[0].button("Run provisional working analysis"):
        st.session_state["recent_ltbi_run_route"] = TECHNICAL_DEMONSTRATION_ROUTE
        st.info(
            "Provisional route selected. The analysis will use the existing "
            "0% compatibility placeholder, temporarily representing all baseline "
            "infection as remote. Outputs will be provisional and not reference, "
            "reviewed or clinician-ready results."
        )
    decision_cols[1].page_link(
        "pages/6_Evidence_Assumptions.py",
        label="Review or enter the assumption",
    )
    if st.session_state.get("recent_ltbi_run_route") == TECHNICAL_DEMONSTRATION_ROUTE:
        ltbi_dev_compatibility_requested = True
        st.warning(
            "Provisional working-analysis route is selected. Every output from this run "
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
    if report.get("isValid") is True:
        st.success("Inputs are valid.")
    elif report.get("isValid") is False:
        st.error("Inputs have validation errors.")
    issue_rows = []
    for group in ("errors", "warnings"):
        for issue in report.get(group) or []:
            if isinstance(issue, dict):
                issue_rows.append(
                    {
                        "Severity": group[:-1],
                        "Field": issue.get("fieldLabel") or issue.get("field") or "",
                        "Message": issue.get("message") or "",
                    }
                )
    if issue_rows:
        st.dataframe(arrow_safe_dataframe(issue_rows), use_container_width=True, hide_index=True)

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
                    "Choose Run provisional working analysis, or review and enter the "
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
        st.page_link("pages/3_Results.py", label="Open Results")
    except Exception as exc:
        message = f"Analysis failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

if st.session_state.get("results_bundle"):
    st.page_link("pages/3_Results.py", label="Open Results")
