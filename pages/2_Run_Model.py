from __future__ import annotations

import streamlit as st

from app.display import arrow_safe_dataframe
from app.parameter_workspace import MODEL_METHOD_LABELS
from app.run_analysis_controls import prepare_run_config_for_recent_ltbi_route
from app.run_progress import StreamlitProgressDisplay, finalising_status, initialising_status
from app.state import (
    get_backend,
    init_session_state,
    is_supported_sa_health_analysis_basis,
    mark_run_completed,
    record_message,
    sanitize_reference_only_state,
    sync_backend_status,
)


def _page_link(path: str, *, label: str) -> None:
    try:
        st.page_link(path, label=label)
    except Exception:
        key = "nav_fallback_" + "".join(ch if ch.isalnum() else "_" for ch in f"{path}_{label}")
        st.button(label, disabled=True, key=key)


init_session_state()
st.session_state["apy_backend_name"] = "python_apy"
sanitize_reference_only_state()
backend = get_backend()

st.title("Run Analysis")

config = st.session_state.get("config")
if not config:
    st.info("Set up the analysis before running it.")
    st.stop()
notice = st.session_state.pop("reference_only_migration_notice", "")
if notice:
    st.warning(notice)

status = backend.status()
sync_backend_status(status)
if status.get("error"):
    st.error(status["error"])

method_label = MODEL_METHOD_LABELS.get(
    str(config.get("analysisMethod") or "expected_value"),
    "Quick deterministic preview - single expected-value calculation",
)
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
    run_type = "Quick deterministic preview"

st.subheader("Current run")
summary_rows = [
    {
        "Setting": "Analysis type",
        "Value": method_label,
    },
    {"Setting": "Run type", "Value": run_type},
    {
        "Setting": "Epidemiological basis",
        "Value": "Fixed SA Health report assumptions for future TB",
    },
    {"Setting": "Screening test", "Value": config.get("testType")},
    {"Setting": "Preventive treatment", "Value": config.get("regimen")},
    {"Setting": "Coverage", "Value": config.get("screenCoverage")},
]
if is_stochastic:
    summary_rows.insert(2, {"Setting": "Repetitions", "Value": f"{reps:,}"})
    summary_rows.insert(3, {"Setting": "Random seed", "Value": seed})
st.dataframe(arrow_safe_dataframe(summary_rows), use_container_width=True, hide_index=True)
st.caption(
    "Both options use the fixed SA Health report assumptions for future TB. "
    "The quick preview is a deterministic approximation. The 2,000-run analysis "
    "reproduces the report method and shows variation across simulated communities."
)
if is_stochastic and reps < 2000:
    st.warning("Preview analyses are useful for checking setup but do not reproduce the SA Health reference.")

if st.session_state.get("results_stale") or st.session_state.get("dirty_config"):
    st.warning("Inputs have changed. Run the analysis again before interpreting results.")
elif st.session_state.get("results_bundle"):
    st.success("Current results correspond to this configuration.")
else:
    st.info("No results have been generated for this setup.")

run_label = "Run analysis"

if st.button(run_label, type="primary"):
    try:
        progress = StreamlitProgressDisplay()
        progress.update(initialising_status())
        run_config = prepare_run_config_for_recent_ltbi_route(
            config,
            selected_route=None,
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
            _page_link("pages/0_Start.py", label="Return to Set up")
            st.stop()
        bundle = backend.run_scenario_bundle(
            run_config,
            validation_report=validation_report,
            progress_callback=progress.callback,
        )
        if not is_supported_sa_health_analysis_basis(bundle):
            bundle_metadata = bundle.get("metadata", {}) if isinstance(bundle, dict) else {}
            ledger_metadata = (
                (((bundle or {}).get("technical") or {}).get("eventLedger") or {}).get("metadata") or {}
                if isinstance(bundle, dict)
                else {}
            )
            st.error(
                "Analysis completed, but the result metadata did not satisfy the SA Health "
                "compatibility provenance contract. Results were not activated."
            )
            st.write(
                {
                    "bundleAnalysisBasis": bundle_metadata.get("analysisBasis"),
                    "bundleNaturalHistorySemantics": bundle_metadata.get("naturalHistorySemantics"),
                    "ledgerAnalysisBasis": ledger_metadata.get("analysisBasis"),
                    "ledgerNaturalHistorySemantics": ledger_metadata.get("naturalHistorySemantics"),
                }
            )
            st.stop()
        progress.update(finalising_status())
        st.session_state["results_bundle"] = bundle
        st.session_state["validation_report"] = validation_report
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
    _page_link("pages/3_Results.py", label="Open Results")
