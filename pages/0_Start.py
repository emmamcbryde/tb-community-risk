from __future__ import annotations

import streamlit as st

from app.state import get_backend, init_session_state, record_message, sync_backend_status
from engine.apy.scenario import (
    DEFAULT_POPULATION_PRESET_ID,
    build_scenario_contract,
    config_updates_from_scenario,
)


init_session_state()

st.title("LTBI Screening Decision Tool")
st.write(
    "Compare the direct benefits, harms, resource requirements and health-system "
    "consequences of targeted LTBI screening and preventive treatment."
)
st.info(
    "This tool supports planning and sequencing decisions. It does not recommend "
    "denying care to any person or group."
)

col_demo, col_new, col_results = st.columns(3)

with col_demo:
    if st.button("Start with APY demonstration", type="primary", use_container_width=True):
        try:
            st.session_state["apy_backend_name"] = "python_apy"
            backend = get_backend()
            config = backend.default_config()
            scenario = build_scenario_contract(DEFAULT_POPULATION_PRESET_ID)
            config.update(config_updates_from_scenario(scenario))
            st.session_state["config"] = config
            st.session_state["validation_report"] = None
            st.session_state["results_bundle"] = None
            st.session_state["economics_results"] = None
            st.session_state["dirty_config"] = False
            st.session_state["dirty_economics"] = False
            st.session_state["results_stale"] = False
            st.session_state.pop("recent_ltbi_run_route", None)
            sync_backend_status(backend.status())
            st.success("APY demonstration loaded. Continue to Define Strategy.")
            st.page_link("pages/1_Scenario.py", label="Open Define Strategy")
        except Exception as exc:
            record_message("error", f"Could not load APY demonstration: {exc}")
            st.error("Could not load the APY demonstration.")

with col_new:
    if st.button("Create a new analysis", use_container_width=True):
        try:
            st.session_state["apy_backend_name"] = "python_apy"
            backend = get_backend()
            st.session_state["config"] = backend.default_config()
            st.session_state["validation_report"] = None
            st.session_state["results_bundle"] = None
            st.session_state["economics_results"] = None
            st.session_state["dirty_config"] = False
            st.session_state["dirty_economics"] = False
            st.session_state["results_stale"] = False
            st.session_state.pop("recent_ltbi_run_route", None)
            sync_backend_status(backend.status())
            st.success("New analysis created. Continue to Define Strategy.")
            st.page_link("pages/1_Scenario.py", label="Open Define Strategy")
        except Exception as exc:
            record_message("error", f"Could not create analysis: {exc}")
            st.error("Could not create a new analysis.")

with col_results:
    if st.session_state.get("results_bundle"):
        st.page_link("pages/3_Results.py", label="Continue to current results")
    else:
        st.button("Continue to current results", disabled=True, use_container_width=True)

with st.expander("Research and Development", expanded=False):
    st.write(
        "Technical implementation details, model version identifiers and diagnostics "
        "are available in the Research and Development section."
    )
    st.page_link("pages/8_Technical_Settings.py", label="Open Technical settings")
