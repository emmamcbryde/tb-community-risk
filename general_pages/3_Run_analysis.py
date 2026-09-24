from __future__ import annotations

import time

import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.state import (
    analysis_backend,
    cached_results_for,
    current_engine_config,
    get_analysis,
    get_intervention,
    get_profile,
    init_general_state,
    page_link,
    results_status,
    set_analysis,
    store_results,
)
from app.general.terminology import DETERMINISTIC_LABEL, STOCHASTIC_LABEL, display_text
from app.run_progress import StreamlitProgressDisplay, finalising_status, initialising_status
from engine.profiles.demonstration import GENERAL_DEFAULT_STOCHASTIC_SIMULATIONS
from engine.profiles.population_profile import Provenance


# Measured on the reference machine (see docs/performance_benchmark.md).
SECONDS_PER_SIMULATION_PER_10K = 0.55
DETERMINISTIC_SECONDS = 30
CALIBRATION_SECONDS = 40
ECONOMICS_SECONDS_PER_SIMULATION = 0.18
CONFIRM_ABOVE_MINUTES = 3
METHODS = {"agent_based": STOCHASTIC_LABEL, "expected_value": DETERMINISTIC_LABEL}

init_general_state()
profile = get_profile()
intervention = get_intervention()
analysis = get_analysis()

st.title("Run analysis")

if st.session_state.get("general_analysis_type") not in METHODS:
    st.session_state["general_analysis_type"] = analysis["analysisMethod"]
method = st.radio("Analysis type", list(METHODS), format_func=METHODS.get, key="general_analysis_type", horizontal=True)
updated = dict(analysis, analysisMethod=method)
population = int(profile.population_size.value or 0)
if method == "agent_based":
    if st.session_state.get("general_n_sims") is None:
        st.session_state["general_n_sims"] = int(analysis["nReps"])
    updated["nReps"] = int(
        st.number_input(
            "Number of simulated populations",
            min_value=1,
            max_value=5000,
            step=100,
            key="general_n_sims",
            help=f"Default {GENERAL_DEFAULT_STOCHASTIC_SIMULATIONS:,}. Fewer simulations run faster but describe variation less precisely.",
        )
    )
    if st.session_state.get("general_seed") is None:
        st.session_state["general_seed"] = int(analysis["seed"])
    updated["seed"] = int(st.number_input("Random seed", min_value=0, max_value=2_147_483_647, step=1, key="general_seed"))
    minutes = (SECONDS_PER_SIMULATION_PER_10K * updated["nReps"] * population / 10_000 + CALIBRATION_SECONDS) / 60
    st.info(
        f"Stochastic analysis: {updated['nReps']:,} simulated populations of {population:,} people. "
        f"Estimated run time about {max(minutes, 1):.0f} minutes on the reference computer, plus about "
        f"{max(ECONOMICS_SECONDS_PER_SIMULATION * updated['nReps'] / 60, 0.1):.0f} minutes for health economics; "
        "keep this page open while it runs."
    )
else:
    minutes = (DETERMINISTIC_SECONDS + CALIBRATION_SECONDS) / 60
    st.info(
        "Deterministic expected-value preview: a single calculation, usually about one minute. "
        "It does not show variation between simulated populations."
    )
if updated != analysis:
    set_analysis(updated)
    st.rerun()

summary = [
    {"Setting": "Profile", "Value": profile.name},
    {"Setting": "Population size", "Value": f"{population:,}"},
    {"Setting": "Analysis type", "Value": METHODS[method]},
    {"Setting": "Screening test", "Value": intervention["testType"]},
    {"Setting": "Preventive treatment", "Value": intervention["regimen"]},
    {"Setting": "Screening coverage", "Value": f"{intervention['screenCoverage'] * 100:.0f}%"},
]
if method == "agent_based":
    summary.insert(3, {"Setting": "Simulated populations", "Value": f"{updated['nReps']:,}"})
st.dataframe(arrow_safe_dataframe(summary), use_container_width=True, hide_index=True)
st.caption(
    "Results describe direct effects for people screened and treated; transmission benefits are not included. "
    "Future TB risk uses the demonstration epidemiological assumptions."
)
if profile.incidence.provenance in {Provenance.WHO_SNAPSHOT, Provenance.LOCAL_UPLOAD}:
    st.caption(
        "The applied incidence data describe the TB burden and trend only; they do not yet change the model's "
        "epidemiology, so results are not a country-specific estimate."
    )

config, error = current_engine_config()
if config is None:
    st.error(f"The profile cannot be analysed yet: {error}")
    page_link("general_pages/1_Set_up_population.py", label="Return to Set up population")
    st.stop()

run_message = st.session_state.pop("general_run_message", "")
if run_message:
    st.success(run_message)
status = results_status()
if status == "current":
    st.success("Current results already match these inputs; running again would reproduce them.")
elif status == "stale":
    st.warning("Inputs that affect health outcomes have changed since the last run. Run the analysis again.")

confirmed = True
if minutes > CONFIRM_ABOVE_MINUTES:
    confirmed = st.checkbox(
        f"I understand this analysis takes about {minutes:.0f} minutes.",
        key="general_confirm_long_run",
    )
run_clicked = st.button("Run analysis", type="primary", disabled=not confirmed or status == "current")

if run_clicked and results_status() == "current":
    st.info("Results are already current for these inputs; no new run was started.")
elif run_clicked:
    cached = cached_results_for(config)
    if cached is not None:
        store_results(cached["bundle"], cached["config"])
        st.session_state["general_run_message"] = "Reused the completed analysis for these identical inputs."
        st.rerun()
    else:
        backend = analysis_backend()
        try:
            progress = StreamlitProgressDisplay()
            progress.update(initialising_status())
            report = backend.validate_config(config)
            if not report.get("isValid"):
                st.error("The analysis settings are not valid. Review the inputs and try again.")
                st.stop()
            started = time.perf_counter()
            with st.spinner("Running analysis..."):
                bundle = backend.run_scenario_bundle(config, validation_report=report, progress_callback=progress.callback)
            progress.update(finalising_status())
            store_results(bundle, config)
            st.session_state["general_last_run_seconds"] = time.perf_counter() - started
            st.session_state["general_run_message"] = f"Analysis completed in {st.session_state['general_last_run_seconds']:.0f} seconds."
        except Exception as exc:  # show a plain message; details stay in the technical log
            st.session_state["general_last_error"] = repr(exc)
            st.error(
                display_text(
                    f"The analysis could not be completed: {exc}",
                    fallback="The analysis could not be completed. Details are recorded under Evidence and technical information.",
                )
            )

        if st.session_state.get("general_run_message"):
            st.rerun()

if st.session_state.get("general_results_bundle"):
    page_link("general_pages/4_Results.py", label="Open Results")
