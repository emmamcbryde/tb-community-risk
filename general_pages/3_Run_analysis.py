from __future__ import annotations

import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.state import (
    RESULTS_KEY,
    STALE_KEY,
    analysis_backend,
    get_analysis,
    get_intervention,
    get_profile,
    init_general_state,
    page_link,
    set_analysis,
    store_results,
)
from app.general.terminology import DETERMINISTIC_LABEL, STOCHASTIC_LABEL, display_text
from app.run_progress import StreamlitProgressDisplay, finalising_status, initialising_status
from engine.profiles.demonstration import GENERAL_DEFAULT_STOCHASTIC_SIMULATIONS
from engine.profiles.engine_mapping import ProfileMappingError, build_engine_config
from engine.profiles.population_profile import validate_profile


SECONDS_PER_SIMULATION_PER_10K = 0.9
METHODS = {"agent_based": STOCHASTIC_LABEL, "expected_value": DETERMINISTIC_LABEL}

init_general_state()
profile = get_profile()
intervention = get_intervention()
analysis = get_analysis()

st.title("Run analysis")

if st.session_state.get("general_analysis_type") not in METHODS:
    st.session_state["general_analysis_type"] = analysis["analysisMethod"]
method = st.radio(
    "Analysis type",
    list(METHODS),
    format_func=METHODS.get,
    key="general_analysis_type",
    horizontal=True,
)
updated = dict(analysis, analysisMethod=method)
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
    population = int(profile.population_size.value or 0)
    minutes = SECONDS_PER_SIMULATION_PER_10K * updated["nReps"] * population / 10_000 / 60
    st.caption(
        f"{updated['nReps']:,} simulated populations of {population:,} people. "
        f"Approximate run time: {max(minutes, 0.1):.0f} minute(s) on a typical computer."
    )
else:
    st.caption("A single expected-value calculation for a quick preview. It does not show variation between simulated populations.")
if updated != analysis:
    set_analysis(updated)
    st.rerun()

summary = [
    {"Setting": "Profile", "Value": profile.name},
    {"Setting": "Population size", "Value": f"{int(profile.population_size.value or 0):,}"},
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

blocking = [issue for issue in validate_profile(profile) if issue["severity"] in {"error", "blocking"}]
if blocking:
    st.error("The population profile needs correction before running:")
    for issue in blocking:
        st.write(f"- {issue['message']}")
    page_link("general_pages/1_Set_up_population.py", label="Return to Set up population")
    st.stop()

if st.session_state.get(RESULTS_KEY):
    if st.session_state.get(STALE_KEY):
        st.warning("Inputs have changed since the last run. Run the analysis again before interpreting results.")
    else:
        st.success("Current results match these inputs.")

if st.button("Run analysis", type="primary"):
    try:
        config = build_engine_config(profile, intervention=intervention, analysis=updated)
    except ProfileMappingError as exc:
        st.error(f"The profile cannot be analysed yet: {exc}")
        st.stop()
    backend = analysis_backend()
    try:
        progress = StreamlitProgressDisplay()
        progress.update(initialising_status())
        report = backend.validate_config(config)
        if not report.get("isValid"):
            st.error("The analysis settings are not valid. Review the inputs and try again.")
            st.stop()
        bundle = backend.run_scenario_bundle(config, validation_report=report, progress_callback=progress.callback)
        progress.update(finalising_status())
        store_results(bundle, config)
        st.success("Analysis completed.")
    except Exception as exc:  # show a plain message; details stay in the technical log
        st.session_state["general_last_error"] = repr(exc)
        st.error(
            display_text(
                f"The analysis could not be completed: {exc}",
                fallback="The analysis could not be completed. Details are recorded under Evidence and technical information.",
            )
        )

if st.session_state.get(RESULTS_KEY):
    page_link("general_pages/4_Results.py", label="Open Results")
