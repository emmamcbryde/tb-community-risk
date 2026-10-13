from __future__ import annotations

from pathlib import Path
import json
import sys
from typing import Any

import pandas as pd
import streamlit as st

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.import_bootstrap import ensure_repo_root

ensure_repo_root(REPO_ROOT)

from app.explicit_tbi_workflow import (
    DETERMINISTIC_LABEL,
    EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
    STOCHASTIC_LABEL,
    TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION,
    WORKFLOW_PAGES,
    active_tb_validation_for_state,
    build_explicit_workflow_config,
    default_workflow_state,
    discard_incompatible_result,
    health_signature,
    percent_to_proportion,
    proportion_to_percent,
    provenance_text,
    results_to_csv,
    run_explicit_workflow_analysis,
    source_for_value,
    validate_workflow_state,
    validation_to_csv,
)


STATE_KEY = "explicit_tbi_workflow_state"
RESULT_KEY = "explicit_tbi_results"
STALE_KEY = "explicit_tbi_results_stale"
PAGE_KEY = "explicit_tbi_page"
WIDGET_VERSION_KEY = "explicit_tbi_widget_version"


st.set_page_config(
    page_title="Explicit recent/remote TBI development workflow",
    layout="wide",
)


def _init_state() -> None:
    st.session_state[STATE_KEY] = dict(
        st.session_state.get(STATE_KEY) or default_workflow_state()
    )
    st.session_state[RESULT_KEY] = discard_incompatible_result(
        st.session_state.get(RESULT_KEY)
    )
    st.session_state.setdefault(STALE_KEY, False)
    st.session_state.setdefault(PAGE_KEY, "Set up")
    st.session_state.setdefault(WIDGET_VERSION_KEY, 0)
    result = st.session_state.get(RESULT_KEY)
    if isinstance(result, dict):
        st.session_state[STALE_KEY] = (
            result.get("healthSignature") != health_signature(st.session_state[STATE_KEY])
        )


def _restore_defaults() -> None:
    st.session_state[STATE_KEY] = default_workflow_state()
    st.session_state[RESULT_KEY] = None
    st.session_state[STALE_KEY] = False
    st.session_state[PAGE_KEY] = "Set up"
    st.session_state[WIDGET_VERSION_KEY] = int(
        st.session_state.get(WIDGET_VERSION_KEY, 0)
    ) + 1


def _set_page(page: str) -> None:
    st.session_state[PAGE_KEY] = page


def _state() -> dict[str, Any]:
    return st.session_state[STATE_KEY]


def _widget_key(name: str) -> str:
    return f"explicit_tbi_{st.session_state[WIDGET_VERSION_KEY]}_{name}"


def _percent_input(label: str, key: str, help_text: str) -> None:
    state = _state()
    value = st.number_input(
        label,
        min_value=0.0,
        max_value=100.0,
        value=proportion_to_percent(float(state[key])),
        step=0.1,
        format="%.6f",
        help=help_text,
        key=_widget_key(key),
    )
    state[key] = percent_to_proportion(value)
    st.caption(
        f"Value used by model: {state[key]:.8f}; unit: proportion of total modelled population; "
        f"source: {source_for_value(state, key)}"
    )


def _setup_page() -> None:
    state = _state()
    old_signature = health_signature(state)
    st.subheader("Population")
    cols = st.columns(3)
    with cols[0]:
        state["population"] = int(
            st.number_input(
                "Total modelled population",
                min_value=1,
                max_value=1_000_000,
                value=int(state["population"]),
                step=100,
                key=_widget_key("population"),
            )
        )
    with cols[1]:
        state["analysisHorizonYears"] = float(
            st.number_input(
                "Analysis horizon",
                min_value=0.5,
                max_value=100.0,
                value=float(state["analysisHorizonYears"]),
                step=1.0,
                key=_widget_key("analysis_horizon"),
            )
        )
        st.caption("Unit: years")
    with cols[2]:
        baseline_pct = st.number_input(
            "Baseline active TB",
            min_value=0.0,
            max_value=100.0,
            value=proportion_to_percent(float(state["baselineActiveTBProportion"])),
            step=0.01,
            format="%.6f",
            help="Prevalent active TB at baseline is sequenced before TBI states.",
            key=_widget_key("baseline_active_tb"),
        )
        state["baselineActiveTBProportion"] = percent_to_proportion(baseline_pct)
        st.caption(
            f"Value used by model: {state['baselineActiveTBProportion']:.8f}; "
            f"source: {source_for_value(state, 'baselineActiveTBProportion')}"
        )

    st.subheader("Baseline TBI history")
    tbi_cols = st.columns(2)
    with tbi_cols[0]:
        _percent_input(
            "Recently infected within the past 5 years",
            "recentTBITarget",
            "At least one TB infection acquired in the five years immediately before baseline.",
        )
    with tbi_cols[1]:
        _percent_input(
            "Remote infection only — infected 5 or more years ago and not reinfected within the past 5 years",
            "remoteOnlyTBITarget",
            "At least one infection more than five years before baseline and no infection in the most recent five years.",
        )
    implied_uninfected = max(
        0.0,
        1.0
        - float(state["baselineActiveTBProportion"])
        - float(state["recentTBITarget"])
        - float(state["remoteOnlyTBITarget"]),
    )
    st.metric("Implied uninfected prevalence", f"{100 * implied_uninfected:.3f}%")

    with st.expander("Effective age distribution", expanded=False):
        try:
            cfg = build_explicit_workflow_config(state)
            calibration = cfg and __import__(
                "engine.apy.calibration_policy",
                fromlist=["resolve_calibration_for_config"],
            ).resolve_calibration_for_config(cfg)
            params = calibration["parameters"]
            table = pd.DataFrame(
                {
                    "Age": params["exactAgeValues"],
                    "Population share": params["exactAgeProb"],
                }
            )
            st.dataframe(table, hide_index=True, use_container_width=True)
            st.caption("The open-ended 85+ source category uses the explicit age85PlusMax provenance from the runner contract.")
        except Exception as exc:
            st.warning(f"Age distribution cannot be displayed yet: {exc}")

    st.subheader("Intervention")
    int_cols = st.columns(3)
    with int_cols[0]:
        state["screeningCoverage"] = percent_to_proportion(
            st.number_input(
                "Screening coverage",
                min_value=0.0,
                max_value=100.0,
                value=proportion_to_percent(float(state["screeningCoverage"])),
                step=5.0,
                key=_widget_key("screening_coverage"),
            )
        )
        state["screeningWindowYears"] = float(
            st.number_input(
                "Screening period",
                min_value=0.01,
                max_value=20.0,
                value=float(state["screeningWindowYears"]),
                step=0.5,
                key=_widget_key("screening_period"),
            )
        )
        st.caption(
            "Coverage means the proportion of eligible people screened at least "
            "once over the complete screening programme window."
        )
    with int_cols[1]:
        state["testType"] = st.selectbox(
            "Test",
            ["IGRA", "TST"],
            index=["IGRA", "TST"].index(str(state["testType"])),
            key=_widget_key("test_type"),
        )
        state["testSensitivity"] = percent_to_proportion(
            st.number_input(
                "Test sensitivity",
                min_value=0.0,
                max_value=100.0,
                value=proportion_to_percent(float(state["testSensitivity"])),
                step=1.0,
                key=_widget_key("test_sensitivity"),
            )
        )
        state["testSpecificity"] = percent_to_proportion(
            st.number_input(
                "Test specificity",
                min_value=0.0,
                max_value=100.0,
                value=proportion_to_percent(float(state["testSpecificity"])),
                step=1.0,
                key=_widget_key("test_specificity"),
            )
        )
    with int_cols[2]:
        regimens = ["3HP", "4R", "3HR", "6H", "9H"]
        state["regimen"] = st.selectbox(
            "Preventive-treatment regimen",
            regimens,
            index=regimens.index(str(state["regimen"])),
            key=_widget_key("regimen"),
        )
        state["pStartTPT"] = percent_to_proportion(
            st.number_input(
                "Treatment initiation",
                min_value=0.0,
                max_value=100.0,
                value=proportion_to_percent(float(state["pStartTPT"])),
                step=5.0,
                key=_widget_key("p_start_tpt"),
            )
        )
        state["regimenPComplete"] = percent_to_proportion(
            st.number_input(
                "Treatment completion",
                min_value=0.0,
                max_value=100.0,
                value=proportion_to_percent(float(state["regimenPComplete"])),
                step=5.0,
                key=_widget_key("regimen_complete"),
            )
        )
        state["regimenEffFull"] = percent_to_proportion(
            st.number_input(
                "TPT efficacy after completion",
                min_value=0.0,
                max_value=100.0,
                value=proportion_to_percent(float(state["regimenEffFull"])),
                step=5.0,
                key=_widget_key("regimen_eff_full"),
            )
        )

    st.subheader("Analysis controls")
    state["analysisMethod"] = st.radio(
        "Analysis choice",
        [DETERMINISTIC_LABEL, STOCHASTIC_LABEL],
        index=[DETERMINISTIC_LABEL, STOCHASTIC_LABEL].index(state["analysisMethod"]),
        key=_widget_key("analysis_method"),
    )
    if state["analysisMethod"] == STOCHASTIC_LABEL:
        acols = st.columns(2)
        with acols[0]:
            state["stochasticRepetitions"] = int(
                st.number_input(
                    "Repetitions",
                    min_value=1,
                    max_value=5000,
                    value=int(state["stochasticRepetitions"]),
                    step=1,
                    key=_widget_key("stochastic_repetitions"),
                )
            )
        with acols[1]:
            state["seed"] = int(
                st.number_input(
                    "Random seed",
                    min_value=0,
                    max_value=2_147_483_647,
                    value=int(state["seed"]),
                    step=1,
                    key=_widget_key("seed"),
                )
            )
        st.caption("Stochastic runs simulate multiple populations; larger repetition counts take longer.")

    with st.expander("Optional active-TB observations", expanded=False):
        state["activeTBObservationEnabled"] = st.checkbox(
            "Include active-TB observation for validation",
            value=bool(state["activeTBObservationEnabled"]),
            key=_widget_key("active_tb_observation_enabled"),
        )
        st.info(
            "These observations are used to compare observed and model-expected active TB. "
            "They do not automatically recalibrate the progression curve."
        )
        if state["activeTBObservationEnabled"]:
            obs_cols = st.columns(3)
            with obs_cols[0]:
                state["activeTBObservedCases"] = float(
                    st.number_input(
                        "Observed active-TB cases",
                        min_value=0.0,
                        value=float(state["activeTBObservedCases"]),
                        step=1.0,
                        key=_widget_key("active_tb_observed_cases"),
                    )
                )
                state["activeTBDenominator"] = float(
                    st.number_input(
                        "Observation denominator",
                        min_value=1.0,
                        value=float(state["activeTBDenominator"]),
                        step=100.0,
                        key=_widget_key("active_tb_denominator"),
                    )
                )
            with obs_cols[1]:
                state["activeTBObservationDurationYears"] = float(
                    st.number_input(
                        "Observation duration",
                        min_value=0.01,
                        value=float(state["activeTBObservationDurationYears"]),
                        step=0.5,
                        key=_widget_key("active_tb_duration"),
                    )
                )
                classifications = [
                    "baseline/prevalent active TB",
                    "screen-detected active TB",
                    "prospective incident active TB",
                    "retrospective notifications/incidence",
                    "mixed or insufficient evidence",
                ]
                state["activeTBObservationClassification"] = st.selectbox(
                    "Observation classification",
                    classifications,
                    index=classifications.index(
                        state["activeTBObservationClassification"]
                    ),
                    key=_widget_key("active_tb_classification"),
                )
            with obs_cols[2]:
                state["activeTBAscertainment"] = st.text_input(
                    "Ascertainment description",
                    value=str(state["activeTBAscertainment"]),
                    key=_widget_key("active_tb_ascertainment"),
                )
                state["activeTBAscertainmentProbability"] = percent_to_proportion(
                    st.number_input(
                        "Ascertainment probability",
                        min_value=0.000001,
                        max_value=100.0,
                        value=proportion_to_percent(
                            float(state["activeTBAscertainmentProbability"])
                        ),
                        step=1.0,
                        key=_widget_key("active_tb_ascertainment_probability"),
                    )
                )
                state["activeTBObservationSource"] = st.text_input(
                    "Observation source",
                    value=str(state["activeTBObservationSource"]),
                    key=_widget_key("active_tb_source"),
                )

    errors = validate_workflow_state(state)
    if errors:
        st.error(errors[0])
    else:
        st.success("Inputs are valid for the isolated explicit pathway.")
    if st.session_state.get(RESULT_KEY) and health_signature(state) != old_signature:
        st.session_state[STALE_KEY] = True

    bottom = st.columns(2)
    with bottom[0]:
        st.button("Restore working defaults", on_click=_restore_defaults, use_container_width=True)
    with bottom[1]:
        st.button(
            "Continue to Run analysis",
            on_click=_set_page,
            args=("Run analysis",),
            type="primary",
            use_container_width=True,
        )


def _run_page() -> None:
    state = _state()
    errors = validate_workflow_state(state)
    st.subheader("Pre-run summary")
    implied_uninfected = (
        1.0
        - float(state["baselineActiveTBProportion"])
        - float(state["recentTBITarget"])
        - float(state["remoteOnlyTBITarget"])
    )
    rows = [
        ("Population", f"{int(state['population']):,}"),
        ("Recent TBI target", f"{100 * state['recentTBITarget']:.3f}%"),
        ("Remote-only TBI target", f"{100 * state['remoteOnlyTBITarget']:.3f}%"),
        ("Implied uninfected", f"{100 * implied_uninfected:.3f}%"),
        ("Baseline active TB", f"{100 * state['baselineActiveTBProportion']:.3f}%"),
        ("Test and regimen", f"{state['testType']} / {state['regimen']}"),
        ("Coverage", f"{100 * state['screeningCoverage']:.1f}%"),
        ("Horizon", f"{state['analysisHorizonYears']:g} years"),
        ("Analysis method", state["analysisMethod"]),
    ]
    if state["analysisMethod"] == STOCHASTIC_LABEL:
        rows.extend(
            [
                ("Repetitions", f"{int(state['stochasticRepetitions']):,}"),
                ("Random seed", str(int(state["seed"]))),
            ]
        )
    st.dataframe(pd.DataFrame(rows, columns=["Setting", "Value"]), hide_index=True, use_container_width=True)
    if errors:
        st.error(errors[0])
        st.button("Return to Set up", on_click=_set_page, args=("Set up",), use_container_width=True)
        return
    if st.button("Run explicit TBI analysis", type="primary", use_container_width=True):
        with st.spinner("Running isolated explicit recent/remote epidemiology analysis..."):
            result = run_explicit_workflow_analysis(state)
        st.session_state[RESULT_KEY] = result
        st.session_state[STALE_KEY] = False
        st.success(f"Analysis completed in {result['elapsedSeconds']:.2f} seconds.")
    if st.session_state.get(RESULT_KEY):
        st.button("Open Results", on_click=_set_page, args=("Results",), type="primary", use_container_width=True)


def _results_page() -> None:
    state = _state()
    result = discard_incompatible_result(st.session_state.get(RESULT_KEY))
    if not result:
        st.info("Run the analysis to create explicit epidemiological results.")
        st.button("Open Run analysis", on_click=_set_page, args=("Run analysis",), use_container_width=True)
        return
    if st.session_state.get(STALE_KEY):
        st.warning("These results are stale because epidemiological inputs changed after the last run.")
    st.subheader("Decision-facing epidemiological results")
    st.dataframe(pd.DataFrame(result["summaryRows"]), hide_index=True, use_container_width=True)
    if result["mode"] == "deterministic":
        st.caption("Simulation intervals do not apply to a single deterministic expected-value calculation.")
    else:
        st.caption(f"Simulation intervals use {int(state['stochasticRepetitions']):,} simulated populations.")
    summary = pd.DataFrame(result["summaryRows"])
    averted = float(
        summary.loc[
            summary["Metric"] == "Direct active TB cases averted among the modelled population",
            summary.columns[1],
        ].iloc[0]
    )
    per_1000 = 1000.0 * averted / float(state["population"])
    st.info(f"Direct active TB cases averted among the modelled population: {averted:.2f}, or {per_1000:.2f} per 1,000 people modelled.")

    with st.expander("Outcomes by recent versus remote-only infection", expanded=False):
        if result["mode"] == "stochastic":
            rows = result["runnerResult"]["explicitEpidemiology"]["replicates"][0]["byEffectiveState"]
            st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
        else:
            det = result["deterministicResult"]["byHorizon"][f"{state['analysisHorizonYears']:g}"]
            st.dataframe(
                pd.DataFrame(
                    [
                        {
                            "State": "recent",
                            "Comparator active TB": det["recentComparatorActiveTBProportion"] * state["population"],
                            "Direct cases averted": det["recentDirectActiveTBCasesAvertedProportion"] * state["population"],
                        },
                        {
                            "State": "remote_only",
                            "Comparator active TB": det["remoteOnlyComparatorActiveTBProportion"] * state["population"],
                            "Direct cases averted": det["remoteOnlyDirectActiveTBCasesAvertedProportion"] * state["population"],
                        },
                    ]
                ),
                hide_index=True,
                use_container_width=True,
            )
    with st.expander("Outcomes by age", expanded=False):
        st.dataframe(pd.DataFrame(result["ageRows"]), hide_index=True, use_container_width=True)
    with st.expander("Timing relative to screening and effective TPT", expanded=False):
        if result["mode"] == "deterministic":
            det = result["deterministicResult"]["byHorizon"][f"{state['analysisHorizonYears']:g}"]
            timing = {
                "TB before screening": det["incidentActiveTBBeforeScreeningProportion"] * state["population"],
                "TB after screening but before effective TPT": det["incidentActiveTBAfterScreenBeforeEffectiveTPTProportion"] * state["population"],
                "Subsequent incident TB": det["subsequentIncidentActiveTBProportion"] * state["population"],
            }
        else:
            intervention = result["runnerResult"]["explicitEpidemiology"]["replicates"][0]["intervention"]
            timing = {
                "TB before screening": intervention["incidentActiveTBBeforeScreening"],
                "TB after screening but before effective TPT": intervention["incidentActiveTBAfterScreenBeforeEffectiveTPT"],
                "Subsequent incident TB": intervention["subsequentIncidentActiveTB"],
            }
        st.dataframe(pd.DataFrame(timing.items(), columns=["Timing", "Count"]), hide_index=True, use_container_width=True)
    with st.expander("Observed-versus-expected active-TB validation", expanded=False):
        validation = active_tb_validation_for_state(state)
        if not validation:
            st.info("No active-TB observation supplied.")
        else:
            st.write(validation["message"])
            st.dataframe(pd.DataFrame([validation]), hide_index=True, use_container_width=True)

    with st.expander("Calibration and age-specific display", expanded=False):
        calibration = result["calibration"]["explicitRecentRemoteCalibration"]
        st.dataframe(
            pd.DataFrame(
                [
                    {"Derived quantity": "Fitted recent exposure hazard", "Value": calibration["fittedRecentHazard"]},
                    {"Derived quantity": "Fitted remote exposure hazard", "Value": calibration["fittedRemoteHazard"]},
                    {"Derived quantity": "Achieved recent prevalence", "Value": calibration["achievedRecentPrevalenceTotalPopulation"]},
                    {"Derived quantity": "Achieved remote-only prevalence", "Value": calibration["achievedRemoteOnlyPrevalenceTotalPopulation"]},
                ]
            ),
            hide_index=True,
            use_container_width=True,
        )
    with st.expander("Downloads", expanded=False):
        st.download_button(
            "Download configuration JSON",
            data=json.dumps(result["config"], indent=2, sort_keys=True),
            file_name="explicit_tbi_configuration.json",
            mime="application/json",
        )
        st.download_button(
            "Download epidemiological results CSV",
            data=results_to_csv(result),
            file_name="explicit_tbi_epidemiological_results.csv",
            mime="text/csv",
        )
        st.download_button(
            "Download active-TB validation CSV",
            data=validation_to_csv(active_tb_validation_for_state(state)),
            file_name="explicit_tbi_active_tb_validation.csv",
            mime="text/csv",
        )
        st.download_button(
            "Download provenance and limitations text",
            data=provenance_text(result),
            file_name="explicit_tbi_provenance.txt",
            mime="text/plain",
        )
    st.button(
        "Continue to Evidence and technical information",
        on_click=_set_page,
        args=("Evidence and technical information",),
        type="primary",
        use_container_width=True,
    )


def _evidence_page() -> None:
    st.subheader("Definitions and model boundary")
    st.markdown(
        """
- `Recently infected within the past 5 years` means at least one infection in the five years immediately before baseline.
- `Remote infection only` means infection at least five years before baseline and no infection in the most recent five years.
- A recent reinfection uses the working reset policy for the central run; the no-reset sensitivity remains isolated.
- The progression curve is survivor-conditioned: future risk uses the incremental cumulative hazard from the current infection clock.
- The post-25-year hazard is an explicit extrapolation from the final evidence anchor.
- Preventive treatment is modelled as preventing only future progression after treatment becomes effective.
- Active-TB observations are validation-only in this workflow.
- The model reports direct active TB cases averted among the modelled population; it does not include endogenous transmission.
- Cost, disability-weight and quality-adjusted calculations are outside this workflow.
"""
    )
    st.info(
        "Specification: docs/explicit_recent_remote_tbi_milestone1_spec.md. "
        "Evidence dossier: docs/explicit_recent_remote_tbi_parameter_decision_dossier.md."
    )
    st.caption(
        f"Pathway: {EXPLICIT_RECENT_REMOTE_PATHWAY_ID}; TPT policy: "
        f"{TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION}."
    )


_init_state()

st.title("Explicit recent/remote TBI development workflow")
st.caption("Development and research workflow for isolated epidemiological review.")

page = st.sidebar.radio(
    "Workflow",
    WORKFLOW_PAGES,
    index=WORKFLOW_PAGES.index(st.session_state[PAGE_KEY]),
    key=PAGE_KEY,
)

if page == "Set up":
    _setup_page()
elif page == "Run analysis":
    _run_page()
elif page == "Results":
    _results_page()
else:
    _evidence_page()
