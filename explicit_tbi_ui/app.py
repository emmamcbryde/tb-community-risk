from __future__ import annotations

from pathlib import Path
import json
import re
import sys
from time import perf_counter
from typing import Any

import pandas as pd
import streamlit as st

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

from app.import_bootstrap import ensure_repo_root

ensure_repo_root(REPO_ROOT)

from app.explicit_tbi_workflow import (
    DETERMINISTIC_LABEL,
    EXPLICIT_RECENT_REMOTE_PATHWAY_ID,
    STOCHASTIC_LABEL,
    TPT_EFFECT_POLICY_PREVENTS_FUTURE_PROGRESSION,
    WORKFLOW_PAGES,
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
    technical_results_json,
    PROGRESSION_RISK_LABELS,
    progression_risk_for_state,
    validate_projection_accounting,
)
from engine.apy.population_tbi_progression_risk import ADJUSTMENT_FIELD


from app.explicit_tbi_formatting import (
    display_table, format_number, format_percent, percentage_input_format,
    format_reconciled_components,
    format_multiplier,
)
from app.explicit_tbi_session import (
    STATE_KEY, RESULT_KEY, STALE_KEY, PAGE_KEY, WIDGET_VERSION_KEY,
    DISPATCH_COUNT_KEY, RENDER_TIMING_KEY, INPUT_ERRORS_KEY, POPULATION_DRAFT_KEY,
    WIDGET_PREFIX, initialize_explicit_session,
)

RENDER_STARTED = perf_counter()


st.set_page_config(
    page_title="TB screening and preventive treatment model",
    layout="wide",
)


def _init_state() -> None:
    if initialize_explicit_session(st.session_state):
        st.caption("Incompatible previous session state was discarded. Set up a new explicit-TBI analysis.")
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
    st.session_state[INPUT_ERRORS_KEY] = {}
    st.session_state.pop(POPULATION_DRAFT_KEY, None)
    st.session_state[WIDGET_VERSION_KEY] = int(
        st.session_state.get(WIDGET_VERSION_KEY, 0)
    ) + 1


def _set_page(page: str) -> None:
    st.session_state[PAGE_KEY] = page
    st.session_state[_widget_key("workflow_navigation")] = page


def _navigate_from_widget(widget_key: str) -> None:
    st.session_state[PAGE_KEY] = st.session_state[widget_key]


def _state() -> dict[str, Any]:
    return st.session_state[STATE_KEY]


def _widget_key(name: str) -> str:
    return f"{WIDGET_PREFIX}{st.session_state[WIDGET_VERSION_KEY]}:{name}"


def _commit_percent(key: str, widget_key: str) -> None:
    _state()[key] = percent_to_proportion(st.session_state[widget_key])
    st.session_state[STALE_KEY] = True


def _commit_population(widget_key: str) -> None:
    text = st.session_state[widget_key].strip()
    st.session_state[POPULATION_DRAFT_KEY] = text
    errors = st.session_state[INPUT_ERRORS_KEY]
    if not re.fullmatch(r"(?:[0-9]+|[0-9]{1,3}(?:,[0-9]{3})+)", text):
        errors["population"] = "Enter a whole number of people, for example 1,500."
        return
    population = int(text.replace(",", ""))
    if not 1 <= population <= 1_000_000:
        errors["population"] = "Total modelled population must be between 1 and 1,000,000 people."
        return
    errors.pop("population", None)
    _state()["population"] = population
    st.session_state[POPULATION_DRAFT_KEY] = format_number(population, decimals=0)
    st.session_state[widget_key] = st.session_state[POPULATION_DRAFT_KEY]
    st.session_state[STALE_KEY] = True


def _percent_input(label: str, key: str, help_text: str = "", *, baseline: bool = False, source: bool = False, minimum: float = 0.0) -> None:
    state = _state()
    widget_key = _widget_key(key)
    st.number_input(
        label,
        min_value=minimum,
        max_value=100.0,
        value=proportion_to_percent(float(state[key])),
        step=0.1,
        format=percentage_input_format(state[key], baseline=baseline),
        help=help_text,
        key=widget_key,
        on_change=_commit_percent,
        args=(key, widget_key),
    )
    if source:
        st.caption(f"Source: {source_for_value(state, key)}")


def _commit_progression_risk_mode(widget_key: str) -> None:
    label = st.session_state[widget_key]
    _state()["progressionRiskMode"] = next(mode for mode, choice in PROGRESSION_RISK_LABELS.items() if choice == label)
    st.session_state[STALE_KEY] = True


def _commit_custom_multiplier(widget_key: str) -> None:
    _state()["customProgressionRiskMultiplier"] = st.session_state[widget_key]
    st.session_state[STALE_KEY] = True


def _progression_risk_control() -> None:
    state = _state()
    st.subheader("Progression risk among people with TBI")
    st.write("The modelled population may have a higher risk of progression to active TB than an otherwise comparable population with TBI. This adjustment multiplies the progression hazard; it does not change the entered TBI prevalence.")
    mode_key = _widget_key("progression_risk_mode")
    choices = list(PROGRESSION_RISK_LABELS.values())
    st.session_state[mode_key] = PROGRESSION_RISK_LABELS[state["progressionRiskMode"]]
    st.radio("Progression-risk scenario", choices,
             key=mode_key, on_change=_commit_progression_risk_mode, args=(mode_key,))
    if state["progressionRiskMode"] == "custom":
        custom_key = _widget_key("custom_progression_multiplier")
        st.number_input("Custom progression-hazard multiplier (×)",
                        value=float(state["customProgressionRiskMultiplier"]), step=0.1,
                        format="%g", key=custom_key,
                        on_change=_commit_custom_multiplier, args=(custom_key,))
    if state["progressionRiskMode"] != "standard":
        st.warning("This is a scenario assumption, not a calibrated or independently estimated relative risk.")


def _setup_page() -> None:
    state = _state()
    old_signature = health_signature(state)
    st.subheader("Population")
    cols = st.columns(3)
    with cols[0]:
        population_key = _widget_key("population")
        st.text_input(
            "Total modelled population (people)",
            value=st.session_state.get(POPULATION_DRAFT_KEY, format_number(state["population"], decimals=0)),
            key=population_key,
            on_change=_commit_population,
            args=(population_key,),
        )
    with cols[1]:
        state["analysisHorizonYears"] = float(
            st.number_input(
                "Analysis horizon",
                min_value=0.5,
                max_value=100.0,
                value=float(state["analysisHorizonYears"]),
                step=1.0,
                format="%g",
                key=_widget_key("analysis_horizon"),
            )
        )
        st.caption("Unit: years; fractional years are supported")
    with cols[2]:
        _percent_input(
            "Baseline active TB (% of total population)", "baselineActiveTBProportion",
            "Prevalent active TB at baseline is sequenced before TBI states.",
            baseline=True, source=True,
        )

    st.subheader("Input known TBI prevalence")
    st.caption("Unit: % of total population. Both inputs are percentages of the total modelled population.")
    tbi_cols = st.columns(2)
    with tbi_cols[0]:
        _percent_input(
            "Recently infected within the past 5 years",
            "recentTBITarget",
            "At least one TB infection acquired in the five years immediately before baseline.",
            source=True,
        )
    with tbi_cols[1]:
        _percent_input(
            "Remote infection only — infected 5 or more years ago and not reinfected within the past 5 years",
            "remoteOnlyTBITarget",
            "At least one infection more than five years before baseline and no infection in the most recent five years.",
            source=True,
        )
    implied_uninfected = max(
        0.0,
        1.0
        - float(state["baselineActiveTBProportion"])
        - float(state["recentTBITarget"])
        - float(state["remoteOnlyTBITarget"]),
    )
    st.metric("Implied uninfected prevalence", format_percent(implied_uninfected))
    _progression_risk_control()

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
            table = table.rename(columns={"Population share": "Population share (% of total population)"})
            st.dataframe(display_table(table, percentage_columns=["Population share (% of total population)"]), hide_index=True, use_container_width=True)
            st.caption("The open-ended 85+ source category uses the explicit age85PlusMax provenance from the runner contract.")
        except Exception as exc:
            st.warning(f"Age distribution cannot be displayed yet: {exc}")

    st.subheader("Intervention")
    int_cols = st.columns(3)
    with int_cols[0]:
        _percent_input("Screening coverage (%)", "screeningCoverage")
        state["screeningWindowYears"] = float(
            st.number_input(
                "Screening period (years)",
                min_value=0.01,
                max_value=20.0,
                value=float(state["screeningWindowYears"]),
                step=0.5,
                format="%g",
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
        _percent_input("Test sensitivity (%)", "testSensitivity")
        _percent_input("Test specificity (%)", "testSpecificity")
    with int_cols[2]:
        regimens = ["3HP", "4R", "3HR", "6H", "9H"]
        state["regimen"] = st.selectbox(
            "Preventive-treatment regimen",
            regimens,
            index=regimens.index(str(state["regimen"])),
            key=_widget_key("regimen"),
        )
        _percent_input("Treatment initiation (%)", "pStartTPT")
        _percent_input("Treatment completion (%)", "regimenPComplete")
        _percent_input("TPT efficacy after completion (%)", "regimenEffFull")

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
                        min_value=0,
                        value=int(state["activeTBObservedCases"]),
                        step=1,
                        key=_widget_key("active_tb_observed_cases"),
                    )
                )
                state["activeTBDenominator"] = float(
                    st.number_input(
                        "Observation denominator",
                        min_value=1,
                        value=int(state["activeTBDenominator"]),
                        step=100,
                        key=_widget_key("active_tb_denominator"),
                    )
                )
            with obs_cols[1]:
                state["activeTBObservationDurationYears"] = float(
                    st.number_input(
                        "Observation duration (years)",
                        min_value=0.01,
                        value=float(state["activeTBObservationDurationYears"]),
                        step=0.5,
                        format="%g",
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
                _percent_input("Ascertainment probability (%)", "activeTBAscertainmentProbability", minimum=0.000001)
                state["activeTBObservationSource"] = st.text_input(
                    "Observation source",
                    value=str(state["activeTBObservationSource"]),
                    key=_widget_key("active_tb_source"),
                )

    errors = list(st.session_state[INPUT_ERRORS_KEY].values()) or validate_workflow_state(state)
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
    errors = list(st.session_state[INPUT_ERRORS_KEY].values()) or validate_workflow_state(state, check_calibration=False)
    st.subheader("Pre-run summary")
    implied_uninfected = (
        1.0
        - float(state["baselineActiveTBProportion"])
        - float(state["recentTBITarget"])
        - float(state["remoteOnlyTBITarget"])
    )
    rows = [
        ("Population", format_number(state['population'], decimals=0)),
        ("Recent TBI target (% of total population)", format_percent(state['recentTBITarget'])),
        ("Remote-only TBI target (% of total population)", format_percent(state['remoteOnlyTBITarget'])),
        ("Implied uninfected (% of total population)", format_percent(implied_uninfected)),
        ("Baseline active TB (% of total population)", format_percent(state['baselineActiveTBProportion'], decimals=2, trim=True)),
        ("Test and regimen", f"{state['testType']} / {state['regimen']}"),
        ("Coverage", format_percent(state['screeningCoverage'])),
        ("Horizon", f"{format_number(state['analysisHorizonYears'])} years"),
        ("Analysis method", state["analysisMethod"]),
    ]
    try:
        risk = progression_risk_for_state(state)
        rows.append(("Progression risk", f"{PROGRESSION_RISK_LABELS[risk['mode']]} ({format_multiplier(risk['multiplier'])} hazard)"))
    except ValueError:
        pass  # The validation error below blocks dispatch without clipping.
    if state["analysisMethod"] == STOCHASTIC_LABEL:
        rows.extend(
            [
                ("Repetitions", format_number(state['stochasticRepetitions'], decimals=0)),
                ("Random seed", str(int(state["seed"]))),
            ]
        )
    st.dataframe(pd.DataFrame(rows, columns=["Setting", "Value"]), hide_index=True, use_container_width=True)
    if errors:
        st.error(errors[0])
        st.button("Return to Set up", on_click=_set_page, args=("Set up",), use_container_width=True)
        return
    if st.button("Run explicit TBI analysis", type="primary", use_container_width=True):
        st.session_state[DISPATCH_COUNT_KEY] += 1
        analysis_started = perf_counter()
        with st.spinner("Running isolated explicit recent/remote epidemiology analysis..."):
            try:
                result = run_explicit_workflow_analysis(state)
            except ValueError as exc:
                st.error(str(exc))
                return
        result["timings"]["streamlitDispatchSeconds"] = perf_counter() - analysis_started
        st.session_state[RESULT_KEY] = result
        st.session_state[STALE_KEY] = False
        st.success(f"Analysis completed in {format_number(result['elapsedSeconds'])} seconds.")
    if st.session_state.get(RESULT_KEY):
        st.button("Open Results", on_click=_set_page, args=("Results",), type="primary", use_container_width=True)


def _results_page() -> None:
    result = discard_incompatible_result(st.session_state.get(RESULT_KEY))
    if not result:
        st.info("Run the analysis to create explicit epidemiological results.")
        st.button("Open Run analysis", on_click=_set_page, args=("Run analysis",), use_container_width=True)
        return
    if st.session_state.get(STALE_KEY):
        st.warning("These results are stale because epidemiological inputs changed after the last run.")
    # Describe the stored run, even if the current inputs have since been edited.
    config = result["config"]
    population = int(config["N"])
    horizon = float(config["followUpHorizonYears"])
    try:
        checked = validate_projection_accounting(config, result["analyticalProjection"], result["comparatorDecomposition"])
        if not checked["isValid"]:
            st.error("Result accounting failed. Run the analysis again: " + "; ".join(checked["errors"]))
            return
    except (KeyError, TypeError, ValueError) as exc:
        st.error(f"Result accounting cannot be verified. Run the analysis again: {exc}")
        return
    st.subheader("Epidemiological results")
    risk = result[ADJUSTMENT_FIELD]
    st.caption(f"Progression risk: {PROGRESSION_RISK_LABELS[risk['mode']]}; resolved hazard multiplier {format_multiplier(risk['multiplier'])}.")
    st.info("No future infections or endogenous transmission are modelled.")
    st.metric(f"Expected future active TB cases over {format_number(horizon)} years among the baseline population",
              format_number(result["comparatorDecomposition"][-1]["expectedComparatorActiveTBCases"]))
    if risk["mode"] != "standard":
        st.caption("The progression-risk adjustment is a provisional scenario assumption, not a calibrated or independently estimated relative risk.")
    st.caption("Counts are people or direct active TB cases over the analysis horizon. Values are approximate when fractional.")
    st.dataframe(display_table(result["summaryRows"]), hide_index=True, use_container_width=True)
    if result["mode"] == "deterministic":
        st.caption("Simulation intervals do not apply to a single deterministic expected-value calculation.")
    else:
        st.caption(f"Simulation intervals use {format_number(config['nReps'], decimals=0)} simulated populations; Mean and low/high 95% simulation intervals are shown.")
    st.subheader("Where expected comparator cases come from")
    st.caption("Analytical expected values for the completed configuration. Baseline active TB is separate and excluded from future incident cases. Values are approximate; totals use full precision.")
    if result["mode"] == "stochastic":
        st.caption("This analytical decomposition and sensitivity comparison are not additional stochastic simulations. Simulation means above may differ.")
    decomposition = result["comparatorDecomposition"]
    rounded = format_reconciled_components([decomposition[i]["expectedComparatorActiveTBCases"] for i in (0, 1, 2)], decomposition[-1]["expectedComparatorActiveTBCases"])
    display_rows = []
    for index, row in enumerate(decomposition):
        display_rows.append({
            "Baseline group": row["baselineGroup"],
            "Number of people": format_number(row["numberOfPeople"]),
            f"Average {format_number(horizon)}-year progression probability": "Not applicable" if index == 3 else ("—" if index == 4 else format_percent(row["averageProgressionProbability"])),
            "Expected comparator active-TB cases": "Excluded from incident projection" if index == 3 else (format_number(row["expectedComparatorActiveTBCases"]) if index == 4 else rounded[index]),
        })
    st.dataframe(pd.DataFrame(display_rows), hide_index=True, use_container_width=True)
    st.subheader("Progression-risk sensitivity comparison")
    st.caption("Analytical expected values: the same population, TBI prevalence, age distribution and intervention. Only the progression-risk adjustment differs.")
    sensitivity = display_table(result["progressionRiskSensitivity"]["rows"])
    sensitivity["Mode"] = [PROGRESSION_RISK_LABELS[row["Mode"]].split(" — ")[0] for row in result["progressionRiskSensitivity"]["rows"]]
    sensitivity["Multiplier"] = [format_multiplier(row['Multiplier']) for row in result["progressionRiskSensitivity"]["rows"]]
    st.dataframe(sensitivity, hide_index=True, use_container_width=True)
    summary = pd.DataFrame(result["summaryRows"])
    value_column = "Expected value" if result["mode"] == "deterministic" else "Mean"
    averted = float(
        summary.loc[
            summary["Metric"] == "Direct active TB cases averted among the modelled population",
            value_column,
        ].iloc[0]
    )
    per_1000 = 1000.0 * averted / population
    qualifier = "expected" if result["mode"] == "deterministic" else "mean"
    st.info(f"Approximately {format_number(averted)} {qualifier} direct active TB cases averted among the modelled population, or {format_number(per_1000)} per 1,000 people modelled.")

    with st.expander("Outcomes by recent versus remote-only infection", expanded=False):
        if result["mode"] == "stochastic":
            rows = result["runnerResult"]["explicitEpidemiology"]["replicates"][0]["byEffectiveState"]
            st.caption("First simulated population: individual counts of people or active TB cases.")
            st.dataframe(display_table(rows), hide_index=True, use_container_width=True)
        else:
            det = result["deterministicResult"]["byHorizon"][f"{horizon:g}"]
            st.dataframe(
                display_table(
                    [
                        {
                            "State": "recent",
                            "Comparator active TB (expected cases)": det["recentComparatorActiveTBProportion"] * population,
                            "Direct cases averted (expected cases)": det["recentDirectActiveTBCasesAvertedProportion"] * population,
                        },
                        {
                            "State": "remote_only",
                            "Comparator active TB (expected cases)": det["remoteOnlyComparatorActiveTBProportion"] * population,
                            "Direct cases averted (expected cases)": det["remoteOnlyDirectActiveTBCasesAvertedProportion"] * population,
                        },
                    ]
                ),
                hide_index=True,
                use_container_width=True,
            )
    with st.expander("Outcomes by age", expanded=False):
        ages = pd.DataFrame(result["ageRows"])
        percent_columns = [column for column in ages if column != "Age group"]
        ages = display_table(ages, percentage_columns=percent_columns)
        ages = ages.rename(columns={column: column + " (%)" for column in percent_columns})
        st.caption(f"Within-age-group percentages; population share is % of total population. Active-TB risk is over {format_number(horizon)} years.")
        st.dataframe(ages, hide_index=True, use_container_width=True)
    with st.expander("Timing relative to screening and effective TPT", expanded=False):
        if result["mode"] == "deterministic":
            det = result["deterministicResult"]["byHorizon"][f"{horizon:g}"]
            timing = {
                "TB before screening": det["incidentActiveTBBeforeScreeningProportion"] * population,
                "TB after screening but before effective TPT": det["incidentActiveTBAfterScreenBeforeEffectiveTPTProportion"] * population,
                "Subsequent incident TB": det["subsequentIncidentActiveTBProportion"] * population,
            }
        else:
            st.caption("First simulated population: individual active TB case counts.")
            intervention = result["runnerResult"]["explicitEpidemiology"]["replicates"][0]["intervention"]
            timing = {
                "TB before screening": intervention["incidentActiveTBBeforeScreening"],
                "TB after screening but before effective TPT": intervention["incidentActiveTBAfterScreenBeforeEffectiveTPT"],
                "Subsequent incident TB": intervention["subsequentIncidentActiveTB"],
            }
        timing_table = pd.DataFrame(timing.items(), columns=["Timing", "Expected cases" if result["mode"] == "deterministic" else "Cases in first simulated population"])
        st.dataframe(display_table(timing_table), hide_index=True, use_container_width=True)
    with st.expander("Observed-versus-expected active-TB validation", expanded=False):
        validation = result.get("activeTBValidation")
        if not validation:
            st.info("No active-TB observation supplied.")
        else:
            st.write(validation["message"])
            validation_rows = [{
                "Observed cases": validation["observedCases"],
                "Expected cases": validation["expectedCases"],
                "Source population (people)": validation["sourcePopulationDenominator"],
                "Observation duration (years)": validation["observationDurationYears"],
                "Observed (% of source population)": validation["observedRate"],
                "Expected (% of source population)": validation["expectedRate"],
                "Difference (cases)": validation["absoluteDifference"],
                "Observed / expected": validation["ratioObservedToExpected"] if validation["ratioObservedToExpected"] is not None else "N/A",
            }]
            st.dataframe(display_table(validation_rows, percentage_columns=["Observed (% of source population)", "Expected (% of source population)"]), hide_index=True, use_container_width=True)
    with st.expander("Downloads", expanded=False):
        st.caption("Machine-readable downloads retain full precision.")
        st.download_button("Download full technical results JSON", data=technical_results_json(result), file_name="explicit_tbi_full_results.json", mime="application/json")
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
            data=validation_to_csv(result.get("activeTBValidation")),
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
    result = discard_incompatible_result(st.session_state.get(RESULT_KEY))
    if result:
        st.subheader("Technical diagnostics")
        st.caption("Fitted hazards are annual exposure hazards per person-year. Technical values retain full precision.")
        calibration = result["calibration"]["explicitRecentRemoteCalibration"]
        st.json(calibration)
        st.caption(f"Configuration hash: {result['configHash']}")
        st.json({"analysisTimingsSeconds": result.get("timings"), "previousScriptRenderSeconds": st.session_state.get(RENDER_TIMING_KEY), "analysisDispatchCount": st.session_state[DISPATCH_COUNT_KEY]})
        st.json(result["eventLedger"])
        st.json(result[ADJUSTMENT_FIELD])


_init_state()

st.title("TB screening and preventive treatment model")
st.caption("Development status: explicit-TBI epidemiology for planning and sequencing; this model is not for denying care.")

navigation_key = _widget_key("workflow_navigation")
st.session_state[navigation_key] = st.session_state[PAGE_KEY]
st.sidebar.radio(
    "Workflow",
    WORKFLOW_PAGES,
    key=navigation_key,
    on_change=_navigate_from_widget,
    args=(navigation_key,),
)
page = st.session_state[PAGE_KEY]

if page == "Set up":
    _setup_page()
elif page == "Run analysis":
    _run_page()
elif page == "Results":
    _results_page()
else:
    _evidence_page()

st.session_state[RENDER_TIMING_KEY] = perf_counter() - RENDER_STARTED
