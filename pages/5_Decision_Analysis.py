from __future__ import annotations

import pandas as pd
import numpy as np
import streamlit as st

from app.display import arrow_safe_dataframe, safe_download_stem
from app.icon_arrays import build_100_person_visual_data, render_100_person_summary
from app.state import init_session_state
from engine.apy.decision_analysis import run_scenario_comparison
from engine.apy.early_review import run_early_screening_review
from engine.apy.evidence import assess_apy_reference_readiness
from engine.apy.scenario import DIRECT_EFFECTS_SCOPE_STATEMENT
from engine.apy.sensitivity import (
    load_sensitivity_specs,
    run_one_way_sensitivity,
    run_threshold_analysis,
)


init_session_state()

st.title("Explore Decisions")
st.caption("Compare strategies, explore supplied sensitivity ranges and review early screening results.")

MODEL_TYPE_LABELS = {
    "expected_value": "Expected outcomes",
    "agent_based": "Simulated community variation",
}
STRATEGY_LABELS = {
    "prevent": "Prioritise people most likely to avoid active TB",
    "ltbi": "Prioritise people most likely to have LTBI",
    "cure": "Prioritise people most likely to complete effective treatment",
    "random": "No risk-based prioritisation",
}
STRATEGY_VALUES = {label: value for value, label in STRATEGY_LABELS.items()}

config = st.session_state.get("config")
economics_config = st.session_state.get("economics_config")
results_bundle = st.session_state.get("results_bundle")

if not config:
    st.info("Create an analysis before exploring decisions.")
    st.stop()

readiness = assess_apy_reference_readiness(config, economics_config or {})
ledger = (results_bundle or {}).get("technical", {}).get("eventLedger", {}) if isinstance(results_bundle, dict) else {}
ledger_valid = bool(ledger and ledger.get("validation", {}).get("isValid") is True)
econ_results = st.session_state.get("economics_results")
econ_valid = bool(econ_results and econ_results.get("validation", {}).get("structurallyValid"))

overview_rows = [
    {"field": "population", "value": config.get("populationPresetId")},
    {"field": "comparator", "value": "current practice / no additional systematic LTBI screening"},
    {"field": "intervention", "value": "targeted LTBI screening and preventive treatment"},
    {"field": "direct-effects scope", "value": DIRECT_EFFECTS_SCOPE_STATEMENT},
    {"field": "overall clinician-ready", "value": readiness.get("overallClinicianReady")},
    {"field": "event ledger valid", "value": ledger_valid},
    {"field": "economics structurally valid", "value": econ_valid},
]
st.dataframe(arrow_safe_dataframe(overview_rows), use_container_width=True, hide_index=True)
if not readiness.get("overallClinicianReady"):
    st.warning("Reference evidence remains unresolved or provisional. The page reports modelled consequences only, not recommendations.")

tab_compare, tab_sensitivity, tab_early = st.tabs(
    ["Compare strategies", "Explore sensitivity", "Review early screening results"]
)

with tab_compare:
    st.subheader("Compare Strategies")
    model_type = st.radio(
        "Analysis type",
        ["expected_value", "agent_based"],
        horizontal=True,
        format_func=lambda value: MODEL_TYPE_LABELS.get(value, str(value)),
    )
    col_a, col_b = st.columns(2)
    with col_a:
        base_test = st.selectbox("Base test", ["IGRA", "TST"], index=0)
        base_regimen = st.selectbox("Base regimen", ["3HP", "4R", "3HR", "6H", "9H"], index=0)
        base_strategy_label = st.selectbox("Base prioritisation", list(STRATEGY_VALUES), index=0)
        base_coverage = st.number_input("Base coverage", min_value=0.0, max_value=1.0, value=float(config.get("screenCoverage", 0.3)), step=0.05)
    with col_b:
        alt_test = st.selectbox("Alternative test", ["IGRA", "TST"], index=1)
        alt_regimen = st.selectbox("Alternative regimen", ["3HP", "4R", "3HR", "6H", "9H"], index=0)
        alt_strategy_label = st.selectbox("Alternative prioritisation", list(STRATEGY_VALUES), index=1)
        alt_coverage = st.number_input("Alternative coverage", min_value=0.0, max_value=1.0, value=float(config.get("screenCoverage", 0.3)), step=0.05)
    n_reps = st.number_input("Simulation replicates", min_value=1, value=min(int(config.get("nReps", 10)), 20), step=1)
    if st.button("Run strategy comparison", type="primary"):
        scenarios = [
            {
                "scenarioId": "base_strategy",
                "label": "Base strategy",
                "changes": {
                    "test": base_test,
                    "regimen": base_regimen,
                    "screeningStrategy": STRATEGY_VALUES[base_strategy_label],
                    "screenCoverage": base_coverage,
                },
            },
            {
                "scenarioId": "alternative_strategy",
                "label": "Alternative strategy",
                "changes": {
                    "test": alt_test,
                    "regimen": alt_regimen,
                    "screeningStrategy": STRATEGY_VALUES[alt_strategy_label],
                    "screenCoverage": alt_coverage,
                },
            },
        ]
        with st.spinner("Running strategy comparison..."):
            st.session_state["decision_scenario_comparison"] = run_scenario_comparison(
                config,
                economics_config,
                scenarios,
                model_type=model_type,
                n_reps=int(n_reps) if model_type == "agent_based" else None,
                master_seed=int(config.get("seed", 1)),
            )
    comparison = st.session_state.get("decision_scenario_comparison")
    if comparison:
        st.dataframe(arrow_safe_dataframe(comparison["scenarioSummaries"]), use_container_width=True, hide_index=True)
        scenarios = comparison.get("scenarios") or []
        if scenarios:
            visual_scenarios = scenarios[:2]
            st.subheader("100-person strategy summaries")
            st.caption("Each strategy shows outcomes per 100 eligible people using its own event ledger.")
            visual_cols = st.columns(len(visual_scenarios))
            for col, scenario in zip(visual_cols, visual_scenarios):
                with col:
                    visual_rows = build_100_person_visual_data(scenario.get("eventLedger"))
                    if visual_rows:
                        render_100_person_summary(
                            visual_rows,
                            title=str(scenario.get("label", "Selected strategy")),
                        )
        if comparison.get("pairedComparisons"):
            st.dataframe(arrow_safe_dataframe(comparison["pairedComparisons"]), use_container_width=True, hide_index=True)
        if comparison.get("pairedReplicateComparisons"):
            with st.expander("Paired replicate diagnostics"):
                st.dataframe(arrow_safe_dataframe(comparison["pairedReplicateComparisons"]), use_container_width=True, hide_index=True)
                st.dataframe(arrow_safe_dataframe(comparison.get("pairedDifferenceSummaries", [])), use_container_width=True, hide_index=True)
        if comparison.get("commonSeedNonpairedDiagnostics"):
            with st.expander("Common-seed non-paired diagnostics"):
                st.caption("Baseline-generating inputs differ; these rows are not paired clinical effect intervals.")
                st.dataframe(arrow_safe_dataframe(comparison["commonSeedNonpairedDiagnostics"]), use_container_width=True, hide_index=True)
        st.caption("Expected outcomes may be fractional. Simulation percentiles describe finite-population variation, not confidence intervals.")
        st.download_button(
            "Download strategy comparison CSV",
            data=pd.DataFrame(comparison["scenarioSummaries"]).to_csv(index=False).encode("utf-8"),
            file_name=f"{safe_download_stem(config.get('scenarioLabel'), 'scenario_comparison')}.csv",
            mime="text/csv",
        )

with tab_sensitivity:
    st.subheader("One-Way Sensitivity")
    specs = load_sensitivity_specs()
    ready_specs = [spec for spec in specs if spec.get("lowValue") is not None and spec.get("highValue") is not None]
    unresolved_specs = [spec for spec in specs if spec.get("lowValue") is None or spec.get("highValue") is None]
    if unresolved_specs:
        st.info("APY reference sensitivity ranges are unresolved; no ranges are invented.")
        st.dataframe(
            arrow_safe_dataframe(
                [
                    {
                        "parameterId": spec["parameterId"],
                        "status": spec["reviewStatus"],
                        "unresolvedReason": spec["unresolvedReason"],
                    }
                    for spec in unresolved_specs
                ]
            ),
            use_container_width=True,
            hide_index=True,
        )
    if ready_specs and st.button("Run one-way sensitivity"):
        outcomes = ["active_tb_cases_prevented", "infection_effectively_treated_total", "incrementalCost", "dalysAverted", "nmb"]
        st.session_state["decision_sensitivity"] = run_one_way_sensitivity(config, economics_config, ready_specs, outcomes)
    sensitivity = st.session_state.get("decision_sensitivity")
    if sensitivity:
        st.dataframe(arrow_safe_dataframe(sensitivity["results"]), use_container_width=True, hide_index=True)
        st.download_button(
            "Download sensitivity results CSV",
            data=pd.DataFrame(sensitivity["results"]).to_csv(index=False).encode("utf-8"),
            file_name=f"{safe_download_stem(config.get('scenarioLabel'), 'sensitivity')}.csv",
            mime="text/csv",
        )

    st.subheader("Threshold Analysis")
    threshold_metric = st.selectbox("Decision metric", ["active_tb_cases_prevented", "nmb"])
    parameter = st.selectbox("Parameter", ["ltbiPrevalence", "pStartTPT"])
    low = st.number_input("Search low", min_value=0.0, max_value=1.0, value=0.01, step=0.01)
    high = st.number_input("Search high", min_value=0.0, max_value=1.0, value=0.2, step=0.01)
    target = st.number_input("Target metric value", value=0.0, step=0.1)
    if st.button("Run threshold analysis"):
        st.session_state["decision_threshold"] = run_threshold_analysis(
            config,
            economics_config,
            {"parameterId": parameter, "adapter": parameter, "label": parameter},
            threshold_metric,
            {"low": low, "high": high, "target": target, "gridPoints": 7},
        )
    threshold = st.session_state.get("decision_threshold")
    if threshold:
        if not threshold["validation"]["isValid"]:
            st.warning("; ".join(item["message"] for item in threshold["validation"].get("warnings", [])))
        st.caption(f"Monotonicity: {threshold.get('monotonicity')}")
        st.dataframe(arrow_safe_dataframe(threshold.get("grid", [])), use_container_width=True, hide_index=True)
        st.dataframe(arrow_safe_dataframe(threshold.get("crossings", [])), use_container_width=True, hide_index=True)

with tab_early:
    st.subheader("Review Early Screening Results")
    screened = st.number_input("Number screened to date", min_value=0, value=0, step=1)
    positives = st.number_input("Positive tests to date", min_value=0, value=0, step=1)
    planned = st.number_input("Planned total screened", min_value=0, value=int(round(float(config.get("screenCoverage", 0.3)) * int(config.get("N", 0)))), step=1)
    review_time = st.number_input("Review time", min_value=0.0, value=0.0, step=0.25)
    st.markdown("Explicit prior")
    prior_mean = st.number_input("Prior mean LTBI prevalence", min_value=0.001, max_value=0.999, value=0.1, step=0.01)
    prior_ess = st.number_input("Prior effective sample size", min_value=0.001, value=20.0, step=1.0)
    grid_low = st.number_input("Prevalence grid low", min_value=0.0, max_value=1.0, value=0.01, step=0.01)
    grid_high = st.number_input("Prevalence grid high", min_value=0.0, max_value=1.0, value=0.2, step=0.01)
    if st.button("Run early review"):
        grid = [float(x) for x in np.linspace(grid_low, grid_high, 9)]
        st.session_state["decision_early_review"] = run_early_screening_review(
            config,
            economics_config,
            {
                "reviewId": "streamlit_early_review",
                "scenarioId": config.get("scenarioLabel"),
                "screenedToDate": screened,
                "testPositiveToDate": positives,
                "plannedTotalScreened": planned,
                "eligiblePopulation": config.get("N"),
                "reviewTimeYears": review_time,
                "prior": {
                    "type": "beta",
                    "mean": prior_mean,
                    "effectiveSampleSize": prior_ess,
                    "source": "User-supplied Streamlit prior",
                },
                "prevalenceGrid": grid,
                "prevalenceThresholds": [prior_mean],
            },
        )
    early = st.session_state.get("decision_early_review")
    if early:
        if not early["validation"]["isValid"]:
            st.error(early["validation"]["errors"])
        else:
            st.dataframe(arrow_safe_dataframe([early["prior"]["summary"], early["posterior"]["summary"]]), use_container_width=True, hide_index=True)
            with st.expander("Technical information", expanded=False):
                st.caption(f"Calibration policy: {early.get('calibrationPolicy')}")
                st.caption(f"Reference calibration: {early.get('referenceCalibrationHash')}")
                st.caption(f"Prior discretisation: {(early.get('prior') or {}).get('discretisationMethod')}")
            st.caption(early.get("likelihoodNotes", ""))
            st.info(early.get("activeTBSurveillanceJointUpdateNotes", ""))
            if early.get("timingApproximation"):
                st.warning(early.get("timingApproximationDescription"))
            if early.get("economicTimingStatus") == "approximate_not_decision_grade":
                st.warning("Approximate economic components - timing not fully represented. NMB and continuation conclusions are unavailable.")
            st.dataframe(arrow_safe_dataframe(early["posteriorProjectionSummary"]), use_container_width=True, hide_index=True)
            st.download_button(
                "Download early-review posterior CSV",
                data=pd.DataFrame(early["priorPosteriorTable"]).to_csv(index=False).encode("utf-8"),
                file_name=f"{safe_download_stem(config.get('scenarioLabel'), 'early_review_posterior')}.csv",
                mime="text/csv",
            )
