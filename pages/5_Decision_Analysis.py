from __future__ import annotations

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from app.decision_comparison_presentation import (
    common_comparator_status,
    comparison_plane_rows,
    comparison_table_rows,
    head_to_head_row,
    strategy_display_name,
)
from app.display import arrow_safe_dataframe, safe_download_stem
from app.icon_arrays import build_100_person_visual_data, render_100_person_summary
from app.state import (
    ensure_frozen_reference_loaded_if_eligible,
    init_session_state,
    sanitize_reference_only_state,
)
from engine.apy.decision_analysis import run_scenario_comparison
from engine.apy.early_review import run_early_screening_review
from engine.apy.evidence import assess_apy_reference_readiness
from engine.apy.sensitivity import (
    load_sensitivity_specs,
    run_one_way_sensitivity,
    run_threshold_analysis,
)


MODEL_TYPE_LABELS = {
    "expected_value": "Quick deterministic preview - single expected-value calculation",
    "agent_based": "SA Health report analysis - simulated communities",
}
STRATEGY_LABELS = {
    "prevent": "Prioritise people most likely to avoid active TB",
    "ltbi": "Prioritise people most likely to have LTBI",
    "cure": "Prioritise people most likely to complete effective treatment",
    "random": "No risk-based prioritisation",
}
STRATEGY_VALUES = {label: value for value, label in STRATEGY_LABELS.items()}
SHORT_STRATEGY_LABELS = {
    "prevent": "prevent targeting",
    "ltbi": "LTBI targeting",
    "cure": "completion targeting",
    "random": "random targeting",
}


def _comparison_chart(rows: list[dict]) -> alt.Chart:
    frame = pd.DataFrame(rows)
    x_col = "DALYs averted compared with business as usual"
    y_col = "Incremental cost compared with business as usual (AUD)"
    points = (
        alt.Chart(frame)
        .mark_point(filled=True, size=130)
        .encode(
            x=alt.X(f"{x_col}:Q", title=x_col),
            y=alt.Y(f"{y_col}:Q", title=y_col),
            color=alt.Color("Strategy:N", legend=alt.Legend(title="Strategy")),
            shape=alt.Shape("Classification:N", legend=alt.Legend(title="Classification")),
            tooltip=[
                alt.Tooltip("Strategy:N"),
                alt.Tooltip(f"{x_col}:Q", format=".2f"),
                alt.Tooltip(f"{y_col}:Q", format=",.0f"),
                alt.Tooltip("Arithmetic ICER:N"),
                alt.Tooltip("Classification:N"),
                alt.Tooltip("Quadrant:N"),
            ],
        )
    )
    labels = (
        alt.Chart(frame)
        .mark_text(align="left", dx=8, dy=-6, fontSize=12)
        .encode(x=f"{x_col}:Q", y=f"{y_col}:Q", text="Strategy:N")
    )
    horizontal = alt.Chart(pd.DataFrame({y_col: [0]})).mark_rule(color="#6b7280").encode(y=f"{y_col}:Q")
    vertical = alt.Chart(pd.DataFrame({x_col: [0]})).mark_rule(color="#6b7280").encode(x=f"{x_col}:Q")
    return (horizontal + vertical + points + labels).properties(height=380)


def _run_strategy_comparison(
    *,
    config: dict,
    economics_config: dict | None,
    model_type: str,
    n_reps: int,
    strategy_a: dict,
    strategy_b: dict,
) -> dict:
    scenarios = [
        {
            "scenarioId": "strategy_a",
            "label": strategy_a["label"],
            "changes": strategy_a["changes"],
        },
        {
            "scenarioId": "strategy_b",
            "label": strategy_b["label"],
            "changes": strategy_b["changes"],
        },
    ]
    return run_scenario_comparison(
        config,
        economics_config,
        scenarios,
        model_type=model_type,
        n_reps=int(n_reps) if model_type == "agent_based" else None,
        master_seed=int(config.get("seed", 1)),
    )


init_session_state()
sanitize_reference_only_state()
ensure_frozen_reference_loaded_if_eligible()

st.title("Explore Decisions")
st.caption(
    "Compare two screening strategies using the same business-as-usual population, "
    "current economic assumptions and fixed SA Health report assumptions for future TB."
)

config = st.session_state.get("config")
economics_config = st.session_state.get("economics_config")

if not config:
    st.info("Create an analysis before exploring decisions.")
    st.stop()

readiness = assess_apy_reference_readiness(config, economics_config or {})
if not readiness.get("overallClinicianReady"):
    st.warning("Some evidence inputs remain provisional. Compare modelled consequences, not recommendations.")

st.subheader("Strategy comparison")
model_type = st.radio(
    "Analysis type",
    ["expected_value", "agent_based"],
    horizontal=True,
    format_func=lambda value: MODEL_TYPE_LABELS.get(value, str(value)),
)

col_a, col_b = st.columns(2)
with col_a:
    st.markdown("Strategy A")
    base_test = st.selectbox("Test A", ["IGRA", "TST"], index=0)
    base_regimen = st.selectbox("Regimen A", ["3HP", "4R", "3HR", "6H", "9H"], index=0)
    base_strategy_label = st.selectbox("Prioritisation A", list(STRATEGY_VALUES), index=0)
    base_coverage = st.number_input(
        "Coverage A",
        min_value=0.0,
        max_value=1.0,
        value=float(config.get("screenCoverage", 0.3)),
        step=0.05,
    )
with col_b:
    st.markdown("Strategy B")
    alt_test = st.selectbox("Test B", ["IGRA", "TST"], index=1)
    alt_regimen = st.selectbox("Regimen B", ["3HP", "4R", "3HR", "6H", "9H"], index=0)
    alt_strategy_label = st.selectbox("Prioritisation B", list(STRATEGY_VALUES), index=1)
    alt_coverage = st.number_input(
        "Coverage B",
        min_value=0.0,
        max_value=1.0,
        value=float(config.get("screenCoverage", 0.3)),
        step=0.05,
    )

n_reps = st.number_input("Simulation repetitions", min_value=1, value=min(int(config.get("nReps", 10)), 20), step=1)
strategy_a = {
    "label": strategy_display_name(
        test=base_test,
        regimen=base_regimen,
        coverage=base_coverage,
        prioritisation=SHORT_STRATEGY_LABELS[STRATEGY_VALUES[base_strategy_label]],
    ),
    "changes": {
        "test": base_test,
        "regimen": base_regimen,
        "screeningStrategy": STRATEGY_VALUES[base_strategy_label],
        "screenCoverage": base_coverage,
    },
}
strategy_b = {
    "label": strategy_display_name(
        test=alt_test,
        regimen=alt_regimen,
        coverage=alt_coverage,
        prioritisation=SHORT_STRATEGY_LABELS[STRATEGY_VALUES[alt_strategy_label]],
    ),
    "changes": {
        "test": alt_test,
        "regimen": alt_regimen,
        "screeningStrategy": STRATEGY_VALUES[alt_strategy_label],
        "screenCoverage": alt_coverage,
    },
}

if st.button("Run strategy comparison", type="primary"):
    with st.spinner("Running strategy comparison..."):
        st.session_state["decision_scenario_comparison"] = _run_strategy_comparison(
            config=config,
            economics_config=economics_config,
            model_type=model_type,
            n_reps=int(n_reps),
            strategy_a=strategy_a,
            strategy_b=strategy_b,
        )

comparison = st.session_state.get("decision_scenario_comparison")
if comparison:
    status = common_comparator_status(comparison)
    if not status["valid"]:
        st.warning(status["message"])
    else:
        st.caption(status["message"])
        display_rows = comparison_table_rows(comparison)
        st.dataframe(arrow_safe_dataframe(display_rows), use_container_width=True, hide_index=True)
        plane_rows = comparison_plane_rows(comparison)
        st.altair_chart(_comparison_chart(plane_rows), use_container_width=True)
        st.caption(
            "Business as usual is the shared comparator at the origin. Strategy points show incremental cost "
            "and DALYs averted compared with that comparator."
        )
        pairwise = head_to_head_row(comparison)
        if pairwise:
            st.markdown("Strategy B compared with Strategy A")
            st.dataframe(arrow_safe_dataframe([pairwise]), use_container_width=True, hide_index=True)
        scenarios = comparison.get("scenarios") or []
        if scenarios:
            with st.expander("100-person strategy summaries", expanded=False):
                visual_cols = st.columns(min(len(scenarios), 2))
                for col, scenario in zip(visual_cols, scenarios[:2]):
                    with col:
                        visual_rows = build_100_person_visual_data(scenario.get("eventLedger"))
                        if visual_rows:
                            render_100_person_summary(
                                visual_rows,
                                title=str(scenario.get("label", "Selected strategy")),
                            )
        st.download_button(
            "Download strategy comparison CSV",
            data=pd.DataFrame(display_rows).to_csv(index=False).encode("utf-8"),
            file_name=f"{safe_download_stem(config.get('scenarioLabel'), 'strategy_comparison')}.csv",
            mime="text/csv",
        )
    with st.expander("Technical reproducibility information", expanded=False):
        st.dataframe(arrow_safe_dataframe(comparison.get("scenarioSummaries", [])), use_container_width=True, hide_index=True)
        if comparison.get("pairedComparisons"):
            st.dataframe(arrow_safe_dataframe(comparison["pairedComparisons"]), use_container_width=True, hide_index=True)
        if comparison.get("pairedReplicateComparisons"):
            st.dataframe(arrow_safe_dataframe(comparison["pairedReplicateComparisons"]), use_container_width=True, hide_index=True)
            st.dataframe(arrow_safe_dataframe(comparison.get("pairedDifferenceSummaries", [])), use_container_width=True, hide_index=True)
        if comparison.get("commonSeedNonpairedDiagnostics"):
            st.dataframe(arrow_safe_dataframe(comparison["commonSeedNonpairedDiagnostics"]), use_container_width=True, hide_index=True)
    st.caption("Preview outcomes may be fractional. Simulation percentiles describe finite-population variation, not confidence intervals.")

with st.expander("Sensitivity and early-review tools", expanded=False):
    tab_sensitivity, tab_early = st.tabs(["Explore sensitivity", "Review early screening results"])
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
        planned = st.number_input(
            "Planned total screened",
            min_value=0,
            value=int(round(float(config.get("screenCoverage", 0.3)) * int(config.get("N", 0)))),
            step=1,
        )
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
                st.caption("Likelihood uses the aggregate-binomial approximation; detailed calibration metadata is retained in exports.")
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
