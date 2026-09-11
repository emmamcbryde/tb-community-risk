from __future__ import annotations

from pathlib import Path

import streamlit as st

from app.display import (
    arrow_safe_dataframe,
    safe_download_stem,
)
from app.icon_arrays import build_100_person_visual_data, render_100_person_summary
from app.results_page_display import (
    detailed_rows_for_display,
    key_metric_rows_for_display,
)
from app.results_workbook import build_results_workbook
from app.state import init_session_state
from engine.apy.scenario import DIRECT_EFFECTS_SCOPE_STATEMENT


init_session_state()

st.title("Results")

bundle = st.session_state.get("results_bundle")
if not bundle:
    st.info("Run the analysis to create results.")
    st.stop()

metadata = bundle.get("metadata", {})
headline = bundle.get("headline", {})
technical = bundle.get("technical", {})
downloads = bundle.get("downloads", {})
economics_config = st.session_state.get("economics_config")
scenario_label = metadata.get("scenarioLabel")

if st.session_state.get("results_stale"):
    st.warning("These results are stale because analysis inputs changed after the last run.")
else:
    st.success("Results are current for the saved analysis inputs.")

if scenario_label:
    st.markdown(f"**Scenario:** {scenario_label}")

scope_statement = (
    technical.get("interfaceConfig", {})
    .get("scenario", {})
    .get("scopeStatement", DIRECT_EFFECTS_SCOPE_STATEMENT)
)
st.caption(scope_statement)

dynamic_metric_rows = []
dynamic_comparison = technical.get("dynamicComparison", {})
if isinstance(dynamic_comparison, dict):
    dynamic_metric_rows = dynamic_comparison.get("metricRows") or []

key_rows = key_metric_rows_for_display(
    headline.get("keyMetricsRows"),
    dynamic_metric_rows,
)
detail_rows = detailed_rows_for_display(
    headline.get("summaryRows"),
    headline.get("keyMetricsRows"),
)

st.subheader("Key metrics")
if key_rows:
    st.dataframe(
        arrow_safe_dataframe(key_rows),
        use_container_width=True,
        hide_index=True,
    )
else:
    st.info("Key metrics are unavailable for these results.")

st.caption(
    "Median, low 95% and high 95% summarise the distribution across repeated "
    "simulated populations. They are not confidence intervals."
)

st.subheader("Detailed summary table")
if detail_rows:
    st.dataframe(
        arrow_safe_dataframe(detail_rows),
        use_container_width=True,
        hide_index=True,
    )
else:
    st.json(headline, expanded=False)

visual_rows = build_100_person_visual_data(technical.get("eventLedger"))
if visual_rows:
    st.caption(
        "Values may include decimals because they are averages across repeated "
        "simulated populations. Unless otherwise stated, these values are per "
        "100 eligible people."
    )
    render_100_person_summary(
        visual_rows,
        title="What this means per 100 eligible people",
    )

with st.expander("Export results", expanded=False):
    if st.session_state.get("results_stale"):
        st.warning("Workbook download is disabled until the analysis is rerun with the current inputs.")
    else:
        workbook_bytes = build_results_workbook(
            config=technical.get("interfaceConfig", {}),
            bundle=bundle,
            backend_status=st.session_state.get("backend_status"),
            economics_results=st.session_state.get("economics_results"),
            economics_config=economics_config,
            results_stale=False,
            dirty_economics=bool(st.session_state.get("dirty_economics")),
            decision_analysis_results={
                "scenarioComparison": st.session_state.get("decision_scenario_comparison"),
                "sensitivity": st.session_state.get("decision_sensitivity"),
                "threshold": st.session_state.get("decision_threshold"),
                "earlyReview": st.session_state.get("decision_early_review"),
            },
        )
        st.download_button(
            "Download consolidated results workbook",
            data=workbook_bytes,
            file_name=f"{safe_download_stem(scenario_label, 'APY_results')}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )

    if downloads.get("available"):
        for label, key in (("Summary CSV", "summaryCsv"), ("Key metrics CSV", "keyMetricsCsv")):
            path_value = downloads.get(key)
            if not path_value:
                continue
            path = Path(str(path_value))
            if path.is_file():
                st.download_button(
                    label,
                    data=path.read_bytes(),
                    file_name=path.name,
                    mime="text/csv",
                )

st.warning(
    "Outputs remain provisional where evidence inputs are unresolved. Review "
    "Evidence & Assumptions before using results as final policy evidence."
)
st.page_link("pages/4_Economics.py", label="Continue to Health Economics")
