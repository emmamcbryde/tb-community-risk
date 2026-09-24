from __future__ import annotations

import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.state import (
    RESULTS_CONFIG_KEY,
    RESULTS_KEY,
    STALE_KEY,
    bundle_outcome_records,
    init_general_state,
    page_link,
)
from app.icon_arrays import build_100_person_visual_data, render_100_person_summary
from app.results_page_display import (
    detailed_rows_for_display,
    format_interval_cells_for_display,
    key_metric_rows_for_display,
)


init_general_state()
st.title("Results")

bundle = st.session_state.get(RESULTS_KEY)
if not bundle:
    st.info("Run the analysis to create results.")
    page_link("general_pages/3_Run_analysis.py", label="Open Run analysis")
    st.stop()

config = st.session_state.get(RESULTS_CONFIG_KEY) or {}
metadata = bundle.get("metadata", {})
headline = bundle.get("headline", {})
technical = bundle.get("technical", {})
model_type = metadata.get("modelType") or metadata.get("analysisMethod")
link = config.get("generalProfileLink") or {}

if st.session_state.get(STALE_KEY):
    st.warning("These results are out of date because inputs changed after the last run.")
else:
    st.success("Results are current for the saved inputs.")

if model_type == "agent_based":
    st.caption(f"Stochastic analysis: {metadata.get('nReps')} simulated populations, seed {metadata.get('seed')}.")
else:
    st.caption("Deterministic expected-value preview: a single calculation without simulation variation.")
if link:
    snapshot = link.get("incidenceSnapshotId") or "none"
    st.caption(f"Profile {link.get('profileId')} · incidence snapshot {snapshot}.")
st.caption(
    "Results estimate direct benefits, harms and costs for people screened and treated. "
    "Transmission-mediated benefits are not yet included."
)

dynamic_rows = ((technical.get("dynamicComparison") or {}).get("metricRows")) or []
key_rows = key_metric_rows_for_display(headline.get("keyMetricsRows"), dynamic_rows, model_type=model_type)
detail_rows = detailed_rows_for_display(headline.get("summaryRows"), headline.get("keyMetricsRows"), model_type=model_type)

st.subheader("Key results")
if key_rows:
    st.dataframe(arrow_safe_dataframe(format_interval_cells_for_display(key_rows)), use_container_width=True, hide_index=True)
    if model_type == "agent_based":
        st.caption(
            "Median and 95% range across simulated populations. These ranges describe simulation variation; "
            "they are not confidence intervals."
        )
else:
    st.info("Key results are unavailable for this run.")

with st.expander("Detailed summary"):
    if detail_rows:
        st.dataframe(arrow_safe_dataframe(format_interval_cells_for_display(detail_rows)), use_container_width=True, hide_index=True)
    else:
        st.write("No detailed summary is available.")

visual_rows = build_100_person_visual_data(bundle_outcome_records(bundle))
if visual_rows:
    render_100_person_summary(visual_rows, title="What this means per 100 eligible people")

st.warning(
    "Results rely on demonstration assumptions that are not specific to any country. "
    "Review Inputs still needing local evidence before using results for planning."
)
page_link("general_pages/5_Health_economics.py", label="Continue to Health economics")
