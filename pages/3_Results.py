from __future__ import annotations

from pathlib import Path

import streamlit as st

from app.display import (
    arrow_safe_dataframe,
    economics_assumptions_json,
    economics_summary_csv,
    safe_download_stem,
)
from app.icon_arrays import build_100_person_visual_data, render_100_person_summary
from app.results_workbook import build_results_workbook
from app.state import init_session_state
from engine.apy.scenario import DIRECT_EFFECTS_SCOPE_STATEMENT


init_session_state()

st.title("Results")
st.caption("Review analysis outputs, assumptions, downloads and readiness warnings.")

bundle = st.session_state.get("results_bundle")
if not bundle:
    st.info("Run the analysis to create results.")
    st.stop()

metadata = bundle.get("metadata", {})
headline = bundle.get("headline", {})
technical = bundle.get("technical", {})
downloads = bundle.get("downloads", {})
economics = st.session_state.get("economics_results")
economics_config = st.session_state.get("economics_config")
scenario_label = metadata.get("scenarioLabel")

if st.session_state.get("results_stale"):
    st.warning("These results are stale because analysis inputs changed after the last run.")

scope_statement = (
    technical.get("interfaceConfig", {})
    .get("scenario", {})
    .get("scopeStatement", DIRECT_EFFECTS_SCOPE_STATEMENT)
)
st.info(scope_statement)

with st.expander("Technical information", expanded=False):
    metadata_rows = [
        {"field": key, "value": value}
        for key, value in metadata.items()
    ]
    st.dataframe(arrow_safe_dataframe(metadata_rows), width="content", hide_index=True)

st.subheader("Headline")
if headline.get("keyMetricsRows"):
    st.markdown("Key metrics")
    st.dataframe(
        arrow_safe_dataframe(headline["keyMetricsRows"]),
        width="stretch",
    )
if headline.get("summaryRows"):
    st.markdown("Summary")
    st.dataframe(
        arrow_safe_dataframe(headline["summaryRows"]),
        width="stretch",
    )
if not headline.get("keyMetricsRows") and not headline.get("summaryRows"):
    st.json(headline, expanded=False)

visual_rows = build_100_person_visual_data(technical.get("eventLedger"))
if visual_rows:
    render_100_person_summary(
        visual_rows,
        title="What this means per 100 eligible people",
    )

st.subheader("Economics")
if not economics:
    st.info("Economics has not been run for these results yet.")
    if economics_config:
        st.markdown("Downloads")
        st.download_button(
            "Download economics assumptions JSON",
            data=economics_assumptions_json(economics_config),
            file_name=f"{safe_download_stem(scenario_label, 'economics_assumptions')}.json",
            mime="application/json",
        )
    st.page_link("pages/4_Economics.py", label="Open Evidence & Assumptions")
else:
    if st.session_state.get("dirty_economics") or st.session_state.get("results_stale"):
        st.warning("Economic results are stale. Open Evidence & Assumptions and rerun the analysis.")
    else:
        st.success("Economics results are available.")

    summary_rows = economics.get("summaryRows") or []
    if summary_rows:
        st.markdown("Summary")
        st.dataframe(arrow_safe_dataframe(summary_rows), width="stretch")
    else:
        st.info("No economics summary rows were returned.")

    st.markdown("Downloads")
    if summary_rows:
        st.download_button(
            "Download economics summary CSV",
            data=economics_summary_csv(economics),
            file_name=f"{safe_download_stem(scenario_label, 'economics_summary')}.csv",
            mime="text/csv",
        )

    if economics_config:
        st.download_button(
            "Download economics assumptions JSON",
            data=economics_assumptions_json(economics_config),
            file_name=f"{safe_download_stem(scenario_label, 'economics_assumptions')}.json",
            mime="application/json",
        )

    st.markdown("Status")
    status = economics.get("status", {})
    status_rows = [
        {"field": "last_economics_run_at", "value": st.session_state.get("last_economics_run_at")},
        {"field": "isComplete", "value": status.get("isComplete")},
        {"field": "missingInputs", "value": status.get("missingInputs")},
        {"field": "notCalculated", "value": status.get("notCalculated")},
    ]
    st.dataframe(arrow_safe_dataframe(status_rows), width="stretch", hide_index=True)
    st.page_link("pages/4_Economics.py", label="Open Evidence & Assumptions")

with st.expander("Additional technical information", expanded=False):
    technical_summary = {
        "available": technical.get("available"),
        "source": technical.get("source"),
        "exampleCohortMeta": technical.get("exampleCohortMeta"),
        "rawMeta": technical.get("rawMeta"),
    }
    st.dataframe(
        arrow_safe_dataframe(
            [{"field": key, "value": value}
            for key, value in technical_summary.items()
            ]
        ),
        width="stretch",
        hide_index=True,
    )
    st.markdown("Calibration")
    st.json(technical.get("calibration", {}), expanded=False)
    st.markdown("Interface config")
    st.json(technical.get("interfaceConfig", {}), expanded=False)
    st.markdown("Dynamic comparison")
    dynamic_comparison = technical.get("dynamicComparison", {})
    if isinstance(dynamic_comparison, dict):
        metric_rows = dynamic_comparison.get("metricRows") or []
        overview = [
            {"field": key, "value": value}
            for key, value in dynamic_comparison.items()
            if key != "metricRows"
        ]
        st.dataframe(arrow_safe_dataframe(overview), width="stretch", hide_index=True)
        if metric_rows:
            st.dataframe(arrow_safe_dataframe(metric_rows), width="stretch", hide_index=True)
    else:
        st.json(dynamic_comparison, expanded=False)

st.subheader("Downloads")
if st.session_state.get("results_stale"):
    st.warning("Excel workbook download is disabled until the analysis is rerun with the current inputs.")
else:
    workbook_bytes = build_results_workbook(
        config=technical.get("interfaceConfig", {}),
        bundle=bundle,
        backend_status=st.session_state.get("backend_status"),
        economics_results=economics,
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
    download_rows = [
        {"field": key, "value": value}
        for key, value in downloads.items()
        if key not in {"available", "payload"}
    ]
    if download_rows:
        st.dataframe(arrow_safe_dataframe(download_rows), width="stretch", hide_index=True)

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
else:
    st.info("No downloads are available for the current results.")
