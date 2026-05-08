from __future__ import annotations

from pathlib import Path

import streamlit as st

from app.display import arrow_safe_dataframe
from app.state import init_session_state


init_session_state()

st.title("Results")
st.caption("Displays the portable app-facing results bundle.")

bundle = st.session_state.get("results_bundle")
if not bundle:
    st.info("Run the model to create a results bundle.")
    st.stop()

metadata = bundle.get("metadata", {})
headline = bundle.get("headline", {})
technical = bundle.get("technical", {})
downloads = bundle.get("downloads", {})

if st.session_state.get("results_stale"):
    st.warning("These results are stale because scenario inputs changed after the last run.")

st.subheader("Metadata")
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

st.subheader("Technical Details")
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

with st.expander("Calibration"):
    st.json(technical.get("calibration", {}), expanded=False)

with st.expander("Interface config"):
    st.json(technical.get("interfaceConfig", {}), expanded=False)

st.subheader("Downloads")
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
    st.info("No downloads are available in the current results bundle.")
