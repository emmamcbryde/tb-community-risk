from __future__ import annotations

import streamlit as st

from app.display import arrow_safe_dataframe
from app.state import (
    get_backend,
    get_backend_name,
    init_session_state,
    matlab_backend_enabled,
    record_message,
    set_backend_name,
    sync_backend_status,
)


init_session_state()

st.title("Technical settings")
st.caption("Research and Development controls for local validation and implementation diagnostics.")

with st.expander("Implementation selection", expanded=False):
    labels = {
        "python_apy": "Python APY v9 port",
        "matlab": "MATLAB v9 reference",
    }
    options = ["python_apy"]
    if matlab_backend_enabled():
        options.append("matlab")
    current = get_backend_name()
    selected = st.selectbox(
        "APY backend",
        options,
        index=options.index(current) if current in options else 0,
        format_func=lambda value: labels.get(value, value),
    )
    if selected != current:
        try:
            set_backend_name(selected)
            st.rerun()
        except Exception as exc:
            record_message("error", f"Could not change implementation: {exc}")
            st.error(str(exc))
    if not matlab_backend_enabled():
        st.info("MATLAB reference backend is unavailable unless explicitly enabled in a local validation environment.")

with st.expander("Backend diagnostics", expanded=False):
    backend = get_backend()
    status = backend.status()
    sync_backend_status(status)
    st.dataframe(
        arrow_safe_dataframe([{"field": key, "value": value} for key, value in status.items()]),
        hide_index=True,
        width="stretch",
    )

with st.expander("Technical metadata", expanded=False):
    bundle = st.session_state.get("results_bundle") or {}
    metadata = bundle.get("metadata", {}) if isinstance(bundle, dict) else {}
    st.dataframe(
        arrow_safe_dataframe([{"field": key, "value": value} for key, value in metadata.items()]),
        hide_index=True,
        width="stretch",
    )
