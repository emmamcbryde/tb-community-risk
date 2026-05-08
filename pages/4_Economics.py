from __future__ import annotations

import streamlit as st

from app.state import get_backend, init_session_state


init_session_state()
backend = get_backend()

st.title("Economics")
st.caption("Initial economics state shell backed by MATLAB.")

if st.button("Load economics defaults", type="primary"):
    try:
        st.session_state["economics_config"] = backend.default_economics_config()
        st.session_state["economics_results"] = None
        st.session_state["dirty_economics"] = False
        st.success("Default economics config loaded.")
    except Exception as exc:
        st.error(f"Could not load economics defaults: {exc}")

econ_config = st.session_state.get("economics_config")
if econ_config:
    st.subheader("Current Economics Config")
    st.json(econ_config, expanded=False)
else:
    st.info("Load MATLAB economics defaults to initialize economics state.")

st.warning(
    "Running economics requires a backend results object. The first scaffold does "
    "not persist raw MATLAB-native results in session state."
)
