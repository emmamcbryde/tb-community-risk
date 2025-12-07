import streamlit as st
from engine.dynamic.exec_dynamic import run_dynamic_model
from engine.dynamic.dynamic_params import load_dynamic_parameters

def render_dynamic_ui():
    st.header("📈 Dynamic TB Model")

    # Load parameters
    params = load_dynamic_parameters()

    # User options
    years = st.slider("Simulation years:", 1, 50, 20)

    if st.button("Run dynamic simulation"):
        try:
            results = run_dynamic_model(params, years)
            st.success("Simulation complete!")
            st.line_chart(results["incidence"])
        except Exception as e:
            st.error(f"Simulation failed: {e}")
