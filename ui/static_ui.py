import streamlit as st
from engine.static.exec_static import run_static_model
from engine.static.static_params import load_static_parameters

def render_static_ui():
    st.header("📊 Static TB Model")

    params = load_static_parameters()

    if st.button("Run static model"):
        try:
            results = run_static_model(params)
            st.success("Static model complete!")
            st.write(results)
        except Exception as e:
            st.error(f"Static model failed: {e}")
