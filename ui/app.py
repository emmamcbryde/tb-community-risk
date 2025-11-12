# ensure we can import from parent folder (project root)
import os, sys
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import streamlit as st
import pandas as pd
import numpy as np




from engine.model import simulate_community

# --- Streamlit page setup ---
st.set_page_config(page_title="TB Community Risk Model", layout="wide")

st.title("🧮 Community-level TB Risk Model")
st.markdown("Estimate tuberculosis risk in a community based on population traits and interventions.")

# --- Sidebar Inputs ---
st.sidebar.header("Community Parameters")

population = st.sidebar.number_input("Population size", min_value=100, value=10000, step=100)
indigenous_pct = st.sidebar.slider("Indigenous population (%)", 0, 100, 60)
smoker_pct = st.sidebar.slider("Smoker population (%)", 0, 100, 30)
diabetes_pct = st.sidebar.slider("Diabetes (%)", 0, 100, 10)
renal_pct = st.sidebar.slider("Renal impairment (%)", 0, 100, 5)
immune_pct = st.sidebar.slider("Immunosuppressed (%)", 0, 100, 3)
recent_infection = st.sidebar.slider("Recent infection (%)", 0, 50, 5)
time_horizon = st.sidebar.slider("Time horizon (years)", 1, 30, 20)

st.sidebar.header("Intervention Strategy")
testing = st.sidebar.selectbox("Testing method", ["None", "TST", "IGRA"])
treatment = st.sidebar.selectbox("Treatment regimen", ["None", "3HP", "4R", "6H", "9H"])
coverage = st.sidebar.slider("Intervention coverage (%)", 0, 100, 80)


st.sidebar.header("Run Simulation")
run_button = st.sidebar.button("Simulate Community")

# --- Run Model ---
if run_button:
    st.info("Running simulation...")

    inputs = {
        "population": population,
        "indigenous_pct": indigenous_pct,
        "smoker_pct": smoker_pct,
        "diabetes_pct": diabetes_pct,
        "renal_pct": renal_pct,
        "immune_pct": immune_pct,
        "recent_infection": recent_infection,
        "time_horizon": time_horizon,
        "testing": testing,
        "treatment": treatment,
        "coverage": coverage / 100,  # convert to fraction
    }

    try:
        df, summary = simulate_community(inputs, file_path="data/parameters.xlsx")

        # --- Display Results ---
        st.subheader("📊 Results Summary")
        c1, c2, c3 = st.columns(3)
        c1.metric("Peak incidence (per 100,000)", f"{summary['Peak incidence']:.1f}")
        c2.metric("Total TB cases", f"{summary['Total cases']:.0f}")
        c3.metric("Total TB deaths", f"{summary['Total TB deaths']:.0f}")

        st.subheader("📈 Annual Incidence Over Time")
        chart_df = df.set_index("Year")[["Incidence_strategy", "Incidence_baseline"]]
        st.line_chart(chart_df)

        st.caption("_Solid line = chosen strategy, dashed line = baseline (no intervention)_")

        st.caption("Model uses parameters from data/parameters.xlsx.")
        st.success("✅ Simulation complete!")

    except Exception as e:
        st.error(f"Simulation failed: {e}")

else:
    st.info("Adjust community parameters in the sidebar and click **Simulate Community** to run the model.")
