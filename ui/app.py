import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(page_title="TB Community Risk Model", layout="wide")

st.title("🧮 Community-level TB Risk Model")

st.sidebar.header("Community Parameters")

# Example inputs
population = st.sidebar.number_input("Population size", min_value=100, value=10000, step=100)
indigenous_pct = st.sidebar.slider("Indigenous population (%)", 0, 100, 60)
smoker_pct = st.sidebar.slider("Smoker population (%)", 0, 100, 30)
recent_infection = st.sidebar.slider("Recent infection (%)", 0, 50, 5)
time_horizon = st.sidebar.slider("Time horizon (years)", 1, 30, 20)

st.sidebar.header("Intervention Strategy")
testing = st.sidebar.selectbox("Testing method", ["None", "TST", "IGRA"])
treatment = st.sidebar.selectbox("Treatment regimen", ["None", "3HP", "4R", "6H", "9H"])

# Simple placeholder model
base_incidence = 200  # per 100k, baseline
risk_factor = (1 + indigenous_pct/100 * 0.5 + smoker_pct/100 * 0.2 + recent_infection/50)
intervention_factor = 1.0

if testing != "None" and treatment != "None":
    intervention_factor = 0.7  # assume 30% reduction in incidence

annual_incidence = base_incidence * risk_factor * intervention_factor

# Display results
st.subheader("Results")
st.metric("Estimated annual incidence (per 100,000)", f"{annual_incidence:.1f}")

# Simulated trajectory
years = np.arange(0, time_horizon + 1)
incidence_over_time = annual_incidence * np.exp(-0.05 * years)
df = pd.DataFrame({"Year": years, "Incidence": incidence_over_time})

st.line_chart(df.set_index("Year"))

st.caption("Prototype demo — calculations will later be replaced by the full Markov engine.")
