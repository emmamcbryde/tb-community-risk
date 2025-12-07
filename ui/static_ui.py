import streamlit as st
import numpy as np
import pandas as pd
from engine.static.exec_static import run_static_model
from engine.static.static_params import load_static_parameters

# --- Default global age distribution (fallback) ---
def default_age_distribution():
    ages = list(range(0, 101))
    prop = np.array([1/101] * 101)
    return pd.DataFrame({"AgeGroup": ages, "Proportion": prop})

# --- Country-based population loader ---
def load_population_data(country_code="AUS", file_path="data/population_age_latest.csv"):
    df = pd.read_csv(file_path)

    df_country = df[df["iso_code"].str.upper() == country_code.upper()]
    if df_country.empty:
        st.warning(f"No data found for {country_code}. Using default global.")
        return default_age_distribution()

    # Convert 1-year ages → proportions
    total_pop = df_country["population"].sum()
    df_country["Proportion"] = df_country["population"] / total_pop

    # Create 5-year bins
    bin_edges = list(range(0,105,5))
    bin_labels = [f"{bin_edges[i]}–{bin_edges[i+1]-1}" for i in range(len(bin_edges)-1)]
    bin_labels.append("100+")

    df_country["AgeBin"] = pd.cut(
        df_country["age"], 
        bins=bin_edges + [200],
        labels=bin_labels,
        right=False
    )

    # Group by bins → FIXED HERE
    df_age_groups = df_country.groupby("AgeBin", as_index=False)[["population", "Proportion"]].sum()

    return df_age_groups


# --- MAIN STATIC UI (only definition) ---
def render_static_ui():

    st.header("📊 Static TB Incidence Model")

    st.sidebar.header("Static Model Inputs")

    # --- Community parameters ---
    population = st.sidebar.number_input("Population size", min_value=100, value=10000)
    user_incidence = st.sidebar.number_input("Baseline TB incidence (per 100k/yr)", 0, 500, 30)
    time_horizon = st.sidebar.slider("Time horizon (years)", 1, 30, 20)

    # --- Risk factors ---
    smoker_pct = st.sidebar.slider("Smoker population (%)", 0, 100, 30)
    diabetes_pct = st.sidebar.slider("Diabetes (%)", 0, 100, 10)
    renal_pct = st.sidebar.slider("Renal impairment (%)", 0, 100, 5)
    immune_pct = st.sidebar.slider("Immunosuppressed (%)", 0, 100, 3)

    # --- Age distribution selection ---
    st.sidebar.subheader("Age distribution")
    country = st.sidebar.text_input("Enter ISO3 code", "AUS")

    age_df = load_population_data(country)

    st.subheader("Age distribution in 5-year bins")
    st.dataframe(age_df)

    # --- Run model ---
    if st.sidebar.button("Run Static Simulation"):
        try:
            params = load_static_parameters()
            results = run_static_model(params)
            st.success("Static model complete!")
            st.write(results)
        except Exception as e:
            st.error(f"Static model failed: {e}")
