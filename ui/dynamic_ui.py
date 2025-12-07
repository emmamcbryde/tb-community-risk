import streamlit as st
import pandas as pd
import numpy as np

from engine.dynamic.exec_dynamic import run_dynamic_model
from engine.dynamic.dynamic_params import load_dynamic_parameters
from engine.infection_backcast import calc_ari_from_incidence, infection_prob_by_age_split


# --- Default fallback age distribution ---
def default_age_distribution():
    ages = list(range(0, 101))
    prop = np.array([1/101] * 101)
    return pd.DataFrame({"AgeGroup": ages, "Proportion": prop})


# --- Country-based population loader (shared logic with static model) ---
def load_population_data(country_code="AUS", file_path="data/population_age_latest.csv"):

    df = pd.read_csv(file_path)
    df_country = df[df["iso_code"].str.upper() == country_code.upper()]

    if df_country.empty:
        st.warning(f"No age distribution found for {country_code}. Using default global.")
        return default_age_distribution()

    # Convert to proportions
    total_pop = df_country["population"].sum()
    df_country["Proportion"] = df_country["population"] / total_pop

    # Create 5-year bins
    bin_edges = list(range(0,105,5))
    bin_labels = [f"{bin_edges[i]}–{bin_edges[i+1]-1}" for i in range(len(bin_edges)-1)]
    bin_labels.append("100+")

    df_country["AgeBin"] = pd.cut(
        df_country["age"],
        bins=bin_edges+[200],
        labels=bin_labels,
        right=False
    )

    df_age_groups = df_country.groupby("AgeBin", as_index=False)[["population","Proportion"]].sum()

    return df_age_groups



def render_dynamic_ui():

    st.header("📈 Dynamic LTBI → TB Model")

    st.sidebar.header("Dynamic Model Inputs")

    # --- Community parameters ---
    population = st.sidebar.number_input("Population size", min_value=100, value=10000)
    user_incidence = st.sidebar.number_input("Baseline TB incidence (per 100k/yr)", 0, 500, 30)
    time_horizon = st.sidebar.slider("Time horizon (years)", 1, 30, 20)

    # --- Infection parameters ---
    beta = st.sidebar.number_input(
        "Transmission rate β (infections per active TB case per year)",
        min_value=0.0, max_value=50.0, value=8.0, step=0.1
    )

    # --- Risk factors ---
    smoker_pct = st.sidebar.slider("Smoker population (%)", 0, 100, 30)
    diabetes_pct = st.sidebar.slider("Diabetes (%)", 0, 100, 10)
    renal_pct = st.sidebar.slider("Renal impairment (%)", 0, 100, 5)
    immune_pct = st.sidebar.slider("Immunosuppressed (%)", 0, 100, 3)
    # --- Testing and Treatment Inputs ---
    testing_method = st.sidebar.selectbox("Testing method", ["None", "TST", "IGRA"])
    treatment_method = st.sidebar.selectbox("Treatment regimen", ["None", "1HP", "3HP", "4R", "6H", "9H"])

    coverage_testing = st.sidebar.slider("Testing coverage", 0.0, 1.0, 0.5)
    coverage_treatment = st.sidebar.slider("Treatment coverage", 0.0, 1.0, 0.7)

    rollout_years = st.sidebar.slider("Rollout duration (years)", 1, 5, 3)

    pre_det_months = st.sidebar.number_input("Mean time to TB diagnosis (pre-intervention months)", 1.0, 60.0, 12.0, 0.5)
    post_det_months = st.sidebar.number_input("Mean time to TB diagnosis (post-intervention months)", 1.0, 60.0, 6.0, 0.5)

    delta_pre = 12.0 / pre_det_months
    delta_post = 12.0 / post_det_months

    # --- Age distribution ---
    st.sidebar.subheader("Age distribution")
    method = st.sidebar.radio(
        "Choose age distribution source:",
        ["Country ISO code (recommended)", "Upload custom CSV", "Default global"]
    )

    if method == "Country ISO code (recommended)":
        country = st.sidebar.text_input("Enter ISO3 code", "AUS")
        age_df = load_population_data(country)

    elif method == "Upload custom CSV":
        file = st.sidebar.file_uploader("Upload age distribution CSV", type="csv")
        if file:
            df = pd.read_csv(file)
            if "AgeGroup" in df.columns and "Proportion" in df.columns:
                age_df = df
            else:
                st.error("CSV must contain AgeGroup and Proportion columns.")
                age_df = default_age_distribution()
        else:
            age_df = default_age_distribution()

    else:
        age_df = default_age_distribution()

    # --- Show age structure ---
    st.subheader("Age Distribution in 5-year bins")
    st.dataframe(age_df)

    # --- Convert to 1-year ages if needed for LTBI backcast ---
    # Use the proportions to backcast LTBI prevalence
    # (Assuming equal distribution within bins)
    ages = list(range(0,101))
    age_counts = {a: population/101 for a in ages}  # Replace with fine-grain expansion later if needed

    # --- LTBI Back-calculation ---
    inc_hist = { -k: user_incidence for k in ages }
    ari_hist = calc_ari_from_incidence(inc_hist)
    ltbi_ever, ltbi_recent, ltbi_remote = infection_prob_by_age_split(ages, ari_hist)

    # --- Run simulation button ---
    if st.sidebar.button("Run Dynamic Simulation"):
        st.info("Running dynamic model...")

        params = load_dynamic_parameters()

        # Add dynamic inputs to params dict
        params["beta"] = beta
        params["smoker_pct"] = smoker_pct
        params["diabetes_pct"] = diabetes_pct
        params["renal_pct"] = renal_pct
        params["immune_pct"] = immune_pct
        params["ltbi_ever"] = ltbi_ever
        params["ltbi_recent"] = ltbi_recent
        params["age_counts"] = age_counts

        try:
            results = run_dynamic_model(
                params=params,
                years=time_horizon
            )

            st.success("Simulation complete!")
            st.subheader("Incidence Over Time")
            st.line_chart(results["incidence"])

        except Exception as e:
            st.error(f"Dynamic model failed: {e}")
