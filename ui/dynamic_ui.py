import streamlit as st
import pandas as pd
import numpy as np

from engine.dynamic.exec_dynamic import run_dynamic_model
from engine.dynamic.dynamic_params import load_dynamic_parameters
from engine.infection_backcast import calc_ari_from_incidence, infection_prob_by_age_split

# =====================================================
# Default fallback: equal distribution 0–100 years
# =====================================================
def default_age_distribution():
    ages = list(range(0, 101))
    prop = np.array([1 / 101] * 101)
    return pd.DataFrame({"AgeGroup": ages, "Proportion": prop})


# =====================================================
# Country-level loader: population_age_latest.csv
# =====================================================
def load_population_data(country_code="AUS", file_path="data/population_age_latest.csv"):
    df = pd.read_csv(file_path)

    df_country = df[df["iso_code"].str.upper() == country_code.upper()]
    if df_country.empty:
        st.warning(f"No age distribution found for {country_code}. Using default.")
        return default_age_distribution()

    # Convert to proportions (for plotting)
    total_pop = df_country["population"].sum()
    df_country["Proportion"] = df_country["population"] / total_pop

    # Produce 5-year bins for display (optional)
    bin_edges = list(range(0, 105, 5))
    bin_labels = [f"{bin_edges[i]}–{bin_edges[i+1]-1}" for i in range(len(bin_edges)-1)]
    bin_labels.append("100+")

    df_country["AgeBin"] = pd.cut(
        df_country["age"],
        bins=bin_edges + [200],
        labels=bin_labels,
        right=False
    )

    df_age_groups = df_country.groupby("AgeBin", as_index=False)[["population", "Proportion"]].sum()

    return df_age_groups, df_country  # second df has age = 0–100 rows


# =====================================================
# Dynamic Model UI
# =====================================================
def render_dynamic_ui():

    st.header("📈 Dynamic LTBI → TB Model")

    st.sidebar.header("Dynamic Model Inputs")

    # -----------------------------
    # Community & epidemiology
    # -----------------------------
    population = st.sidebar.number_input("Population size", min_value=100, value=10000)
    user_incidence = st.sidebar.number_input("Baseline TB incidence (per 100k/yr)", 0, 500, 30)
    time_horizon = st.sidebar.slider("Time horizon (years)", 1, 30, 20)

    # -----------------------------
    # Transmission
    # -----------------------------
    beta = st.sidebar.number_input(
        "Transmission rate β (infections per active case per year)",
        min_value=0.0, max_value=50.0, value=8.0, step=0.1
    )

    # -----------------------------
    # Risk factors
    # -----------------------------
    smoker_pct = st.sidebar.slider("Smoker population (%)", 0, 100, 30)
    diabetes_pct = st.sidebar.slider("Diabetes (%)", 0, 100, 10)
    renal_pct = st.sidebar.slider("Renal impairment (%)", 0, 100, 5)
    immune_pct = st.sidebar.slider("Immunosuppressed (%)", 0, 100, 3)

    # -----------------------------
    # Testing & treatment options
    # -----------------------------
    testing_method = st.sidebar.selectbox("Testing method", [ "TST", "IGRA", "None"])
    treatment_method = st.sidebar.selectbox(
        "Treatment regimen", ["1HP", "3HP", "4R", "6H", "9H", "None"]
    )

    coverage_testing = st.sidebar.slider("Testing coverage", 0.0, 1.0, 0.5)
    coverage_treatment = st.sidebar.slider("Treatment coverage", 0.0, 1.0, 0.7)
    rollout_years = st.sidebar.slider("Rollout duration (years)", 1, 5, 3)

    # TB case detection (affects infectious period)
    pre_det_months = st.sidebar.number_input(
        "Mean time to TB diagnosis (pre-intervention months)",
        1.0, 60.0, 12.0, 0.5
    )
    post_det_months = st.sidebar.number_input(
        "Mean time to TB diagnosis (post-intervention months)",
        1.0, 60.0, 6.0, 0.5
    )
    delta_pre = 12.0 / pre_det_months
    delta_post = 12.0 / post_det_months

    # -----------------------------
    # Age distribution selection
    # -----------------------------
    st.sidebar.subheader("Age distribution")
    method = st.sidebar.radio(
        "Choose age distribution source:",
        ["Country ISO code (recommended)", "Upload custom CSV", "Default global"]
    )

    # 1. Country ISO code
    if method == "Country ISO code (recommended)":
        country = st.sidebar.text_input("Enter ISO3 code", "AUS")
        age_df_display, df_country = load_population_data(country)

    # 2. Custom CSV upload
    elif method == "Upload custom CSV":
        file = st.sidebar.file_uploader("Upload age distribution CSV", type="csv")
        if file:
            df = pd.read_csv(file)
            if "AgeGroup" in df.columns and "Proportion" in df.columns:
                age_df_display = df
                # Construct 1-year counts based on proportions
                df_country = pd.DataFrame({
                    "age": df["AgeGroup"],
                    "population": (df["Proportion"] * population).values
                })
            else:
                st.error("CSV must contain AgeGroup and Proportion columns.")
                age_df_display = default_age_distribution()
                df_country = pd.DataFrame({
                    "age": range(0, 101),
                    "population": [population/101] * 101
                })
        else:
            age_df_display = default_age_distribution()
            df_country = pd.DataFrame({
                "age": range(0, 101),
                "population": [population/101] * 101
            })

    # 3. Default global
    else:
        age_df_display = default_age_distribution()
        df_country = pd.DataFrame({
            "age": range(0, 101),
            "population": [population/101] * 101
        })

    # Show age distribution
    st.subheader("Age Distribution in 5-year bins")
    st.dataframe(age_df_display)

    # Build age_counts
    # National population counts → rescaled to chosen population
    total_pop_country = df_country["population"].sum()
    age_counts = {
        int(row["age"]): population * (row["population"] / total_pop_country)
        for _, row in df_country.iterrows()
    }

    ages = sorted(age_counts.keys())

    # -----------------------------
    # LTBI backcast
    # -----------------------------
    inc_hist = {-k: user_incidence for k in ages}
    ari_hist = calc_ari_from_incidence(inc_hist)
    ltbi_ever, ltbi_recent, ltbi_remote = infection_prob_by_age_split(ages, ari_hist)

    # -----------------------------
    # RUN MODEL
    # -----------------------------
    if st.sidebar.button("Run Dynamic Simulation"):
        st.info("Running dynamic model...")

        params_base = load_dynamic_parameters()
        params_int = params_base.copy()

        # Common parameters to both baseline and intervention
        for p in [params_base, params_int]:
            p["beta"] = beta
            p["smoker_pct"] = smoker_pct
            p["diabetes_pct"] = diabetes_pct
            p["renal_pct"] = renal_pct
            p["immune_pct"] = immune_pct
            p["ltbi_ever"] = ltbi_ever
            p["ltbi_recent"] = ltbi_recent
            p["age_counts"] = age_counts
            p["delta_pre"] = delta_pre
            p["delta_post"] = delta_post

        # ----------------------------------
        # Baseline scenario: NO intervention
        # ----------------------------------
        params_base["coverage_testing"] = 0.0
        params_base["coverage_treatment"] = 0.0
        params_base["treatment_method"] = "None"
        params_base["testing_method"] = "None"
        params_base["rollout_years"] = 0
        params_base["initial_incidence_per_100k"] = user_incidence

        # ----------------------------------
        # Intervention scenario (user inputs)
        # ----------------------------------
        params_int["coverage_testing"] = coverage_testing
        params_int["coverage_treatment"] = coverage_treatment
        params_int["treatment_method"] = treatment_method
        params_int["testing_method"] = testing_method
        params_int["rollout_years"] = rollout_years
        params_int["initial_incidence_per_100k"] = user_incidence
        try:
            baseline = run_dynamic_model(params_base, years=time_horizon, intervention=False)
            intervention = run_dynamic_model(params_int, years=time_horizon, intervention=True)

            # Construct output dataframe
            df_out = pd.DataFrame({
                "Year": baseline["time"],
                "Baseline_incidence": baseline["incidence"] * 100000 / sum(age_counts.values()),
                "Intervention_incidence": intervention["incidence"] * 100000 / sum(age_counts.values())
            })
            df_out["Cases_averted"] = df_out["Baseline_incidence"] - df_out["Intervention_incidence"]

            st.success("Simulation complete!")
            st.subheader("📈 Annual Incidence Over Time (per 100,000)")
            st.line_chart(df_out[["Baseline_incidence", "Intervention_incidence"]])

            st.subheader("🔍 Cases Averted")
            st.write(df_out)

        except Exception as e:
            st.error(f"Dynamic model failed: {e}")
