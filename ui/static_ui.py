import streamlit as st
import pandas as pd
import numpy as np
import altair as alt

from engine.dynamic.exec_dynamic import run_dynamic_model
from engine.infection_backcast import (
    calc_ari_from_incidence,
    infection_prob_by_age_split,
)


# =====================================================
# Default fallback age distribution
# =====================================================
def default_age_distribution():
    ages = list(range(0, 101))
    prop = np.array([1 / 101] * 101)
    return pd.DataFrame({"AgeGroup": ages, "Proportion": prop})


# =====================================================
# Load country-specific OWID population structure
# =====================================================
def load_population_data(
    country_code="AUS", file_path="data/population_age_latest.csv"
):

    df = pd.read_csv(file_path)
    df_country = df[df["iso_code"].str.upper() == country_code.upper()]

    if df_country.empty:
        st.warning(f"No population structure for {country_code}. Using default global.")
        return default_age_distribution(), pd.DataFrame(
            {
                "age": range(0, 101),
                "population": [1] * 101,
            }
        )

    # Convert to proportions for display
    total_pop_country = df_country["population"].sum()
    df_country["Proportion"] = df_country["population"] / total_pop_country

    # Build 5-year bins for display
    bin_edges = list(range(0, 105, 5))
    bin_labels = [
        f"{bin_edges[i]}–{bin_edges[i+1]-1}" for i in range(len(bin_edges) - 1)
    ]
    bin_labels.append("100+")

    df_country["AgeBin"] = pd.cut(
        df_country["age"], bins=bin_edges + [200], labels=bin_labels, right=False
    )

    df_age_groups = df_country.groupby("AgeBin", as_index=False)[
        ["population", "Proportion"]
    ].sum()

    return df_age_groups, df_country


# =====================================================
# Main Dynamic Model UI
# =====================================================
def render_static_ui():
    st.subheader("Static TB Model")
    st.warning(
        "The static model is a simplified approximation of the dynamic model. "
        "It does not simulate transmission feedback and should be interpreted "
        "as an approximate projection."
        "This model estimates TB incidence based on LTBI, risk factors, and interventions.\n\n"
        "• Baseline = no new LTBI test-and-treat\n"
        "• Intervention = selected LTBI testing and treatment\n"
        "• Outputs are annualised and shown as rates (per 100,000) and counts"
    )

    # --------------------------------------------------
    # Core epidemiological inputs
    # --------------------------------------------------
    population = st.sidebar.number_input("Population size", min_value=50, value=10000)
    user_incidence = st.sidebar.number_input(
        "Baseline annual incidence (per 100k)", 0, 500, 30
    )
    time_horizon = st.sidebar.slider("Time horizon (years)", 1, 30, 20)

    beta = st.sidebar.number_input(
        "Transmission rate β", min_value=0.0, max_value=50.0, value=8.0, step=0.1
    )

    # --------------------------------------------------
    # Risk factors
    # --------------------------------------------------
    smoker_pct = st.sidebar.slider(
        "Smoker (%)",
        0,
        100,
        30,
        help="Proportion of the population with current tobacco smoking.",
    )
    alcohol_pct = st.sidebar.slider(
        "Excess alcohol use (%)",
        0,
        100,
        15,
        help="Proportion of the population with excess alcohol consumption.",
    )
    diabetes_pct = st.sidebar.slider(
        "Diabetes (%)",
        0,
        100,
        10,
        help="Proportion of the population with diagnosed diabetes mellitus.",
    )
    renal_pct = st.sidebar.slider(
        "Renal impairment (%)",
        0,
        100,
        5,
        help="Proportion of the population with moderate–severe chronic kidney disease.",
    )
    HIV_treated_pct = st.sidebar.slider(
        "HIV Treated with antiretrovirals (%)",
        0,
        100,
        3,
        help="Proportion of the population with HIV under treatment.",
    )
    HIV_untreated_pct = st.sidebar.slider(
        "HIV Untreated (%)",
        0,
        100,
        3,
        help="Proportion of the population with HIV not under treatment.",
    )

    # --------------------------------------------------
    # LTBI Test & Treat (PULSE model)
    # --------------------------------------------------
    testing_method = st.sidebar.selectbox("Testing method", ["TST", "IGRA", "None"])
    treatment_method = st.sidebar.selectbox(
        "Treatment regimen", ["1HP", "3HP", "4R", "6H", "9H", "None"]
    )

    ltbi_coverage = st.sidebar.slider(
        "LTBI Test & Treat total coverage (fraction of population)", 0.0, 1.0, 0.5
    )
    rollout_years = st.sidebar.slider("Rollout duration (years)", 1, 10, 5)

    # --------------------------------------------------
    # Diagnosis delays
    # --------------------------------------------------
    pre_det_months = st.sidebar.number_input(
        "Diagnosis delay (months) before intervention", 1.0, 60.0, 12.0, 0.5
    )
    post_det_months = st.sidebar.number_input(
        "Diagnosis delay (months) after intervention", 1.0, 60.0, 6.0, 0.5
    )

    delta_pre = 12.0 / pre_det_months
    delta_post = 12.0 / post_det_months

    # --------------------------------------------------
    # Historical incidence pattern
    # --------------------------------------------------
    st.sidebar.subheader("Historical Incidence Pattern")

    hist_pattern = st.sidebar.selectbox(
        "Choose pattern:",
        [
            "Constant",
            "Falling 3%/year",
            "Rising 3%/year",
            "Upload CSV (year, incidence)",
        ],
    )

    uploaded_inc_df = None
    if hist_pattern == "Upload CSV (year, incidence)":
        inc_file = st.sidebar.file_uploader("Upload incidence CSV", type="csv")
        if inc_file:
            try:
                tmp = pd.read_csv(inc_file)
                if {"year", "incidence"}.issubset(tmp.columns):
                    uploaded_inc_df = tmp.sort_values("year")
                    st.success("Incidence history loaded.")
                else:
                    st.error("CSV must contain: year, incidence")
            except Exception as e:
                st.error(f"Could not read file: {e}")

    # --------------------------------------------------
    # AGE DISTRIBUTION INPUT
    # --------------------------------------------------
    st.sidebar.subheader("Age Distribution")

    age_method = st.sidebar.radio(
        "Choose method:",
        ["Country ISO code (recommended)", "Upload custom CSV", "Default global"],
    )

    if age_method == "Country ISO code (recommended)":
        country = st.sidebar.text_input("ISO3 code", "AUS")
        age_df_display, df_country = load_population_data(country)

    elif age_method == "Upload custom CSV":
        file = st.sidebar.file_uploader(
            "Upload CSV with AgeGroup,Proportion", type="csv"
        )
        if file:
            df = pd.read_csv(file)
            if {"AgeGroup", "Proportion"}.issubset(df.columns):
                age_df_display = df
                df_country = pd.DataFrame(
                    {
                        "age": df["AgeGroup"].astype(int),
                        "population": df["Proportion"] * population,
                    }
                )
            else:
                st.error("CSV must include AgeGroup and Proportion. Using default.")
                age_df_display = default_age_distribution()
                df_country = pd.DataFrame(
                    {"age": range(0, 101), "population": [population / 101] * 101}
                )
        else:
            age_df_display = default_age_distribution()
            df_country = pd.DataFrame(
                {"age": range(0, 101), "population": [population / 101] * 101}
            )

    else:
        age_df_display = default_age_distribution()
        df_country = pd.DataFrame(
            {"age": range(0, 101), "population": [population / 101] * 101}
        )

    st.subheader("📊 Age Distribution (5-year bins)")
    st.dataframe(age_df_display)

    # Scale OWID pop to user-selected population
    total_pop_country = df_country["population"].sum()
    age_counts = {
        int(row["age"]): population * (row["population"] / total_pop_country)
        for _, row in df_country.iterrows()
    }

    ages = sorted(age_counts.keys())

    # --------------------------------------------------
    # BUILD INCIDENCE HISTORY FOR LTBI BACK-CALC
    # --------------------------------------------------
    # --------------------------------------------------
    # BUILD INCIDENCE HISTORY FOR LTBI BACK-CALCULATION
    # --------------------------------------------------

    # ALWAYS define a default first
    inc_hist = {-k: user_incidence for k in ages}

    if hist_pattern == "Constant":
        inc_hist = {-k: user_incidence for k in ages}

    elif hist_pattern == "Falling 3%/year":
        inc_hist = {-k: user_incidence * (1.03**k) for k in ages}

    elif hist_pattern == "Rising 3%/year":
        inc_hist = {-k: user_incidence * (0.97**k) for k in ages}

    elif hist_pattern == "Upload CSV (year, incidence)":
        if uploaded_inc_df is not None:

            years = uploaded_inc_df["year"].values
            incs = uploaded_inc_df["incidence"].values

            year_min = years[0]
            year_max = years[-1]
            inc_min = np.min(incs)
            inc_max = np.max(incs)

            # Geometric trend (ignore zeros)
            ratios = []
            for i in range(1, len(incs)):
                if incs[i - 1] > 0 and incs[i] > 0:
                    ratios.append(incs[i] / incs[i - 1])
            trend = np.exp(np.mean(np.log(ratios))) if ratios else 1.0

            inc_map = dict(zip(years, incs))
            inc_hist = {}

            for a in ages:
                target_year = year_max - a

                if year_min <= target_year <= year_max:
                    nearest = min(inc_map.keys(), key=lambda y: abs(y - target_year))
                    inc_hist[-a] = inc_map[nearest]

                elif target_year > year_max:
                    k = target_year - year_max
                    extrap = incs[-1] * (trend**k)
                    inc_hist[-a] = min(extrap, inc_max)

                else:
                    k = year_min - target_year
                    extrap = incs[0] * (trend ** (-k))
                    inc_hist[-a] = max(extrap, inc_min)

        else:
            # User selected CSV mode but has not uploaded yet
            st.warning("Incidence CSV not uploaded yet – using constant incidence.")
            inc_hist = {-k: user_incidence for k in ages}

    # --------------------------------------------------
    # LTBI BACK-CALCULATION
    # --------------------------------------------------
    ari_hist = calc_ari_from_incidence(inc_hist)
    ltbi_ever, ltbi_recent, ltbi_remote = infection_prob_by_age_split(ages, ari_hist)

    # --------------------------------------------------
    # LTBI BY AGE STACKED CHART (≤ 60 YEARS)
    # --------------------------------------------------
    ltbi_age_df = pd.DataFrame(
        {
            "Age": ages,
            "LTBI_recent": 100 * pd.Series(ltbi_recent),
            "LTBI_remote": 100 * (pd.Series(ltbi_ever) - pd.Series(ltbi_recent)),
        }
    )

    ltbi_age_df = ltbi_age_df[ltbi_age_df["Age"] <= 60]
    ltbi_age_df = ltbi_age_df.melt(id_vars="Age", var_name="Type", value_name="Percent")

    st.subheader("📉 LTBI Prevalence by Age (stacked %, ages 0–60)")

    chart = (
        alt.Chart(ltbi_age_df)
        .mark_area()
        .encode(
            x="Age:Q", y="Percent:Q", color="Type:N", tooltip=["Age", "Type", "Percent"]
        )
    )
    st.altair_chart(chart, use_container_width=True)

    # --------------------------------------------------
    # RUN SIMULATIONS
    # --------------------------------------------------
    if st.sidebar.button("Run Dynamic Simulation"):

        st.info("Running baseline and intervention...")

        params_base = load_dynamic_parameters()
        params_int = params_base.copy()

        # shared parameters
        for p in (params_base, params_int):
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
            p["pre_det_months"] = pre_det_months
            p["initial_incidence_per_100k"] = user_incidence

        # baseline = no intervention
        params_base["treatment_method"] = "None"
        params_base["testing_method"] = "None"
        params_base["ltbi_coverage"] = 0.0
        params_base["rollout_years"] = 0

        # intervention
        params_int["treatment_method"] = treatment_method
        params_int["testing_method"] = testing_method
        params_int["ltbi_coverage"] = ltbi_coverage
        params_int["rollout_years"] = rollout_years

        try:
            baseline = run_dynamic_model(
                params_base, years=time_horizon, intervention=False
            )
            intervention = run_dynamic_model(
                params_int, years=time_horizon, intervention=True
            )

            total_pop = sum(age_counts.values())
            base_I = baseline["incidence"]
            int_I = intervention["incidence"]

            df_out = pd.DataFrame(
                {
                    "Year": baseline["time"],
                    "Baseline_inc_count": base_I,
                    "Intervention_inc_count": int_I,
                    "Baseline_inc_per100k": base_I * 100000 / total_pop,
                    "Intervention_inc_per100k": int_I * 100000 / total_pop,
                }
            )

            df_out["Cases_averted_count"] = (
                df_out["Baseline_inc_count"] - df_out["Intervention_inc_count"]
            )
            df_out["Cases_averted_per100k"] = (
                df_out["Baseline_inc_per100k"] - df_out["Intervention_inc_per100k"]
            )

            df_out = df_out.round(1)

            st.success("Simulation complete.")

            st.subheader("📈 Annual Incidence per 100,000")
            st.line_chart(df_out[["Baseline_inc_per100k", "Intervention_inc_per100k"]])

            st.subheader("📈 Annual Incidence (counts)")
            st.line_chart(df_out[["Baseline_inc_count", "Intervention_inc_count"]])

            st.subheader("🔍 Cases Averted")
            st.write(df_out[["Year", "Cases_averted_count", "Cases_averted_per100k"]])

        except Exception as e:
            st.error(f"Dynamic model failed: {e}")
