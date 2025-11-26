import os, sys
import streamlit as st
import pandas as pd
import numpy as np
import altair as alt

DEBUG = True

# --- Fix working directory ---
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
os.chdir(project_root)

print("Working directory set to:", os.getcwd())

# --- Engine imports ---
from engine.model import simulate_community
from engine.dynamic_model import simulate_dynamic_ltbi  # <-- YOU MUST HAVE THIS FILE
from engine.infection_backcast import calc_ari_from_incidence, infection_prob_by_age


# --- STUBS for missing functions (prevent NameErrors) ---
def default_age_distribution():
    """Temporary fallback until real implementation supplied."""
    ages = list(range(0, 101))
    prop = np.array([1 / 101] * 101)
    return pd.DataFrame({"AgeGroup": ages, "Proportion": prop})


def load_custom_csv(file):
    df = pd.read_csv(file)
    if "AgeGroup" not in df.columns:
        raise ValueError("Custom age file must contain AgeGroup column.")
    if "Count" in df.columns:
        df["Proportion"] = df["Count"] / df["Count"].sum()
    elif "Proportion" not in df.columns:
        raise ValueError("Custom age file must contain Count or Proportion.")
    return df[["AgeGroup", "Proportion"]]


def fetch_age_dist_worldbank(country_code, year):
    """Placeholder until real API implemented."""
    raise NotImplementedError("World Bank API not implemented yet.")


# --- Fallback loader already provided ---
def load_fallback_csv(country_code="AUS", file_path="data/population_age_latest.csv"):
    df = pd.read_csv(file_path)
    if "iso_code" in df.columns:
        df = df[df["iso_code"].str.upper() == country_code.upper()]
    if df.empty:
        raise ValueError(f"No data for {country_code} in fallback CSV.")
    df = df[df["year"] == df["year"].max()]
    df = df.groupby("age", as_index=False)["population"].sum()
    df["Proportion"] = df["population"] / df["population"].sum()
    return df.rename(columns={"age": "AgeGroup"})[["AgeGroup", "Proportion"]]


# --- Streamlit page setup ---
st.set_page_config(page_title="TB Community Risk Model", layout="wide")
st.title("🧮 Community-level TB Risk Model")
st.markdown(
    "Estimate tuberculosis risk based on community characteristics and interventions."
)

# --- Sidebar: Community Parameters ---
st.sidebar.header("Community Parameters")
population = st.sidebar.number_input(
    "Population size", min_value=100, value=10000, step=100
)
smoker_pct = st.sidebar.slider("Smoker population (%)", 0, 100, 30)
diabetes_pct = st.sidebar.slider("Diabetes (%)", 0, 100, 10)
renal_pct = st.sidebar.slider("Renal impairment (%)", 0, 100, 5)
immune_pct = st.sidebar.slider("Immunosuppressed (%)", 0, 100, 3)
recent_infection = st.sidebar.slider("Recent infection (%)", 0, 50, 5)
user_incidence = st.sidebar.number_input(
    "Baseline TB incidence (per 100k/yr)", 0, 500, 30
)
# ARI and LTBI backcast
hist_pattern = st.sidebar.selectbox(
    "Historical incidence pattern (last 20 years)",
    ["Constant", "Declining 3%/yr", "Increasing 3%/yr"],
)
time_horizon = st.sidebar.slider("Time horizon (years)", 1, 30, 20)

# --- Sidebar: Intervention Strategy ---
st.sidebar.header("Intervention Strategy")
testing = st.sidebar.selectbox("Testing method", ["None", "TST", "IGRA"])
treatment = st.sidebar.selectbox(
    "Treatment regimen", ["None", "1HP", "3HP", "4R", "6H", "9H"]
)
rollout_years = st.sidebar.slider("Treatment rollout duration (years)", 1, 5, 3)

coverage_testing = st.sidebar.slider("Testing coverage", 0.0, 1.0, 0.5)
coverage_treatment = st.sidebar.slider("Treatment coverage", 0.0, 1.0, 0.7)
ltbi_prev_input = st.sidebar.slider(
    "LTBI prevalence (population average)", 0.0, 1.0, 0.02
)
# TB case detection parameters (in months)
pre_det_months = st.sidebar.number_input(
    "Mean time to TB diagnosis pre-intervention (months)",
    min_value=0.5, max_value=60.0, value=12.0, step=0.5
)
post_det_months = st.sidebar.number_input(
    "Mean time to TB diagnosis post-intervention (months)",
    min_value=0.5, max_value=60.0, value=6.0, step=0.5
)

delta_pre = 12.0 / pre_det_months   # per year
delta_post = 12.0 / post_det_months



# --- Sidebar: Model Mode ---
st.sidebar.header("Model mode")
model_mode = st.sidebar.radio(
    "Choose model type:", ["Static incidence model", "Dynamic LTBI→TB model"]
)

# --- Sidebar: Age distribution ---
st.sidebar.header("Population age structure")
age_method = st.sidebar.radio(
    "Choose age distribution source:",
    [
        "Fetch from World Bank/OWID",
        "Use local CSV",
        "Default global",
        "Upload custom CSV",
    ],
)

# Load age distribution
if age_method == "Fetch from World Bank/OWID":
    country = st.sidebar.text_input("Enter ISO3 code", "AUS")
    year = st.sidebar.number_input("Year", 1950, 2050, 2023)

    try:
        # ⭐ Move EVERYTHING inside the try block
        age_df = fetch_age_dist_worldbank(country, year)
        st.success(f"Fetched age distribution for {country}.")
    except Exception as e:
        st.error(f"Could not fetch OWID/World Bank data: {e}")

        if DEBUG:
            import traceback

           # st.code(traceback.format_exc())

        # ⭐ GUARANTEED fallback path
        try:
            age_df = load_fallback_csv(country)
            st.info("Loaded local CSV fallback.")
        except:
            age_df = default_age_distribution()
            st.warning("Using default global distribution.")


elif age_method == "Use local CSV":
    country = st.sidebar.text_input("Enter ISO3 code", "AUS")
    try:
        age_df = load_fallback_csv(country)
        st.success("Loaded fallback age data.")
    except Exception as e:
        st.error(f"Local CSV failed: {e}")
        age_df = default_age_distribution()

elif age_method == "Upload custom CSV":
    file = st.sidebar.file_uploader("Upload age distribution CSV", type="csv")
    if file:
        try:
            age_df = load_custom_csv(file)
            st.success("Custom age distribution loaded.")
        except Exception as e:
            st.error(f"Custom CSV failed: {e}")
            age_df = default_age_distribution()
    else:
        age_df = default_age_distribution()

else:
    age_df = default_age_distribution()

# Display age distribution
st.subheader("Population Age Distribution")
st.dataframe(age_df)

# Add counts from proportions
age_df["Count"] = age_df["Proportion"] * population
age_counts = dict(zip(age_df["AgeGroup"], age_df["Count"]))

# --- Simulation Button ---
st.sidebar.header("Run Simulation")

if st.sidebar.button("Simulate Community"):
    st.info("Running simulation...")

    # Gather inputs
    inputs = {
        "population": population,
        "smoker_pct": smoker_pct,
        "diabetes_pct": diabetes_pct,
        "renal_pct": renal_pct,
        "immune_pct": immune_pct,
        "recent_infection": recent_infection,
        "time_horizon": time_horizon,
        "testing": testing,
        "treatment": treatment,
        "user_incidence": user_incidence,
        "coverage_testing": coverage_testing,
        "coverage_treatment": coverage_treatment,
        "ltbi_prev": ltbi_prev_input,
        "rollout_years": rollout_years,
        "delta_pre": delta_pre,
        "delta_post": delta_post,
    }


    max_age = max(age_counts.keys())
    years_back = range(0, max_age + 1)  # go back as far as oldest age group

    if hist_pattern == "Constant":
        inc_hist = {-k: user_incidence for k in years_back}

    elif hist_pattern == "Declining 3%/yr":
        inc_hist = {-k: user_incidence * (1.03**k) for k in years_back}

    elif hist_pattern == "Increasing 3%/yr":
        inc_hist = {-k: user_incidence * (0.97**k) for k in years_back}

    ari_hist = calc_ari_from_incidence(inc_hist)
    ages = list(age_counts.keys())
    ltbi_prev = infection_prob_by_age(ages, ari_hist)
    # Take the most recnt ARI (year 0) as the baseline FOI for the dynamic model
    # inc_hist uses keys 0, -1, -2, ... so max key is "0" = present
    base_ari_from_history = ari_hist[max(ari_hist.keys())]

    # Store in inputs so the dynamic model can use it
    inputs["base_ari"] = base_ari_from_history

    try:
        # --- Model selection ---
        if model_mode == "Static incidence model":
            df, summary = simulate_community(inputs, file_path="data/parameters.xlsx")
        else:
            df, summary = simulate_dynamic_ltbi(
                age_counts=age_counts,
                ltbi_prev=ltbi_prev,
                inputs=inputs,
                file_path="data/parameters.xlsx",
            )
        st.write("Dynamic model output preview:", df.head())

        # --- Results summary ---
        st.subheader("📊 Results Summary")
        # Convert the summary into a table for cleaner presentation
        summary_df = pd.DataFrame.from_dict(summary, orient='index', columns=["Value"])
        st.write("📊 Results Summary", summary_df)

        # --- Incidence plot ---
        df["Year"] = pd.to_numeric(df["Year"])
        df["Incidence_per_100k"] = pd.to_numeric(df["Incidence_per_100k"])

        hover = alt.selection_single(
            fields=["Year"], nearest=True, on="mouseover", empty="none"
        )

        chart = (
            alt.Chart(df)
            .mark_line()
            .encode(x="Year:Q", y="Incidence_per_100k:Q", color="Scenario:N")
        )

        points = (
            chart.mark_circle(size=120)
            .encode(opacity=alt.condition(hover, alt.value(1), alt.value(0)))
            .add_selection(hover)
        )

        rule = (
            alt.Chart(df)
            .mark_rule(color="gray")
            .encode(
                x="Year:Q", opacity=alt.condition(hover, alt.value(0.6), alt.value(0))
            )
        )

        tooltips = chart.mark_rule().encode(
            tooltip=["Year:Q", "Scenario:N", "Incidence_per_100k:Q"],
            opacity=alt.condition(hover, alt.value(1), alt.value(0)),
        )

        final_chart = (chart + points + rule + tooltips).interactive()
        st.subheader("📈 Annual Incidence Over Time")
        st.altair_chart(final_chart, use_container_width=True)

        # --- LTBI by age ---
        st.subheader("📉 LTBI Prevalence by Age")
        ltbi_df = pd.DataFrame(
            {"Age": ages, "LTBI Prevalence": [ltbi_prev[a] for a in ages]}
        )
        st.area_chart(ltbi_df.set_index("Age"))

        st.success("Simulation complete!")

    except Exception as e:
        st.error(f"Simulation failed: {e}")
