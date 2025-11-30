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
from engine.dynamic_model import simulate_dynamic_ltbi
from engine.infection_backcast import (
    calc_ari_from_incidence,
    infection_prob_by_age_split,
)

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


# --- Fallback loader ---
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

hist_pattern = st.sidebar.selectbox(
    "Historical incidence pattern (for LTBI backcast)",
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
infection_per_case = st.sidebar.number_input(
    "Transmission rate β (infections per active TB case per year)",
    min_value=0.0,
    max_value=50.0,
    value=8.0,     # default you choose
    step=0.1
)

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

    # OWID/World Bank not implemented – use fallback silently
    st.info("Using fallback age distribution (OWID/World Bank not implemented).")

    try:
        age_df = load_fallback_csv(country)
        st.success("Loaded fallback age distribution.")
    except Exception:
        age_df = default_age_distribution()
        st.warning("Using default global age distribution.")


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
st.subheader("Population Age Distribution (5-Year Bins)")

# --- Ensure numeric types ---
age_df["AgeGroup"] = pd.to_numeric(age_df["AgeGroup"], errors="coerce")
age_df["Proportion"] = pd.to_numeric(age_df["Proportion"], errors="coerce")
age_df["Count"] = pd.to_numeric(age_df.get("Count", age_df["Proportion"] * population), errors="coerce")
age_df = age_df.dropna(subset=["AgeGroup", "Count"])

# --- Create 5-year bins ---
bin_edges = list(range(0, 105, 5))    # 0,5,10,...100
bin_labels = [f"{bin_edges[i]}–{bin_edges[i+1]-1}" for i in range(len(bin_edges)-1)]
bin_labels.append("100+")

age_df["AgeBin"] = pd.cut(
    age_df["AgeGroup"],
    bins=bin_edges + [200],    # allow ages above 100
    labels=bin_labels,
    right=False
)

# --- Aggregate ---
bin_df = (
    age_df.groupby("AgeBin", as_index=False)["Count"]
    .sum()
)

# ----- KEY FIX: Add a numeric sort index -----
bin_order_map = {label: i for i, label in enumerate(bin_labels)}
bin_df["AgeBin_order"] = bin_df["AgeBin"].map(bin_order_map)

# Sort by order
bin_df = bin_df.sort_values("AgeBin_order")

# --- Plot ---
hist_chart = (
    alt.Chart(bin_df)
    .mark_bar()
    .encode(
        x=alt.X("AgeBin:N",
                sort=list(bin_df["AgeBin"]),        # force correct ordering
                title="Age Group (5-year bins)"
        ),
        y=alt.Y("Count:Q", title="Population"),
        tooltip=["AgeBin:N", "Count:Q"]
    )
)

st.altair_chart(hist_chart, use_container_width=True)



# --- Simulation Button ---
st.sidebar.header("Run Simulation")

if st.sidebar.button("Simulate Community"):
    st.info("Running simulation...")

    # 1. Build historical incidence series
    max_age = max(age_counts.keys())
    years_back = range(0, max_age + 1)

    if hist_pattern == "Constant":
        inc_hist = {-k: user_incidence for k in years_back}
    elif hist_pattern == "Declining 3%/yr":
        inc_hist = {-k: user_incidence * (1.03 ** k) for k in years_back}
    else:  # Increasing 3%/yr
        inc_hist = {-k: user_incidence * (0.97 ** k) for k in years_back}

    # 2. Convert incidence → ARI → LTBI splits
    ari_hist = calc_ari_from_incidence(inc_hist)
    ages = sorted(age_counts.keys())

    ltbi_ever, ltbi_recent, ltbi_remote = infection_prob_by_age_split(
        ages,
        ari_hist,
        window_recent=5
    )

    ltbi_prev = ltbi_ever  # for dynamic model & plotting

    # 3. Build inputs
    inputs = {
        "population": population,
        "smoker_pct": smoker_pct,
        "diabetes_pct": diabetes_pct,
        "renal_pct": renal_pct,
        "immune_pct": immune_pct,
        "time_horizon": time_horizon,
        "testing": testing,
        "treatment": treatment,
        "user_incidence": user_incidence,
        "coverage_testing": coverage_testing,
        "coverage_treatment": coverage_treatment,
        "rollout_years": rollout_years,
        "delta_pre": delta_pre,
        "delta_post": delta_post,
        "age_counts": age_counts,
        "ltbi_ever_by_age": ltbi_ever,
        "ltbi_recent_by_age": ltbi_recent,
        "ltbi_remote_by_age": ltbi_remote,
            # --- Build a baseline-only inputs copy for comparison ---

    }
    baseline_inputs = inputs.copy()
    baseline_inputs["coverage_testing"] = 0.0
    baseline_inputs["coverage_treatment"] = 0.0
    baseline_inputs["testing"] = "None"
    baseline_inputs["treatment"] = "None"
    baseline_inputs["secondary_cases_per_index"] = 0.0  # beta = 0
    baseline_inputs["time_horizon"] = 1  # just enough to get year 0

    try:
        # --- Model selection ---
        if model_mode == "Static incidence model":
            df, summary = simulate_community(inputs, file_path="data/parameters.xlsx")
        else:
            df, summary = simulate_dynamic_ltbi(
                age_counts=age_counts,
                beta=inputs.get("infection_per_case", 1.0),
                inputs=inputs,
                file_path="data/parameters.xlsx",
            )

        st.write("Model output preview:", df.head())

        # --- ltbi_prev summary ---
        st.subheader("📊 ltbi_prev Summary")
        summary_df = pd.DataFrame.from_dict(summary, orient="index", columns=["Value"])
        summary_df = summary_df.round({"Value": 1})
        st.write(summary_df)
        
        # --- Baseline-only comparison: static vs dynamic at Year 0 ---
        df_static_base, _ = simulate_community(
            baseline_inputs, file_path="data/parameters.xlsx"
        )
        df_dynamic_base, _ = simulate_dynamic_ltbi(
            age_counts=age_counts,
            ltbi_prev=ltbi_prev,  # ever-infected dict
            inputs=baseline_inputs,
            beta=baseline_inputs.get("secondary_cases_per_index", 0.0),  # beta=0 baseline
            file_path="data/parameters.xlsx",
        )

        # Extract Year 0 incidence for baseline branches
        static0 = df_static_base.loc[
            (df_static_base["Scenario"] == "Baseline") &
            (df_static_base["Year"] == 0),
            "Incidence_per_100k"
        ].iloc[0]

        dynamic0 = df_dynamic_base.loc[
            (df_dynamic_base["Scenario"] == "Dynamic_baseline") &
            (df_dynamic_base["Year"] == 0),
            beta==baseline_inputs.get("secondary_cases_per_index", 0.0), 
            "Incidence_per_100k"
        ].iloc[0]

        st.subheader("🔍 Baseline Year 0 Incidence (No FOI, No Interventions)")
        st.write(
            pd.DataFrame(
                {
                    "Model": ["Static (hazard)", "Dynamic (mechanistic)"],
                    "Incidence_per_100k": [round(static0, 1), round(dynamic0, 1)],
                }
            )
        )

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
        # Truncate to age <= 60
        ltbi_df = ltbi_df[ltbi_df["Age"] <= 60]
        st.area_chart(ltbi_df.set_index("Age"))

        st.success("Simulation complete!")

    except Exception as e:
        st.error(f"Simulation failed: {e}")
