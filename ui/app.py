import os, sys
# --- ensure project root is in Python path ---
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
os.chdir(project_root)
print("Working directory set to:", os.getcwd())

import streamlit as st
import pandas as pd
import numpy as np
from engine.model import simulate_community


# --- Utility functions ---
def fetch_age_dist_worldbank(country_code="AUS", year=None):
    """
    Try to fetch age distribution from OWID dataset hosted on GitHub.
    """
    owid_url = "https://raw.githubusercontent.com/owid/owid-datasets/master/datasets/Population%20by%20age%20group%20(UN%20WPP)/Population%20by%20age%20group%20(UN%20WPP).csv"
    df = pd.read_csv(owid_url)
    df = df[df["Entity"].str.contains(country_code, case=False, na=False)]
    if df.empty:
        raise ValueError(f"No age data found for {country_code}.")
    if year:
        df = df[df["Year"] == int(year)]
    else:
        df = df[df["Year"] == df["Year"].max()]
    cols = [c for c in df.columns if "population" in c.lower()]
    sub = df[["Entity", "Year"] + cols].melt(id_vars=["Entity", "Year"], 
                                             var_name="AgeGroup", 
                                             value_name="Population")
    sub = sub.groupby("AgeGroup", as_index=False)["Population"].sum()
    sub["Proportion"] = sub["Population"] / sub["Population"].sum()
    return sub[["AgeGroup", "Proportion"]]


def load_fallback_csv(country_code="AUS", file_path="data/population_age_latest.csv"):
    """
    Load fallback age distribution from your local CSV file.
    """
    df = pd.read_csv(file_path)
    if "iso_code" in df.columns:
        df = df[df["iso_code"].str.upper() == country_code.upper()]
    if df.empty:
        raise ValueError(f"No data for {country_code} in local CSV.")
    df = df[df["year"] == df["year"].max()]
    df = df.groupby("age", as_index=False)["population"].sum()
    df["Proportion"] = df["population"] / df["population"].sum()
    return df.rename(columns={"age": "AgeGroup"})[["AgeGroup", "Proportion"]]


def default_age_distribution():
    """
    Global default (UN 2023 est.)
    """
    age_groups = [
        "0-4","5-9","10-14","15-19","20-24","25-29",
        "30-34","35-39","40-44","45-49","50-54",
        "55-59","60-64","65-69","70-74","75-79","80+"
    ]
    props = [
        0.088,0.086,0.085,0.084,0.083,0.080,
        0.078,0.075,0.071,0.066,0.061,
        0.055,0.049,0.041,0.033,0.025,0.020
    ]
    df = pd.DataFrame({"AgeGroup": age_groups, "Proportion": props})
    df["Proportion"] /= df["Proportion"].sum()
    return df


def load_custom_csv(file_obj):
    """
    Allow user-defined CSV upload.
    """
    df = pd.read_csv(file_obj)
    if "Count" in df.columns:
        df["Proportion"] = df["Count"] / df["Count"].sum()
    elif "Proportion" in df.columns:
        df["Proportion"] /= df["Proportion"].sum()
    else:
        raise ValueError("CSV must have 'Count' or 'Proportion' column.")
    return df[["AgeGroup", "Proportion"]]


# --- Streamlit layout ---
st.set_page_config(page_title="TB Community Risk Model", layout="wide")

st.title("🧮 Community-level TB Risk Model")
st.markdown("Estimate TB risk based on community characteristics and interventions.")

# --- Sidebar Inputs ---
st.sidebar.header("Community Parameters")

population = st.sidebar.number_input("Population size", min_value=100, value=10000, step=100)
indigenous_pct = st.sidebar.slider("Indigenous population (%)", 0, 100, 60)
smoker_pct = st.sidebar.slider("Smoker population (%)", 0, 100, 30)
diabetes_pct = st.sidebar.slider("Diabetes (%)", 0, 100, 10)
renal_pct = st.sidebar.slider("Renal impairment (%)", 0, 100, 5)
immune_pct = st.sidebar.slider("Immunosuppressed (%)", 0, 100, 3)
recent_infection = st.sidebar.slider("Recent infection (%)", 0, 50, 5)
user_incidence = st.sidebar.number_input("Baseline TB incidence (per 100k/year)", 0, 500, 30)
time_horizon = st.sidebar.slider("Time horizon (years)", 1, 30, 20)

st.sidebar.header("Intervention Strategy")
testing = st.sidebar.selectbox("Testing method", ["None", "TST", "IGRA"])
treatment = st.sidebar.selectbox("Treatment regimen", ["None", "3HP", "4R", "6H", "9H"])

# --- Age distribution section ---
st.sidebar.header("Population age structure")

method = st.sidebar.radio(
    "Choose age distribution source:",
    ["Fetch from World Bank/OWID", "Use local CSV (data/population_age_latest.csv)", "Default global", "Upload custom CSV"]
)

if method == "Fetch from World Bank/OWID":
    country = st.sidebar.text_input("Enter country name or ISO3 code", "AUS")
    year = st.sidebar.number_input("Year (optional)", min_value=1950, max_value=2050, value=2023)
    try:
        age_df = fetch_age_dist_worldbank(country, year)
        st.success(f"✅ Age distribution fetched for {country} ({year})")
    except Exception as e:
        st.warning(f"⚠️ Could not fetch live data: {e}\nLoading local fallback.")
        try:
            age_df = load_fallback_csv(country)
            st.info("Loaded fallback CSV data.")
        except Exception as e2:
            st.error(f"❌ Local fallback failed: {e2}. Using global default.")
            age_df = default_age_distribution()

elif method == "Use local CSV (data/population_age_latest.csv)":
    country = st.sidebar.text_input("Enter country ISO3 code", "AUS")
    try:
        age_df = load_fallback_csv(country)
        st.success(f"✅ Loaded local CSV for {country}.")
    except Exception as e:
        st.warning(f"⚠️ Local CSV unavailable or invalid: {e}. Using global default.")
        age_df = default_age_distribution()

elif method == "Upload custom CSV":
    file = st.sidebar.file_uploader("Upload your own age distribution (AgeGroup, Count/Proportion)", type=["csv"])
    if file:
        try:
            age_df = load_custom_csv(file)
            st.success("✅ Custom age distribution loaded.")
        except Exception as e:
            st.error(f"❌ Error in custom file: {e}. Using global default.")
            age_df = default_age_distribution()
    else:
        age_df = default_age_distribution()

else:
    st.info("Using default global UN 2023 distribution.")
    age_df = default_age_distribution()

# Show the loaded data
st.subheader("Population Age Distribution")
st.dataframe(age_df)
age_df["Count"] = age_df["Proportion"] * population

# --- Run Simulation ---
st.sidebar.header("Run Simulation")
if st.sidebar.button("Simulate Community"):
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
        "user_incidence": user_incidence,
    }

    try:
        df, summary = simulate_community(inputs, file_path="data/parameters.xlsx")

        st.subheader("📊 Results Summary")
        c1, c2, c3 = st.columns(3)
        c1.metric("Peak incidence (per 100,000)", f"{summary['Peak incidence']:.1f}")
        c2.metric("Total TB cases", f"{summary['Total cases']:.0f}")
        c3.metric("Total TB deaths", f"{summary['Total TB deaths']:.0f}")

        st.subheader("📈 Annual Incidence Over Time")
        st.line_chart(df.set_index("Year")[["Incidence_per_100k"]])

        st.success("✅ Simulation complete!")

    except Exception as e:
        st.error(f"Simulation failed: {e}")
