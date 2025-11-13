
import pandas as pd
import requests
import io

# ------------------------------
# 1. Fetch from World Bank
# ------------------------------
def fetch_age_dist_worldbank(country_code="AUS", year=None):
    """
    Fetch age distribution for a given ISO2/3 country code from the World Bank API.
    Returns a pandas DataFrame with columns: AgeGroup, Population, Proportion.
    """

    # Base dataset: SP.POP.0004.FE (female 0–4), SP.POP.0004.MA (male 0–4), etc.
    # World Bank’s structure uses multiple indicators; we’ll get all ages 0-4 ... 80+ in one call.

    url = f"https://api.worldbank.org/v2/en/country/{country_code}/indicator/SP.POP.TOTL?format=json"
    try:
        # Get total population first (for normalisation)
        r = requests.get(url)
        total = r.json()[1][0]["value"]
    except Exception as e:
        raise ConnectionError(f"Could not fetch population for {country_code}: {e}")

    # Fetch age-structured population from OWID as a reliable fallback
    owid_url = "https://raw.githubusercontent.com/owid/owid-datasets/master/datasets/Population%20by%20age%20group%20(UN%20WPP)/Population%20by%20age%20group%20(UN%20WPP).csv"
    df = pd.read_csv(owid_url)
    df = df[df["Entity"].str.contains(country_code, case=False, na=False)]
    if year:
        df = df[df["Year"] == int(year)]
    else:
        df = df[df["Year"] == df["Year"].max()]

    # Simplify and normalise
    cols = [c for c in df.columns if "population" in c.lower()]
    sub = df[["Entity", "Year"] + cols].melt(
        id_vars=["Entity", "Year"], var_name="AgeGroup", value_name="Population"
    )
    sub = sub.groupby("AgeGroup", as_index=False)["Population"].sum()
    sub["Proportion"] = sub["Population"] / sub["Population"].sum()
    sub = sub[["AgeGroup", "Population", "Proportion"]]
    return sub
    
def fetch_age_dist(country_code="AUS", year=None):
    try:
        # attempt live fetch (as earlier)
        age_df = fetch_age_dist_worldbank(country_code, year)
        return age_df
    except Exception as e:
        import warnings
        warnings.warn(f"Live fetch failed: {e}. Using local cached CSV.")
        csv_path = "data/population_age_latest.csv"
        df = pd.read_csv(csv_path)
        sub = df[(df["iso_code"] == country_code.upper()) & 
                 (df["year"] == (year if year is not None else df["year"].max()))]
        sub = sub.copy()
        sub["Proportion"] = sub["population"] / sub["population"].sum()
        return sub[["age", "Proportion"]].rename(columns={"age":"AgeGroup"})


# ------------------------------
# 2. Manual upload handler
# ------------------------------
def load_custom_age_distribution(file_obj):
    """
    Read uploaded CSV with columns: AgeGroup, Count or Proportion.
    Returns DataFrame with normalised 'Proportion'.
    """
    df = pd.read_csv(file_obj)
    if "Count" in df.columns:
        df["Proportion"] = df["Count"] / df["Count"].sum()
    elif "Proportion" in df.columns:
        df["Proportion"] = df["Proportion"] / df["Proportion"].sum()
    else:
        raise ValueError("CSV must contain either 'Count' or 'Proportion' column.")
    return df[["AgeGroup", "Proportion"]]


# ------------------------------
# 3. Default global fallback
# ------------------------------
def default_age_distribution():
    """
    Returns an approximate global age structure (UN 2023, 5-year bands).
    """
    age_groups = [
        "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
        "30-34", "35-39", "40-44", "45-49", "50-54",
        "55-59", "60-64", "65-69", "70-74", "75-79", "80+"
    ]
    global_props = [
        0.088, 0.086, 0.085, 0.084, 0.083, 0.080,
        0.078, 0.075, 0.071, 0.066, 0.061,
        0.055, 0.049, 0.041, 0.033, 0.025, 0.020
    ]
    df = pd.DataFrame({"AgeGroup": age_groups, "Proportion": global_props})
    df["Proportion"] = df["Proportion"] / df["Proportion"].sum()
    return df