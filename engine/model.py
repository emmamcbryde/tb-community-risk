import numpy as np
import pandas as pd
from engine.params import extract_core_parameters


def simulate_community(inputs, file_path="data/parameters.xlsx"):
    """
    Simulate TB incidence and mortality in a mixed community using parameters from Excel.

    Parameters
    ----------
    inputs : dict
        Dictionary from Streamlit app, e.g.:
        {
          'population': 10000,
          'indigenous_pct': 60,
          'smoker_pct': 30,
          'diabetes_pct': 10,
          'renal_pct': 5,
          'immune_pct': 3,
          'recent_infection': 5,
          'time_horizon': 20,
          'testing': 'IGRA',
          'treatment': '3HP'
        }

    file_path : str
        Path to parameters.xlsx (default: data/parameters.xlsx)

    Returns
    -------
    df : pandas.DataFrame
        Yearly values for incidence, cases, and TB deaths.
    summary : dict
        Summary statistics (peak incidence, total cases, total deaths).
    """

    # === Load parameters ===
    params = extract_core_parameters(file_path)

    # Access key tables (fall back to empty dicts if missing)
    mort = params.get("mortality", {})
    tbmort = params.get("tb_mortality", {})
    comorbid = params.get("comorbidities", {})
    community_risk = params.get("community_risk", {})
    cascade = params.get("cascade", {})

    # === Inputs ===
    N = inputs.get("population", 10000)
    indigenous = inputs.get("indigenous_pct", 0) / 100
    smoker = inputs.get("smoker_pct", 0) / 100
    diabetes = inputs.get("diabetes_pct", 0) / 100
    renal = inputs.get("renal_pct", 0) / 100
    immune = inputs.get("immune_pct", 0) / 100
    recent_inf = inputs.get("recent_infection", 0) / 100
    T = inputs.get("time_horizon", 20)
    testing = inputs.get("testing", "None")
    treatment = inputs.get("treatment", "None")

    # === Core relative risks ===
    base_incidence = 200  # baseline per 100k
    rr_indigenous = community_risk.get("Indigenous", 2.0)
    rr_smoking = comorbid.get("Smoking", 1.5)
    rr_diabetes = comorbid.get("Diabetes", 1.8)
    rr_renal = comorbid.get("Renal disease", 2.0)
    rr_immune = comorbid.get("Immunosuppressed", 4.0)
    rr_recent = 2.0

    # === Combine relative risks ===
    baseline_risk = (
        1
        + indigenous * (rr_indigenous - 1)
        + smoker * (rr_smoking - 1)
        + diabetes * (rr_diabetes - 1)
        + renal * (rr_renal - 1)
        + immune * (rr_immune - 1)
        + recent_inf * (rr_recent - 1)
    )

    # === Intervention effects ===
    intervention_factor = 1.0
    if testing != "None" and treatment != "None":
        intervention_factor = 0.7  # full test+treat
    elif testing != "None":
        intervention_factor = 0.85
    elif treatment != "None":
        intervention_factor = 0.8

    # === Annual simulation ===
    years = np.arange(0, T + 1)
    incidence = base_incidence * baseline_risk * intervention_factor * np.exp(-0.05 * years)
    cases = N * incidence / 1e5

    # Average mortality lookups (fallback defaults)
    avg_mortality = np.mean(list(mort.values())) if mort else 0.0007
    avg_tb_mort = np.mean(list(tbmort.values())) if tbmort else 0.05

    natural_deaths = N * avg_mortality * (years / T)
    tb_deaths = cases * avg_tb_mort

    # === Assemble DataFrame ===
    df = pd.DataFrame({
        "Year": years,
        "Incidence_per_100k": incidence,
        "Cases": cases,
        "TB_Deaths": tb_deaths,
        "Natural_Deaths": natural_deaths
    })

    # === Summary ===
    summary = {
        "Peak incidence": np.max(incidence),
        "Total cases": np.sum(cases),
        "Total TB deaths": np.sum(tb_deaths),
        "Total natural deaths": np.sum(natural_deaths)
    }

    return df, summary
