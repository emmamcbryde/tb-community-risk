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

    # === Intervention effects with coverage ===
    # === Intervention effects with coverage ===
    coverage = inputs.get("coverage", 1.0)

    if testing != "None" and treatment != "None":
        full_effect = 0.7   # 30% reduction if everyone covered
    elif testing != "None":
        full_effect = 0.85  # 15% reduction
    elif treatment != "None":
        full_effect = 0.8   # 20% reduction
    else:
        full_effect = 1.0   # no change

    # Adjust for coverage
    intervention_factor = 1 - coverage * (1 - full_effect)

    # Also produce a baseline (no-intervention) series for comparison
    baseline_factor = 1.0


    # Adjust for coverage
    intervention_factor = 1 - coverage * (1 - full_effect)

    # Also produce a baseline (no-intervention) series for comparison
    baseline_factor = 1.0
 


   # === Annual simulation ===
    years = np.arange(0, T + 1)
    incidence_baseline = base_incidence * baseline_risk * baseline_factor * np.exp(-0.05 * years)
    incidence_strategy = base_incidence * baseline_risk * intervention_factor * np.exp(-0.05 * years)

    cases_baseline = N * incidence_baseline / 1e5
    cases_strategy = N * incidence_strategy / 1e5

    # Average mortality lookups (fallback defaults)
    avg_tb_mort = np.mean(list(tbmort.values())) if tbmort else 0.05

    tb_deaths_baseline = cases_baseline * avg_tb_mort
    tb_deaths_strategy = cases_strategy * avg_tb_mort

    # === Assemble DataFrame ===
    df = pd.DataFrame({
        "Year": years,
        "Incidence_baseline": incidence_baseline,
        "Incidence_strategy": incidence_strategy,
        "Cases_baseline": cases_baseline,
        "Cases_strategy": cases_strategy,
        "TB_Deaths_baseline": tb_deaths_baseline,
        "TB_Deaths_strategy": tb_deaths_strategy,
    })
    summary = {
        "Peak incidence": np.max(incidence_strategy),
        "Total cases": np.sum(cases_strategy),
        "Total TB deaths": np.sum(tb_deaths_strategy),
        "Peak incidence (baseline)": np.max(incidence_baseline),
        "Total cases (baseline)": np.sum(cases_baseline),
        "Total TB deaths (baseline)": np.sum(tb_deaths_baseline),
    }


    return df, summary
