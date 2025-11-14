import numpy as np
import pandas as pd
from engine.params import extract_core_parameters

import numpy as np
import pandas as pd

def simulate_community(inputs, file_path="data/parameters.xlsx"):
    """
    Simulate baseline and intervention TB incidence trajectories.
    Always returns a valid df and summary dict.
    """

    # Inputs
    N = inputs.get("population", 10000)
    T = inputs.get("time_horizon", 20)
    base_incidence = float(inputs.get("user_incidence", 30))
    testing = inputs.get("testing", "None")
    treatment = inputs.get("treatment", "None")

    # Risk factors
    indigenous = inputs.get("indigenous_pct", 0) / 100
    smoker     = inputs.get("smoker_pct", 0) / 100
    diabetes   = inputs.get("diabetes_pct", 0) / 100
    renal      = inputs.get("renal_pct", 0) / 100
    immune     = inputs.get("immune_pct", 0) / 100
    recent_inf = inputs.get("recent_infection", 0) / 100

    # Multipliers
    rr_indigenous = 2.0
    rr_smoker     = 1.5
    rr_diabetes   = 1.8
    rr_renal      = 2.0
    rr_immune     = 4.0
    rr_recent_inf = 2.0

    # Combined baseline risk
    baseline_risk = (
        1
        + indigenous * (rr_indigenous - 1)
        + smoker     * (rr_smoker - 1)
        + diabetes   * (rr_diabetes - 1)
        + renal      * (rr_renal - 1)
        + immune     * (rr_immune - 1)
        + recent_inf * (rr_recent_inf - 1)
    )

    # Intervention factor
    intervention_factor = 1.0
    if testing != "None" and treatment != "None":
        intervention_factor = 0.70
    elif testing != "None":
        intervention_factor = 0.85
    elif treatment != "None":
        intervention_factor = 0.80

    # Time steps
    years = np.arange(0, T + 1)

    # Baseline (slower decline)
    baseline_incidence = base_incidence * baseline_risk * np.exp(-0.02 * years)
    baseline_cases     = N * baseline_incidence / 100_000

    # Intervention scenario (faster decline)
    intervention_incidence = base_incidence * baseline_risk * intervention_factor * np.exp(-0.05 * years)
    intervention_cases     = N * intervention_incidence / 100_000

    # Build DF
    df = pd.DataFrame({
        "Year": np.tile(years, 2),
        "Scenario": np.repeat(["Baseline", "Intervention"], len(years)),
        "Incidence_per_100k": np.concatenate([baseline_incidence, intervention_incidence]),
        "Cases": np.concatenate([baseline_cases, intervention_cases])
    })

    # Summary — IMPORTANT: use consistent key names
    summary = {
        "Peak_incidence_baseline": float(np.max(baseline_incidence)),
        "Peak_incidence_intervention": float(np.max(intervention_incidence)),
        "Total_cases_baseline": float(np.sum(baseline_cases)),
        "Total_cases_intervention": float(np.sum(intervention_cases)),
    }

    return df, summary
