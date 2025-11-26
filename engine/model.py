import numpy as np
import pandas as pd
from engine.params import extract_core_parameters
from engine.intervention import compute_full_effect


def simulate_community(inputs, file_path="data/parameters.xlsx"):

    # Inputs
    N = inputs.get("population", 10000)
    T = inputs.get("time_horizon", 20)
    base_incidence = float(inputs.get("user_incidence", 30))
    testing = inputs.get("testing", "None")
    treatment = inputs.get("treatment", "None")
    rollout_years = float(inputs.get("rollout_years", 3.0))

    # Risk factors
    smoker = inputs.get("smoker_pct", 0) / 100
    diabetes = inputs.get("diabetes_pct", 0) / 100
    renal = inputs.get("renal_pct", 0) / 100
    immune = inputs.get("immune_pct", 0) / 100
    recent_inf = inputs.get("recent_infection", 0) / 100

    # Multipliers
    rr_smoker = 1.5
    rr_diabetes = 1.8
    rr_renal = 2.0
    rr_immune = 4.0
    rr_recent_inf = 10.0

    # Combined baseline risk
    baseline_risk = (
        1
        + smoker * (rr_smoker - 1)
        + diabetes * (rr_diabetes - 1)
        + renal * (rr_renal - 1)
        + immune * (rr_immune - 1)
        + recent_inf * (rr_recent_inf - 1)
    )

    # --- rollout function (user-defined years) ---
    def rollout_factor(t):
        if t <= 0:
            return 0.0
        if t < rollout_years:
            return 1 / rollout_years
        return 0

    # --- Intervention cascade inputs ---
    coverage_testing = float(inputs.get("coverage_testing", 0.5))
    coverage_treatment = float(inputs.get("coverage_treatment", 0.7))
    ltbi_prev = float(inputs.get("ltbi_prev", 0.25))

    # --- Compute full_effect ---
    full_effect, cascade = compute_full_effect(
        testing, treatment, ltbi_prev, coverage_testing, coverage_treatment
    )

    # --- Time series ---
    years = np.arange(0, T + 1)

    intervention_multiplier = np.array(
        [1 - rollout_factor(t) * (1-full_effect) for t in years]
    )

    # Baseline curves
    baseline_incidence = base_incidence * np.exp(-0.0 * years)
    baseline_cases = N * baseline_incidence / 100_000

    # Intervention curves
    intervention_incidence = (
        base_incidence * intervention_multiplier * np.exp(-0.0 * years)
    )
    intervention_cases = N * intervention_incidence / 100_000

    df = pd.DataFrame(
        {
            "Year": np.tile(years, 2),
            "Scenario": np.repeat(["Baseline", "Intervention"], len(years)),
            "Incidence_per_100k": np.concatenate(
                [baseline_incidence, intervention_incidence]
            ),
            "Cases": np.concatenate([baseline_cases, intervention_cases]),
        }
    )

    summary = {
        "Peak_incidence_baseline": float(np.max(baseline_incidence)),
        "Peak_incidence_intervention": float(np.max(intervention_incidence)),
        "Total_cases_baseline": float(np.sum(baseline_cases)),
        "Total_cases_intervention": float(np.sum(intervention_cases)),
        "Cascade": cascade,
    }

    return df, summary
