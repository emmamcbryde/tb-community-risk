import numpy as np
import pandas as pd
from engine.params import extract_core_parameters
from engine.intervention import compute_full_effect


def simulate_community(inputs, file_path="data/parameters.xlsx"):

    # Inputs
    N = float(inputs.get("population", 10000))
    T = int(inputs.get("time_horizon", 20))
    testing = inputs.get("testing", "None")
    treatment = inputs.get("treatment", "None")
    rollout_years = float(inputs.get("rollout_years", 3.0))

    # Risk factors
    smoker = inputs.get("smoker_pct", 0) / 100.0
    diabetes = inputs.get("diabetes_pct", 0) / 100.0
    renal = inputs.get("renal_pct", 0) / 100.0
    immune = inputs.get("immune_pct", 0) / 100.0
    recent_inf = inputs.get("recent_infection", 0) / 100.0

    # Multipliers
    rr_smoker = 1.5
    rr_diabetes = 1.8
    rr_renal = 2.0
    rr_immune = 4.0
    rr_recent_inf = 10.0

    # Combined baseline risk
    baseline_risk = (
        1.0
        + smoker * (rr_smoker - 1.0)
        + diabetes * (rr_diabetes - 1.0)
        + renal * (rr_renal - 1.0)
        + immune * (rr_immune - 1.0)
        + recent_inf * (rr_recent_inf - 1.0)
    )

    # --- rollout function (user-defined years) ---
    def rollout_factor(t):
        if t <= 0:
            return 0.0
        if t < rollout_years:
            return t / rollout_years
        return 1.0

    # --- Intervention cascade inputs ---
    coverage_testing = float(inputs.get("coverage_testing", 0.5))
    coverage_treatment = float(inputs.get("coverage_treatment", 0.7))
    ltbi_prev_global = float(inputs.get("ltbi_prev", 0.25))  # global slider for static model

    # --- Compute full_effect for static model ---
    full_effect, cascade = compute_full_effect(
        testing, treatment, ltbi_prev_global, coverage_testing, coverage_treatment
    )

    years = np.arange(0, T + 1)

    # --------------------------------------------------------------------
    # AGE-STRUCTURED MECHANISTIC BASELINE INCIDENCE USING LTBI SPLIT
    # --------------------------------------------------------------------
    age_counts = inputs.get("age_counts", {})
    ltbi_recent_by_age = inputs.get("ltbi_recent_by_age", {})
    ltbi_remote_by_age = inputs.get("ltbi_remote_by_age", {})

    fast_prog = 0.01    # recent → TB
    slow_prog = 0.001   # remote → TB

    baseline_cases_age = []
    for a, N_a in age_counts.items():
        Lr = ltbi_recent_by_age.get(a, 0.0)
        Ls = ltbi_remote_by_age.get(a, 0.0)
        cases_a = float(N_a) * (Lr * fast_prog + Ls * slow_prog)
        baseline_cases_age.append(cases_a)

    raw_baseline_cases = np.sum(baseline_cases_age)
    baseline_total_cases0 = raw_baseline_cases * baseline_risk
    baseline_incidence0 = baseline_total_cases0 / N * 1e5

    # Keep baseline incidence flat over time (you can add decline if wanted)
    baseline_incidence = np.full_like(years, baseline_incidence0, dtype=float)
    baseline_cases = np.full_like(years, baseline_total_cases0, dtype=float)

    # --------------------------------------------------------------------
    # INTERVENTION MULTIPLIER (ROLL-IN EFFECT)
    # --------------------------------------------------------------------
    intervention_multiplier = np.array(
        [1.0 - rollout_factor(t) * (1.0 - full_effect) for t in years]
    )

    intervention_incidence = baseline_incidence * intervention_multiplier
    intervention_cases = baseline_cases * intervention_multiplier

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
