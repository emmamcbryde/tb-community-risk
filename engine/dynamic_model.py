# engine/dynamic_model.py

import numpy as np
import pandas as pd
from engine.params import extract_core_parameters
from engine.intervention import compute_full_effect


def simulate_dynamic_ltbi(
    age_counts,
    ltbi_prev,
    inputs,
    file_path="data/parameters.xlsx",
):
    """
    Dynamic model with two latency compartments:
      - L_fast : recent infection (<5 years), fast progression
      - L_slow : remote infection, slow progression
    """

    # Approximate annual risk of infection (ARI) from incidence (Stýblo-like rule)
    # ARI is per-person per-year probability of infection
    # 1. Prefer ARI derived from the historical incidence backcast (if provided)
    base_ari = inputs.get("base_ari", None)

    # 2. If not provided, fall back to a simple Stybló-like mapping from user_incidence
    if base_ari is None:
        user_incidence = inputs.get("user_incidence", 30)
        # Rough rule: ARI (fraction) ≈ incidence_per_100k / 5000
        base_ari = user_incidence / 5000.0

    # Safety cap: don't let ARI exceed, say, 10% per year
    base_ari = min(base_ari, 0.10)

    params = extract_core_parameters(file_path)
    mort = params.get("mortality", {})
    tbmort = params.get("tb_mortality", {})
    react = params.get("reactivation", {})
    T = inputs.get("time_horizon", 20)
    # If RRates not loaded, use defaults:
    if not react:
        react = {a: 0.001 for a in range(0, 101)}  # slow progression by default

    # If mortality not loaded, assume flat mortality curve (to avoid zero-pop)
    if not mort:
        mort = {a: 0.01 for a in range(0, 101)}  # 1% annual mortality

    # If TB mortality missing:
    if not tbmort:
        tbmort = {a: 0.20 for a in range(0, 101)}  # 20% annual death in TB

    ages = sorted(age_counts.keys())
    A = len(ages)

    # Compartments: age x time
    S = np.zeros((A, T + 1))
    L_fast = np.zeros((A, T + 1))
    L_slow = np.zeros((A, T + 1))
    I = np.zeros((A, T + 1))
    R = np.zeros((A, T + 1))

    # Initial population by age
    N0 = np.array([age_counts[a] for a in ages])

    # LTBI prevalence at t=0
    LTBI0 = np.array([ltbi_prev.get(a, 0.0) for a in ages])

    # At time 0 → assume 95% LTBI is "slow latency"
    S[:, 0] = (1 - LTBI0) * N0
    L_fast[:, 0] = LTBI0 * N0 * 0.05
    L_slow[:, 0] = LTBI0 * N0 * 0.95
    I[:, 0] = 0.0
    R[:, 0] = 0.0

    # Fast/slow progression hazards
    r_fast = 0.01  # recent infection progression
    r_slow = 0.001  # remote infection progression

    # Recent infection decays into slow infection over ~5 years
    tau_fast_to_slow = 1 / 5.0

    # Treatment rollout over 3 years
    def ltbi_treatment_coverage(t):
        if t < 0:
            return 0
        if t < 3:
            return (t + 1) / 3
        return 1.0

        # Compute intervention effect once at model start

    testing = inputs.get("testing", "None")
    treatment = inputs.get("treatment", "None")
    coverage_testing = inputs.get("coverage_testing", 0.0)
    coverage_treatment = inputs.get("coverage_treatment", 0.0)

    # Use population-average LTBI prevalence (or replace with age-specific later)
    avg_ltbi = float(np.mean([ltbi_prev.get(a, 0.0) for a in ages]))

    full_effect, cascade = compute_full_effect(
        testing,
        treatment,
        ltbi_prev=avg_ltbi,
        coverage_testing=coverage_testing,
        coverage_treatment=coverage_treatment,
    )
    treat_eff = 1 - full_effect  # so “effect” means cure probability added into model

    total_incidence = np.zeros(T + 1)

    for t in range(T):
        for i, a in enumerate(ages):

            # death hazards
            mu = mort.get(a, 0.0)
            mu_tb = tbmort.get(a, 0.05)
            # ---- New infections (force of infection / ARI) ----
            # Simple approximation: ARI declines slowly over time, consistent with static model
            ari_t = base_ari * np.exp(-0.0 * t)

            # New infections: susceptible → recent LTBI (all into L_fast)
            new_inf = S[i, t] * ari_t

            # ---- Progression to active TB ----
            new_TB_fast = L_fast[i, t] * r_fast
            new_TB_slow = L_slow[i, t] * r_slow
            new_TB = new_TB_fast + new_TB_slow

            # ---- Natural ageing of latency ----
            move_fast_to_slow = L_fast[i, t] * tau_fast_to_slow

            # ---- LTBI treatment ----
            covL = ltbi_treatment_coverage(t)
            treated_fast = L_fast[i, t] * covL
            treated_slow = L_slow[i, t] * covL

            # Successful cures
            cure_fast = treated_fast * treat_eff
            cure_slow = treated_slow * treat_eff

            # ---- Deaths ----
            dS = S[i, t] * mu
            dLF = L_fast[i, t] * mu
            dLS = L_slow[i, t] * mu
            dI = I[i, t] * mu_tb
            dR = R[i, t] * mu

            # ---- Updates ----
            # ---- Updates ----
            # Susceptibles lose new_inf, gain successfully treated LTBI
            S[i, t + 1] = S[i, t] - dS - new_inf + cure_fast + cure_slow

            # Recent LTBI gains new infections
            L_fast[i, t + 1] = (
                L_fast[i, t]
                + new_inf
                - new_TB_fast
                - move_fast_to_slow
                - treated_fast
                - dLF
            )

            L_slow[i, t + 1] = (
                L_slow[i, t]
                - new_TB_slow
                - treated_slow
                - dLS
                + move_fast_to_slow  # incoming from fast
            )

            # Active TB
            cureTB = I[i, t] * 0.7
            I[i, t + 1] = I[i, t] + new_TB - dI - cureTB

            R[i, t + 1] = R[i, t] + cureTB - dR

            total_incidence[t] += new_TB

    total_incidence[T] = total_incidence[T - 1]
    total_pop0 = np.sum(N0)
    incidence_per_100k = total_incidence / total_pop0 * 1e5

    years = np.arange(0, T + 1)
    df = pd.DataFrame(
        {
            "Year": years,
            "Scenario": ["Dynamic_LTBI_split"] * (T + 1),
            "Incidence_per_100k": incidence_per_100k,
            "Cases": total_incidence,
        }
    )

    summary = {
        "Peak_incidence": float(np.max(incidence_per_100k)),
        "Total_cases": float(np.sum(total_incidence)),
    }

    return df, summary
