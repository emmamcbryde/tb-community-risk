# engine/dynamic_model.py

import numpy as np
import pandas as pd
from engine.params import extract_core_parameters

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

    params = extract_core_parameters(file_path)
    mort = params.get("mortality", {})
    tbmort = params.get("tb_mortality", {})
    react = params.get("reactivation", {})

    T = inputs.get("time_horizon", 20)

    ages = sorted(age_counts.keys())
    A = len(ages)

    # Compartments: age x time
    S       = np.zeros((A, T+1))
    L_fast  = np.zeros((A, T+1))
    L_slow  = np.zeros((A, T+1))
    I       = np.zeros((A, T+1))
    R       = np.zeros((A, T+1))

    # Initial population by age
    N0 = np.array([age_counts[a] for a in ages])

    # LTBI prevalence at t=0
    LTBI0 = np.array([ltbi_prev.get(a, 0.0) for a in ages])

    # At time 0 → assume ALL LTBI is "slow latency"
    S[:, 0]       = (1 - LTBI0) * N0
    L_fast[:, 0]  = 0.0
    L_slow[:, 0]  = LTBI0 * N0
    I[:, 0]       = 0.0
    R[:, 0]       = 0.0

    # Fast/slow progression hazards
    r_fast = 0.01   # recent infection progression
    r_slow = 0.001  # remote infection progression

    # Recent infection decays into slow infection over ~5 years
    tau_fast_to_slow = 1/5.0

    # Treatment rollout over 3 years
    def ltbi_treatment_coverage(t):
        if t < 0:
            return 0
        if t < 3:
            return (t+1)/3
        return 1.0

    treat_eff = 0.8

    total_incidence = np.zeros(T+1)

    for t in range(T):
        for i, a in enumerate(ages):

            # death hazards
            mu = mort.get(a, 0.0)
            mu_tb = tbmort.get(a, 0.05)

            # ---- Progression to active TB ----
            new_TB_fast = L_fast[i,t] * r_fast
            new_TB_slow = L_slow[i,t] * r_slow
            new_TB = new_TB_fast + new_TB_slow

            # ---- Natural ageing of latency ----
            move_fast_to_slow = L_fast[i,t] * tau_fast_to_slow

            # ---- LTBI treatment ----
            covL = ltbi_treatment_coverage(t)
            treated_fast = L_fast[i,t] * covL
            treated_slow = L_slow[i,t] * covL

            # Successful cures
            cure_fast = treated_fast * treat_eff
            cure_slow = treated_slow * treat_eff

            # ---- Deaths ----
            dS = S[i,t] * mu
            dLF = L_fast[i,t] * mu
            dLS = L_slow[i,t] * mu
            dI = I[i,t] * mu_tb
            dR = R[i,t] * mu

            # ---- Updates ----
            S[i,t+1] = S[i,t] - dS + cure_fast + cure_slow

            L_fast[i,t+1] = (
                L_fast[i,t]
                - new_TB_fast
                - move_fast_to_slow
                - treated_fast
                - dLF
            )

            L_slow[i,t+1] = (
                L_slow[i,t]
                - new_TB_slow
                - treated_slow
                - dLS
                + move_fast_to_slow  # incoming from fast
            )

            # Active TB
            cureTB = I[i,t] * 0.7
            I[i,t+1] = (
                I[i,t]
                + new_TB
                - dI
                - cureTB
            )

            R[i,t+1] = R[i,t] + cureTB - dR

            total_incidence[t] += new_TB

    total_incidence[T] = total_incidence[T-1]
    total_pop0 = np.sum(N0)
    incidence_per_100k = total_incidence / total_pop0 * 1e5

    years = np.arange(0, T+1)
    df = pd.DataFrame({
        "Year": years,
        "Scenario": ["Dynamic_LTBI_split"] * (T+1),
        "Incidence_per_100k": incidence_per_100k,
        "Cases": total_incidence
    })

    summary = {
        "Peak_incidence": float(np.max(incidence_per_100k)),
        "Total_cases": float(np.sum(total_incidence))
    }

    return df, summary
