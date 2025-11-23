# engine/dynamic_model.py

import numpy as np
import pandas as pd
from engine.params import extract_core_parameters
from engine.intervention import REGIMENS


def simulate_dynamic_ltbi(
    age_counts,
    ltbi_prev,
    inputs,
    file_path="data/parameters.xlsx",
):
    """
    Dynamic age-structured LTBI → TB model.
      S       : Susceptible
      L_fast  : Recent infection (<5 years), fast progression
      L_slow  : Remote infection, slow progression
      I       : Active TB
      R       : Recovered
    """

    # --- 1. Force of infection (ARI) ---
    # Prefer ARI derived from historical incidence backcast
    base_ari = inputs.get("base_ari", None)

    if base_ari is None:
        user_incidence = inputs.get("user_incidence", 30)
        # Stýblo-like approximation
        base_ari = user_incidence / 5000.0

    # Cap for safety
    base_ari = min(base_ari, 0.10)

    # --- 2. Load parameters ---
    params = extract_core_parameters(file_path)
    mort = params.get("mortality", {})
    tbmort = params.get("tb_mortality", {})
    react = params.get("reactivation", {})

    T = inputs.get("time_horizon", 20)

    if not react:
        react = {a: 0.001 for a in range(0, 101)}

    if not mort:
        mort = {a: 0.01 for a in range(0, 101)}

    if not tbmort:
        tbmort = {a: 0.20 for a in range(0, 101)}

    ages = sorted(age_counts.keys())
    A = len(ages)

    # --- 3. Regimen cure logic (dynamic model uses LTBI compartments directly) ---
    treatment = inputs.get("treatment", "None")
    coverage_treatment = inputs.get("coverage_treatment", 0.0)
    completion_override = inputs.get("completion_override", None)
    rollout_years = inputs.get("rollout_years", 3),
    reg = REGIMENS.get(treatment, REGIMENS["None"])
    efficacy = reg["efficacy"]

    if completion_override is not None:
        completion_used = completion_override
    else:
        completion_used = reg["completion"]

    # Cure probability per treated LTBI individual
    regimen_cure = efficacy * completion_used

    # --- 4. Initialise compartments ---
    S = np.zeros((A, T + 1))
    L_fast = np.zeros((A, T + 1))
    L_slow = np.zeros((A, T + 1))
    I = np.zeros((A, T + 1))
    R = np.zeros((A, T + 1))

    N0 = np.array([age_counts[a] for a in ages])

    LTBI0 = np.array([ltbi_prev.get(a, 0.0) for a in ages])

    # Split LTBI initial conditions (5% recent, 95% remote)
    S[:, 0] = (1 - LTBI0) * N0
    L_fast[:, 0] = LTBI0 * N0 * 0.05
    L_slow[:, 0] = LTBI0 * N0 * 0.95
    I[:, 0] = 0.0
    R[:, 0] = 0.0

    # Progression hazards
    r_fast = 0.01  # recent infection → active TB
    r_slow = 0.001  # remote infection → active TB
    tau_fast_to_slow = 1 / 5.0  # recent → remote latency transition

    # Treatment rollout: ramps to 100% over 3 years
    def ltbi_treatment_coverage(t, rollout_years):
        if t < 0:
            return 0.0
        if t < rollout_years:
            return (t + 1) / rollout_years
        return 1.0

    # --- 5. Simulation loop ---
    total_incidence = np.zeros(T + 1)

    for t in range(T):
        for i, a in enumerate(ages):

            # --- Mortality ---
            mu = mort.get(a, 0.0)
            mu_tb = tbmort.get(a, 0.05)

            # --- New infections (S → L_fast) ---
            ari_t = base_ari  # constant for now
            new_inf = S[i, t] * ari_t

            # --- Progression (LTBI → TB) ---
            new_TB_fast = L_fast[i, t] * r_fast
            new_TB_slow = L_slow[i, t] * r_slow
            new_TB = new_TB_fast + new_TB_slow

            # --- Recent LTBI → Remote LTBI ---
            move_fast_to_slow = L_fast[i, t] * tau_fast_to_slow

            # --- LTBI treatment ---
            cov_rollout = ltbi_treatment_coverage(t, rollout_years)
            covL = cov_rollout * coverage_treatment

            treated_fast = L_fast[i, t] * covL
            treated_slow = L_slow[i, t] * covL

            # Successful cures from treated LTBI
            cure_fast = treated_fast * regimen_cure
            cure_slow = treated_slow * regimen_cure

            # --- Deaths ---
            dS = S[i, t] * mu
            dLF = L_fast[i, t] * mu
            dLS = L_slow[i, t] * mu
            dI = I[i, t] * mu_tb
            dR = R[i, t] * mu

            # --- Update compartments ---
            S[i, t + 1] = S[i, t] - dS - new_inf + cure_fast + cure_slow

            L_fast[i, t + 1] = (
                L_fast[i, t]
                + new_inf
                - new_TB_fast
                - move_fast_to_slow
                - treated_fast
                - dLF
            )

            L_slow[i, t + 1] = (
                L_slow[i, t] - new_TB_slow - treated_slow - dLS + move_fast_to_slow
            )

            cureTB = I[i, t] * 0.7
            I[i, t + 1] = I[i, t] + new_TB - dI - cureTB

            R[i, t + 1] = R[i, t] + cureTB - dR

            total_incidence[t] += new_TB

    # Keep final point stable
    total_incidence[T] = total_incidence[T - 1]

    total_pop0 = np.sum(N0)
    incidence_per_100k = total_incidence / total_pop0 * 1e5

    # --- Output ---
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
