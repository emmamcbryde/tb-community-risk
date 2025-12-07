import numpy as np
import pandas as pd


def simulate_dynamic(params, years):
    """
    Dynamic TB model with:
    - age-structured seeding,
    - fast & slow latent compartments,
    - risk factor multipliers,
    - force-of-infection transmission.
    """

    # -------- Extract parameters passed from UI ----------
    beta = params["beta"]
    smoker_pct = params["smoker_pct"]
    diabetes_pct = params["diabetes_pct"]
    renal_pct = params["renal_pct"]
    immune_pct = params["immune_pct"]

    ltbi_ever = params["ltbi_ever"]        # dict age → prob ever infected
    ltbi_recent = params["ltbi_recent"]    # dict age → prob infected in last 5 yrs
    age_counts = params["age_counts"]      # dict age → population count

    sigma_fast = 0.01      # annual fast progression rate
    sigma_slow = 0.001     # annual slow progression rate
    omega = 1/5            # fast → slow transition rate
    gamma = 1/2            # 2-year mean infectious period

    # ----------- Combined risk multiplier ---------------
    risk_multiplier = (
        1
        + 0.01 * smoker_pct
        + 0.02 * diabetes_pct
        + 0.03 * renal_pct
        + 0.05 * immune_pct
    )

    sigma_fast_eff = sigma_fast * risk_multiplier
    sigma_slow_eff = sigma_slow * risk_multiplier

    # ----------- Seed initial compartments --------------
    S0 = 0
    L_fast0 = 0
    L_slow0 = 0
    I0 = 0
    R0 = 0

    for a, pop in age_counts.items():
        p_ever = ltbi_ever.get(a, 0)
        p_recent = ltbi_recent.get(a, 0)

        S0 += pop * (1 - p_ever)
        L_fast0 += pop * p_recent
        L_slow0 += pop * (p_ever - p_recent)

    # ----------- Prepare simulation arrays --------------
    t = np.arange(0, years + 1)
    S = np.zeros(len(t))
    L_fast = np.zeros(len(t))
    L_slow = np.zeros(len(t))
    I = np.zeros(len(t))
    R = np.zeros(len(t))

    S[0] = S0
    L_fast[0] = L_fast0
    L_slow[0] = L_slow0
    I[0] = I0
    R[0] = R0

    N = S0 + L_fast0 + L_slow0 + I0 + R0

    # ----------- Euler integration loop -----------------
    for i in range(1, len(t)):

        lambda_t = beta * I[i-1] / N

        dS = -lambda_t * S[i-1]
        dLf = lambda_t * S[i-1] - (sigma_fast_eff + omega) * L_fast[i-1]
        dLs = omega * L_fast[i-1] - sigma_slow_eff * L_slow[i-1]
        dI = sigma_fast_eff * L_fast[i-1] + sigma_slow_eff * L_slow[i-1] - gamma * I[i-1]
        dR = gamma * I[i-1]

        S[i] = S[i-1] + dS
        L_fast[i] = L_fast[i-1] + dLf
        L_slow[i] = L_slow[i-1] + dLs
        I[i] = I[i-1] + dI
        R[i] = R[i-1] + dR

    return {
        "time": t,
        "S": S,
        "L_fast": L_fast,
        "L_slow": L_slow,
        "I": I,
        "R": R,
        "incidence": I
    }

