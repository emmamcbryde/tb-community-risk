import numpy as np
import pandas as pd


def simulate_dynamic(params, years, intervention=True):
    """
    Dynamic TB model with:
    - age-structured seeding,
    - fast & slow latent compartments,
    - risk factor multipliers,
    - force-of-infection transmission,
    - LTBI test & treat as a pulse (fraction of population per year),
      returning L_fast and L_slow to S.
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

    initial_inc = params.get("initial_incidence_per_100k", 0.0)
    pre_det_months = params.get("pre_det_months", 12.0)
    delta_pre = params.get("delta_pre", 0.5)
    delta_post = params.get("delta_post", delta_pre)

    ltbi_coverage = params.get("ltbi_coverage", 0.0)   # total fraction of population ever test+treated
    rollout_years = params.get("rollout_years", 0)

    treatment_method = params.get("treatment_method", "None")

    # Natural history (base) parameters
    sigma_fast = 0.01      # annual fast progression rate
    sigma_slow = 0.001     # annual slow progression rate
    omega = 1/5            # fast → slow transition rate

    # ---------- Risk multipliers on progression ----------
    risk_multiplier = (
        1
        + 0.01 * smoker_pct
        + 0.02 * diabetes_pct
        + 0.03 * renal_pct
        + 0.05 * immune_pct
    )

    sigma_fast_eff = sigma_fast * risk_multiplier
    sigma_slow_eff = sigma_slow * risk_multiplier

    # ---------- Seed initial compartments ----------------
    S0 = 0.0
    L_fast0 = 0.0
    L_slow0 = 0.0

    for a, pop in age_counts.items():
        p_ever = ltbi_ever.get(a, 0.0)
        p_recent = ltbi_recent.get(a, 0.0)

        S0      += pop * (1.0 - p_ever)
        L_fast0 += pop * p_recent
        L_slow0 += pop * (p_ever - p_recent)

    # Total population (from age structure)
    N_total = sum(age_counts.values())

    # Prevalence = Incidence × Duration (years)
    duration_years = pre_det_months / 12.0
    I0 = N_total * (initial_inc / 100000.0) * duration_years

    # Subtract active TB from susceptible pool once
    S0 -= I0
    R0 = 0.0

    # Recompute total
    N = S0 + L_fast0 + L_slow0 + I0 + R0

    # --------- Prepare simulation arrays -----------------
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

    # --------- Intervention: LTBI test & treat pulse -----

    # LTBI regimen efficacies (placeholder – could be moved to params later)
    regimen_efficacy = {
        "None": 0.0,
        "1HP": 0.90,
        "3HP": 0.90,
        "4R": 0.80,
        "6H": 0.70,
        "9H": 0.75,
    }

    eff = regimen_efficacy.get(treatment_method, 0.0)

    # Effective total fraction of population that will ever be cured of LTBI
    # (user coverage × regimen efficacy)
    coverage_total = ltbi_coverage * eff

    if rollout_years > 0:
        annual_fraction_treated = coverage_total / rollout_years
    else:
        annual_fraction_treated = 0.0

    def tau_t(i_year):
        """
        Pulse LTBI treatment: fixed fraction per year during rollout, zero after.
        i_year is the integer year index.
        """
        if not intervention:
            return 0.0
        if rollout_years <= 0:
            return 0.0
        if i_year <= rollout_years:
            return annual_fraction_treated
        return 0.0

    # ---------- Diagnosis improvement (gamma(t)) ----------
    gamma_pre = delta_pre
    gamma_post = delta_post

    def gamma_t(i_year):
        """
        Time-varying diagnosis rate.
        Simple ramp from gamma_pre to gamma_post over rollout_years if intervention,
        otherwise stays at gamma_pre.
        """
        if not intervention:
            return gamma_pre
        if rollout_years <= 0:
            return gamma_post
        frac = min(1.0, i_year / rollout_years)
        return gamma_pre + (gamma_post - gamma_pre) * frac

    # ------------- Euler integration loop ----------------
    for i in range(1, len(t)):

        if N <= 0:
            # Guard against numerical blow-up
            S[i] = S[i-1]
            L_fast[i] = L_fast[i-1]
            L_slow[i] = L_slow[i-1]
            I[i] = I[i-1]
            R[i] = R[i-1]
            continue

        lambda_t = beta * I[i-1] / N

        tau_now = tau_t(i)
        gamma_now = gamma_t(i)

        # --- ODEs with LTBI treatment returning to S ---

        # Susceptible: infection - plus cures from LTBI
        dS = (
            -lambda_t * S[i-1]
            + tau_now * (L_fast[i-1] + L_slow[i-1])
        )

        # Fast latent: infection in – progression – transition to slow – LTBI cure
        dLf = (
            lambda_t * S[i-1]
            - (sigma_fast_eff + omega) * L_fast[i-1]
            - tau_now * L_fast[i-1]
        )

        # Slow latent: from fast – progression – LTBI cure
        dLs = (
            omega * L_fast[i-1]
            - sigma_slow_eff * L_slow[i-1]
            - tau_now * L_slow[i-1]
        )

        # Active TB: progression from both latent compartments – recovery
        dI = (
            sigma_fast_eff * L_fast[i-1]
            + sigma_slow_eff * L_slow[i-1]
            - gamma_now * I[i-1]
        )

        # Recovered: from active TB treatment
        dR = gamma_now * I[i-1]

        # Update states
        S[i]      = S[i-1] + dS
        L_fast[i] = L_fast[i-1] + dLf
        L_slow[i] = L_slow[i-1] + dLs
        I[i]      = I[i-1] + dI
        R[i]      = R[i-1] + dR

        # Recompute N (should remain approx constant)
        N = S[i] + L_fast[i] + L_slow[i] + I[i] + R[i]

    return {
        "time": t,
        "S": S,
        "L_fast": L_fast,
        "L_slow": L_slow,
        "I": I,
        "R": R,
        "incidence": I,   # you’re currently using I as "incidence proxy" in UI
    }
