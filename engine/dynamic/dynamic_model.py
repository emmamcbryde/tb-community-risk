import numpy as np
import pandas as pd


def simulate_dynamic(params, years, intervention=True):
    """
    Dynamic TB model with:
    - Age-structured seeding from LTBI (ltbi_ever, ltbi_recent, age_counts)
    - Fast & slow latent compartments
    - Risk-factor multipliers on progression
    - Force-of-infection transmission
    - Pulse-style LTBI test & treat (fraction per year over rollout period),
      removing L_fast & L_slow to S
    - Diagnosis improvement (gamma) over rollout years
    - Fine time-step integration (dt = 0.1 years) with ANNUALISED output.
    """
    required_keys = [
    "smoker_pct", "alcohol_pct", "diabetes_pct",
    "renal_pct", "HIV_treated_pct", "HIV_untreated_pct"]

    for k in required_keys:
        if k not in params:
            raise ValueError(f"Missing parameter: {k}")


    # -------- Extract parameters --------
    beta = params["beta"]
    smoker_pct   = params["smoker_pct"]
    diabetes_pct = params["diabetes_pct"]
    renal_pct    = params["renal_pct"]
    HIV_untreated_pct = params["HIV_untreated_pct"]
    HIV_treated_pct   = params["HIV_treated_pct"] 
    alcohol_pct   = params["alcohol_pct"]
    ltbi_ever   = params["ltbi_ever"]        # dict age → Prob(ever infected)
    ltbi_recent = params["ltbi_recent"]      # dict age → Prob(recent infection)
    age_counts  = params["age_counts"]       # dict age → population

    initial_inc = params.get("initial_incidence_per_100k", 0.0)
    pre_det_months = params.get("pre_det_months", 12.0)
    delta_pre  = params.get("delta_pre", 0.5)
    delta_post = params.get("delta_post", delta_pre)

    ltbi_coverage    = params.get("ltbi_coverage", 0.0)
    rollout_years    = params.get("rollout_years", 0)
    treatment_method = params.get("treatment_method", "None")

    # -------- Natural history constants --------
    sigma_fast = 0.01   # per year
    sigma_slow = 0.001
    omega      = 1.0/5.0   # fast → slow per year


    # --------------------------------------------------
    # Risk factor relative risks (progression LTBI -> TB)
    # --------------------------------------------------
    RR = {
        "smoker": 1.5,
        "alcohol": 2.0,
        "diabetes": 3.0,
        "renal": 2.5,
        "HIV_treated": 4.0,
        "HIV_untreated": 10.0,
    }

    # Population prevalences (fractions)
    p = {
        "smoker": smoker_pct / 100.0,
        "alcohol": alcohol_pct / 100.0,
        "diabetes": diabetes_pct / 100.0,
        "renal": renal_pct / 100.0,
        "HIV_treated": HIV_treated_pct / 100.0,
        "HIV_untreated": HIV_untreated_pct / 100.0,
    }

    # Combine multiplicatively
    risk_multiplier = 1.0
    for k in RR:
        risk_multiplier *= (1 - p[k]) + p[k] * RR[k]


    sigma_fast_eff = sigma_fast * risk_multiplier
    sigma_slow_eff = sigma_slow * risk_multiplier

    # -------- LTBI regimen efficacy --------
    regimen_efficacy = {
        "None": 0.0,
        "1HP": 0.90,
        "3HP": 0.90,
        "4R": 0.80,
        "6H": 0.70,
        "9H": 0.75,
    }
    eff = regimen_efficacy.get(treatment_method, 0.0)

    # Effective total coverage of LTBI (over rollout period)
    coverage_total = ltbi_coverage * eff

    if rollout_years > 0:
        annual_fraction_treated = coverage_total / rollout_years
    else:
        annual_fraction_treated = 0.0

    # -------- Seeding compartments from LTBI --------
    S0 = 0.0
    L_fast0 = 0.0
    L_slow0 = 0.0

    for a, pop in age_counts.items():
        p_ever   = ltbi_ever.get(a,   0.0)
        p_recent = ltbi_recent.get(a, 0.0)

        S0      += pop * (1.0 - p_ever)
        L_fast0 += pop * p_recent
        L_slow0 += pop * (p_ever - p_recent)

    N_total = sum(age_counts.values())

    # Seed active TB prevalence from incidence × duration
    duration_years = pre_det_months / 12.0
    I0 = N_total * (initial_inc / 100000.0) * duration_years

    S0 -= I0
    R0 = 0.0

    N = S0 + L_fast0 + L_slow0 + I0 + R0

    # Guard
    if N <= 0:
        N = 1.0

    # -------- Time grid (dt = 0.1 years) --------
    dt = 0.1
    n_steps = int(years / dt)
    t = np.linspace(0.0, years, n_steps + 1)

    S      = np.zeros(n_steps + 1)
    L_fast = np.zeros(n_steps + 1)
    L_slow = np.zeros(n_steps + 1)
    I      = np.zeros(n_steps + 1)
    R      = np.zeros(n_steps + 1)

    S[0]      = S0
    L_fast[0] = L_fast0
    L_slow[0] = L_slow0
    I[0]      = I0
    R[0]      = R0

    # -------- Tau(t) and gamma(t) functions --------
    def tau_at_time(time):
        if not intervention:
            return 0.0
        if rollout_years <= 0:
            return 0.0
        if time <= rollout_years:
            return annual_fraction_treated
        return 0.0

    def gamma_at_time(time):
        if not intervention:
            return delta_pre
        if rollout_years <= 0:
            return delta_post
        frac = min(1.0, time / rollout_years)
        return delta_pre + (delta_post - delta_pre) * frac

    # -------- Integrate with Euler (dt = 0.1) --------
    for i in range(1, n_steps + 1):
        time = t[i]

        if N <= 0:
            S[i] = S[i-1]
            L_fast[i] = L_fast[i-1]
            L_slow[i] = L_slow[i-1]
            I[i] = I[i-1]
            R[i] = R[i-1]
            continue

        lambda_t = beta * I[i-1] / N

        tau_now   = tau_at_time(time)
        gamma_now = gamma_at_time(time)

        # d/dt terms
        dS_dt = (
            -lambda_t * S[i-1]
            + tau_now * (L_fast[i-1] + L_slow[i-1])
        )

        dLf_dt = (
            lambda_t * S[i-1]
            - (sigma_fast_eff + omega + tau_now) * L_fast[i-1]
        )

        dLs_dt = (
            omega * L_fast[i-1]
            - (sigma_slow_eff + tau_now) * L_slow[i-1]
        )

        dI_dt = (
            sigma_fast_eff * L_fast[i-1]
            + sigma_slow_eff * L_slow[i-1]
            - gamma_now * I[i-1]
        )

        dR_dt = gamma_now * I[i-1]

        # Euler step
        S[i]      = S[i-1]      + dS_dt  * dt
        L_fast[i] = L_fast[i-1] + dLf_dt * dt
        L_slow[i] = L_slow[i-1] + dLs_dt * dt
        I[i]      = I[i-1]      + dI_dt  * dt
        R[i]      = R[i-1]      + dR_dt  * dt

        N = S[i] + L_fast[i] + L_slow[i] + I[i] + R[i]

    # -------- Annualise output (at exact years) --------
    annual_years = np.arange(0, years + 1)
    annual_I = []

    for y in annual_years:
        idx = int(round(y / dt))
        idx = min(idx, n_steps)
        annual_I.append(I[idx])

    annual_I = np.array(annual_I)

    return {
        "time": annual_years,
        "S": S,
        "L_fast": L_fast,
        "L_slow": L_slow,
        "I": I,
        "R": R,
        # UI treats this as "incidence proxy"; now smooth at annual points
        "incidence": annual_I,
    }
