import numpy as np


def simulate_dynamic(params, years, intervention=True):
    """
    Dynamic TB model with:
    - Age-structured seeding from LTBI (ltbi_ever, ltbi_recent, age_counts)
    - Fast & slow latent compartments
    - Risk-factor multipliers on progression
    - Force-of-infection transmission
    - Pulse-style LTBI test & treat (fraction per year over rollout period), removing L_fast & L_slow to S
    - Diagnosis improvement (gamma) over rollout years
    - Fine time-step integration (dt = 0.1 years) with ANNUALISED outputs.

    Returns:
    - annual_incidence: NEW active TB cases per year (counts), length = years
    - annual_prevalence_I: I compartment at integer years (counts), length = years+1
    """

    # ---------------------------
    # Required parameters
    # ---------------------------
    required_keys = [
        "beta",
        "age_counts",
        "ltbi_ever",
        "ltbi_recent",
        "smoker_pct",
        "alcohol_pct",
        "diabetes_pct",
        "renal_pct",
        "HIV_treated_pct",
        "HIV_untreated_pct",
    ]
    for k in required_keys:
        if k not in params:
            raise ValueError(f"Missing parameter: {k}")

    # ---------------------------
    # Extract parameters
    # ---------------------------
    beta = float(params["beta"])
    beta = float(params.get("beta", 0.0))
    beta_series = params.get("beta_series", None)
    if beta_series is not None:
        beta_series = np.asarray(beta_series, dtype=float).reshape(-1)
        if beta_series.size == 0:
            beta_series = None

    smoker_pct = float(params["smoker_pct"])
    alcohol_pct = float(params["alcohol_pct"])
    diabetes_pct = float(params["diabetes_pct"])
    renal_pct = float(params["renal_pct"])
    HIV_treated_pct = float(params["HIV_treated_pct"])
    HIV_untreated_pct = float(params["HIV_untreated_pct"])

    ltbi_ever = params["ltbi_ever"]  # dict age -> P(ever infected)
    ltbi_recent = params["ltbi_recent"]  # dict age -> P(infected in last 5y)
    age_counts = params["age_counts"]  # dict age -> population count

    initial_inc = float(params.get("initial_incidence_per_100k", 0.0))

    pre_det_months = float(params.get("pre_det_months", 12.0))
    delta_pre = float(params.get("delta_pre", 12.0 / max(pre_det_months, 0.1)))
    delta_post = float(params.get("delta_post", delta_pre))

    ltbi_coverage = float(params.get("ltbi_coverage", 0.0))
    rollout_years = int(params.get("rollout_years", 0))
    treatment_method = params.get("treatment_method", "None")

    # ---------------------------
    # Natural history constants
    # ---------------------------
    sigma_fast = 0.01  # per year
    sigma_slow = 0.001  # per year
    omega = 1.0 / 5.0  # fast -> slow per year (rate 1/5)

    # ---------------------------
    # Risk factor RR multipliers (progression LTBI -> TB)
    # ---------------------------
    RR = {
        "smoker": 1.5,
        "alcohol": 2.0,
        "diabetes": 3.0,
        "renal": 2.5,
        "HIV_treated": 4.0,
        "HIV_untreated": 10.0,
    }

    p = {
        "smoker": smoker_pct / 100.0,
        "alcohol": alcohol_pct / 100.0,
        "diabetes": diabetes_pct / 100.0,
        "renal": renal_pct / 100.0,
        "HIV_treated": HIV_treated_pct / 100.0,
        "HIV_untreated": HIV_untreated_pct / 100.0,
    }

    risk_multiplier = 1.0
    for k in RR:
        risk_multiplier *= (1.0 - p[k]) + p[k] * RR[k]

    sigma_fast_eff = sigma_fast * risk_multiplier
    sigma_slow_eff = sigma_slow * risk_multiplier

    # ---------------------------
    # LTBI regimen efficacy
    # ---------------------------
    regimen_efficacy = {
        "None": 0.0,
        "1HP": 0.90,
        "3HP": 0.90,
        "4R": 0.80,
        "6H": 0.70,
        "9H": 0.75,
    }
    eff = float(regimen_efficacy.get(treatment_method, 0.0))

    # Effective total coverage of LTBI (over rollout period)
    coverage_total = ltbi_coverage * eff
    if rollout_years > 0:
        annual_fraction_treated = coverage_total / rollout_years
    else:
        annual_fraction_treated = 0.0

    # ---------------------------
    # Seed compartments from LTBI by age
    # ---------------------------
    S0 = 0.0
    L_fast0 = 0.0
    L_slow0 = 0.0

    for a, pop in age_counts.items():
        p_ever = float(ltbi_ever.get(a, 0.0))
        p_recent = float(ltbi_recent.get(a, 0.0))
        S0 += pop * (1.0 - p_ever)
        L_fast0 += pop * p_recent
        L_slow0 += pop * (p_ever - p_recent)

    N_total = float(sum(age_counts.values()))

    # Seed I0 using: prevalence ≈ incidence * duration (duration = pre_det_months)
    duration_years = pre_det_months / 12.0
    I0 = N_total * (initial_inc / 100000.0) * duration_years

    # ensure we don't go negative
    S0 = max(S0 - I0, 0.0)
    R0 = 0.0

    # ---------------------------
    # Time grid (dt = 0.1 years)
    # ---------------------------
    dt = 0.1
    n_steps = int(round(years / dt))
    t = np.linspace(0.0, years, n_steps + 1)

    S = np.zeros(n_steps + 1)
    L_fast = np.zeros(n_steps + 1)
    L_slow = np.zeros(n_steps + 1)
    I = np.zeros(n_steps + 1)
    R = np.zeros(n_steps + 1)

    S[0], L_fast[0], L_slow[0], I[0], R[0] = S0, L_fast0, L_slow0, I0, R0

    # Accumulate annual incidence as cases per year interval [y, y+1)
    annual_incidence = np.zeros(years)  # length = years (0..years-1)

    # ---------------------------
    # Intervention functions
    # ---------------------------
    def tau_at_time(time):
        if (not intervention) or rollout_years <= 0:
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

    # ---------------------------
    # Integrate (Euler, dt=0.1)
    # ---------------------------
    for i in range(1, n_steps + 1):
        time_prev = t[i - 1]

        N = S[i - 1] + L_fast[i - 1] + L_slow[i - 1] + I[i - 1] + R[i - 1]
        if N <= 0:
            S[i], L_fast[i], L_slow[i], I[i], R[i] = (
                S[i - 1],
                L_fast[i - 1],
                L_slow[i - 1],
                I[i - 1],
                R[i - 1],
            )
            continue
        beta_now = beta
        if beta_series is not None:
            idx = int(np.floor(t))  # 0..years-1
            idx = max(0, min(idx, len(beta_series) - 1))
            beta_now = float(beta_series[idx])

        # Use beta_now everywhere beta was used:
        # lambda_inf = beta_now * I / N

        lambda_t = beta_now * I[i - 1] / N
        tau_now = tau_at_time(time_prev)
        gamma_now = gamma_at_time(time_prev)

        # New TB cases (incidence flow) this instant
        new_cases_rate = sigma_fast_eff * L_fast[i - 1] + sigma_slow_eff * L_slow[i - 1]

        # Accumulate annual cases
        year_idx = int(np.floor(time_prev))
        if 0 <= year_idx < years:
            annual_incidence[year_idx] += new_cases_rate * dt

        # ODEs
        dS_dt = (-lambda_t * S[i - 1]) + tau_now * (L_fast[i - 1] + L_slow[i - 1])
        dLf_dt = (lambda_t * S[i - 1]) - (sigma_fast_eff + omega + tau_now) * L_fast[
            i - 1
        ]
        dLs_dt = (omega * L_fast[i - 1]) - (sigma_slow_eff + tau_now) * L_slow[i - 1]
        dI_dt = new_cases_rate - gamma_now * I[i - 1]
        dR_dt = gamma_now * I[i - 1]

        # Euler step
        S[i] = S[i - 1] + dS_dt * dt
        L_fast[i] = L_fast[i - 1] + dLf_dt * dt
        L_slow[i] = L_slow[i - 1] + dLs_dt * dt
        I[i] = I[i - 1] + dI_dt * dt
        R[i] = R[i - 1] + dR_dt * dt

        # numeric safety
        S[i] = max(S[i], 0.0)
        L_fast[i] = max(L_fast[i], 0.0)
        L_slow[i] = max(L_slow[i], 0.0)
        I[i] = max(I[i], 0.0)
        R[i] = max(R[i], 0.0)

    # ---------------------------
    # Annual prevalence snapshots at integer years
    # ---------------------------
    annual_years = np.arange(0, years + 1)
    annual_prevalence_I = np.zeros(years + 1)
    for y in annual_years:
        idx = int(round(y / dt))
        idx = min(idx, n_steps)
        annual_prevalence_I[y] = I[idx]

    return {
        "time": annual_years,
        "annual_incidence_time": np.arange(0, years),
        "annual_incidence": annual_incidence,  # counts per year
        "annual_prevalence_I": annual_prevalence_I,  # I at year boundaries
        "S": S,
        "L_fast": L_fast,
        "L_slow": L_slow,
        "I": I,
        "R": R,
    }
