import numpy as np


def simulate_dynamic(params, years, intervention=True):
    """
    Dynamic TB transmission/progression model with:
    - Susceptible (S)
    - Fast latent (L_fast) and slow latent (L_slow)
    - Active TB (I)
    - Recovered (R)

    Key features:
    - Age-structured seeding from LTBI (ltbi_ever, ltbi_recent, age_counts) OR
      optional overriding via params["initial_state"].
    - Fast & slow latent progression (sigma_fast, sigma_slow) with fast->slow transition (omega).
    - Risk-factor multipliers on progression.
    - Force-of-infection transmission (lambda = beta * I / N).
    - Pulse-style LTBI test & treat: during rollout, a fixed fraction/year moves from L compartments to S.
    - Diagnosis improvement: gamma transitions I -> R with linear ramp over rollout years.
    - Euler integration with dt = 0.1 years.
    - ANNUALISED outputs: annual incidence (counts) per year.

    Parameters expected in `params`:
      Required always:
        - age_counts (dict age->count)
        - smoker_pct, alcohol_pct, diabetes_pct, renal_pct, HIV_treated_pct, HIV_untreated_pct
        - beta (float) OR beta_series (array-like)

      Required if initial_state is NOT provided:
        - ltbi_ever (dict age->P ever infected)
        - ltbi_recent (dict age->P infected in last 5y)

      Optional:
        - initial_state: dict with keys S, L_fast, L_slow, I, R (counts at t=0)
        - initial_incidence_per_100k (float): used only when initial_state is absent
        - pre_det_months, delta_pre, delta_post
        - ltbi_coverage (0..1), rollout_years (int)
        - treatment_method (str)
    """

    # ---------------------------
    # Validate essentials
    # ---------------------------
    if "age_counts" not in params:
        raise ValueError("Missing parameter: age_counts")

    risk_keys = [
        "smoker_pct",
        "alcohol_pct",
        "diabetes_pct",
        "renal_pct",
        "HIV_treated_pct",
        "HIV_untreated_pct",
    ]
    for k in risk_keys:
        if k not in params:
            raise ValueError(f"Missing parameter: {k}")

    # beta: support scalar beta and optional beta_series (per-year)
    beta_series = params.get("beta_series", None)
    if beta_series is not None:
        beta_series = np.asarray(beta_series, dtype=float).reshape(-1)
        if beta_series.size == 0:
            beta_series = None

    if beta_series is None:
        if "beta" not in params:
            raise ValueError("Missing parameter: beta")
        beta = float(params["beta"])
    else:
        # allow missing scalar beta (use last value as fallback)
        beta_raw = params.get("beta", float(beta_series[-1]))
        beta_arr = np.asarray(beta_raw, dtype=float).reshape(-1)
        beta = float(beta_arr[-1])

    age_counts = params["age_counts"]

    # LTBI inputs (only required if we are NOT overriding initial conditions)
    initial_state = params.get("initial_state", None)
    if initial_state is None:
        for k in ("ltbi_ever", "ltbi_recent"):
            if k not in params:
                raise ValueError(
                    f"Missing parameter: {k} (required when initial_state is not provided)"
                )
        ltbi_ever = params["ltbi_ever"]
        ltbi_recent = params["ltbi_recent"]
    else:
        ltbi_ever = params.get("ltbi_ever", {})
        ltbi_recent = params.get("ltbi_recent", {})

    # ---------------------------
    # Extract parameters
    # ---------------------------
    smoker_pct = float(params["smoker_pct"])
    alcohol_pct = float(params["alcohol_pct"])
    diabetes_pct = float(params["diabetes_pct"])
    renal_pct = float(params["renal_pct"])
    HIV_treated_pct = float(params["HIV_treated_pct"])
    HIV_untreated_pct = float(params["HIV_untreated_pct"])

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
    omega = 1.0 / 5.0  # fast -> slow per year

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

    coverage_total = ltbi_coverage * eff
    annual_fraction_treated = (
        coverage_total / rollout_years if rollout_years > 0 else 0.0
    )

    # ---------------------------
    # Seed compartments
    # ---------------------------
    N_total = float(sum(float(v) for v in age_counts.values()))

    if initial_state is not None:
        # Overridden initial state at t=0 (counts)
        S0 = float(initial_state.get("S", initial_state.get("S0", 0.0)))
        L_fast0 = float(initial_state.get("L_fast", initial_state.get("L_fast0", 0.0)))
        L_slow0 = float(initial_state.get("L_slow", initial_state.get("L_slow0", 0.0)))
        I0 = float(initial_state.get("I", initial_state.get("I0", 0.0)))
        R0 = float(initial_state.get("R", initial_state.get("R0", 0.0)))

        total0 = S0 + L_fast0 + L_slow0 + I0 + R0
        if total0 <= 0:
            raise ValueError("initial_state has non-positive total population.")

        # Make totals match the user population (age_counts) to keep per-100k conversions consistent
        diff = N_total - total0
        S0 = S0 + diff

        # Numerical safety (diff should be tiny)
        S0 = max(S0, 0.0)

    else:
        # Seed from LTBI by age
        S0 = 0.0
        L_fast0 = 0.0
        L_slow0 = 0.0

        for a, pop in age_counts.items():
            pop = float(pop)
            p_ever = float(ltbi_ever.get(a, 0.0))
            p_recent = float(ltbi_recent.get(a, 0.0))
            S0 += pop * (1.0 - p_ever)
            L_fast0 += pop * p_recent
            L_slow0 += pop * (p_ever - p_recent)

        # Seed I0 using prevalence ≈ incidence * duration
        duration_years = pre_det_months / 12.0
        I0 = N_total * (initial_inc / 100000.0) * duration_years

        # subtract prevalent TB from susceptible to keep totals stable
        S0 = max(S0 - I0, 0.0)
        R0 = 0.0

    # ---------------------------
    # Time grid (dt = 0.1 years)
    # ---------------------------
    years = int(years)
    if years <= 0:
        raise ValueError("years must be a positive integer.")

    dt = 0.1
    n_steps = int(round(years / dt))
    t = np.linspace(0.0, float(years), n_steps + 1)

    S = np.zeros(n_steps + 1)
    L_fast = np.zeros(n_steps + 1)
    L_slow = np.zeros(n_steps + 1)
    I = np.zeros(n_steps + 1)
    R = np.zeros(n_steps + 1)

    S[0], L_fast[0], L_slow[0], I[0], R[0] = S0, L_fast0, L_slow0, I0, R0

    annual_incidence = np.zeros(
        years, dtype=float
    )  # new TB cases per year interval [y, y+1)

    # ---------------------------
    # Intervention functions
    # ---------------------------
    def tau_at_time(time):
        if (not intervention) or rollout_years <= 0:
            return 0.0
        if time < float(rollout_years):
            return annual_fraction_treated
        return 0.0

    def gamma_at_time(time):
        if not intervention:
            return delta_pre
        if rollout_years <= 0:
            return delta_post
        frac = min(1.0, time / float(rollout_years))
        return delta_pre + (delta_post - delta_pre) * frac

    # ---------------------------
    # Integrate (Euler, dt=0.1)
    # ---------------------------
    for i in range(1, n_steps + 1):
        time_prev = float(t[i - 1])

        N = S[i - 1] + L_fast[i - 1] + L_slow[i - 1] + I[i - 1] + R[i - 1]
        if N <= 0:
            # carry forward
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
            beta_idx = int(np.floor(time_prev))  # 0..years-1
            beta_idx = max(0, min(beta_idx, beta_series.size - 1))
            beta_now = float(beta_series[beta_idx])

        lambda_t = beta_now * I[i - 1] / N
        tau_now = tau_at_time(time_prev)
        gamma_now = gamma_at_time(time_prev)

        # New TB cases (incidence flow) at this instant
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
    annual_years = np.arange(0, years + 1, dtype=float)
    annual_prevalence_I = np.zeros(years + 1, dtype=float)
    for y in range(0, years + 1):
        j = int(round(y / dt))
        j = min(j, n_steps)
        annual_prevalence_I[y] = float(I[j])

    return {
        "time": annual_years,
        "annual_incidence_time": np.arange(0, years, dtype=float),
        "annual_incidence": annual_incidence,  # counts per year
        "annual_prevalence_I": annual_prevalence_I,  # I at year boundaries (counts)
        # full trajectories at dt resolution (for chaining / debugging)
        "t": t,
        "S": S,
        "L_fast": L_fast,
        "L_slow": L_slow,
        "I": I,
        "R": R,
    }
