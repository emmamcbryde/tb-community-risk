import numpy as np


def simulate_dynamic(params, years, intervention=True):
    """
    Dynamic TB model with:
    - (Default) LTBI-based seeding from age-structured LTBI probabilities:
        ltbi_ever, ltbi_recent, age_counts
    - Optional override seeding via params["initial_state"] (for stitching calibration->projection)
    - Fast & slow latent compartments:
        L_fast progresses at 0.01 /y, L_slow at 0.001 /y
        people move L_fast -> L_slow at rate 1/5 per year
    - Risk-factor multipliers on progression
    - Force-of-infection transmission S -> L_fast
    - Pulse LTBI test & treat (fraction per year over rollout period), removing L_fast & L_slow to S
    - Diagnosis improvement (gamma) over rollout years
    - Fine time-step integration (dt = 0.1 years) with annualised outputs.

    Returns:
    - annual_incidence: new active TB cases per year (counts), length = years
    - annual_prevalence_I: I compartment at integer years (counts), length = years+1
    - compartment trajectories at dt resolution (S, L_fast, L_slow, I, R)
    - final_state dict (S, L_fast, L_slow, I, R) at end of simulation (counts)
    """

    # ---------------------------
    # Required parameters
    # ---------------------------
    required_keys = [
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
    # Beta handling: scalar beta and optional beta_series (per-year)
    # ---------------------------
    beta_series = params.get("beta_series", None)
    if beta_series is not None:
        beta_series = np.asarray(beta_series, dtype=float).reshape(-1)
        if beta_series.size == 0:
            beta_series = None

    # Scalar fallback beta (used if beta_series is None, or as fallback value)
    beta_fallback = params.get("beta", 0.0)
    beta_fallback_arr = np.asarray(beta_fallback, dtype=float)
    beta_fallback_scalar = (
        float(beta_fallback_arr)
        if beta_fallback_arr.shape == ()
        else float(beta_fallback_arr.ravel()[-1])
    )

    # ---------------------------
    # Extract parameters
    # ---------------------------
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

    # Optional: override initial state to stitch calibration -> projection
    initial_state = params.get("initial_state", None)

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

    # Effective total coverage of LTBI (over rollout period)
    coverage_total = ltbi_coverage * eff
    annual_fraction_treated = (
        coverage_total / rollout_years if rollout_years > 0 else 0.0
    )

    # ---------------------------
    # Seed compartments
    # ---------------------------
    N_target = float(sum(age_counts.values()))

    if isinstance(initial_state, dict):
        # Stitch-mode: start from an explicit state vector
        S0 = float(initial_state.get("S", 0.0))
        L_fast0 = float(initial_state.get("L_fast", 0.0))
        L_slow0 = float(initial_state.get("L_slow", 0.0))
        I0 = float(initial_state.get("I", 0.0))
        R0 = float(initial_state.get("R", 0.0))

        # numeric safety
        S0 = max(S0, 0.0)
        L_fast0 = max(L_fast0, 0.0)
        L_slow0 = max(L_slow0, 0.0)
        I0 = max(I0, 0.0)
        R0 = max(R0, 0.0)

        # Ensure population mass matches target (adjust S)
        N0 = S0 + L_fast0 + L_slow0 + I0 + R0
        if N_target > 0 and abs(N0 - N_target) > 1e-6:
            S0 = max(S0 + (N_target - N0), 0.0)

    else:
        # Default: seed from LTBI by age
        S0 = 0.0
        L_fast0 = 0.0
        L_slow0 = 0.0

        for a, pop in age_counts.items():
            p_ever = float(ltbi_ever.get(a, 0.0))
            p_recent = float(ltbi_recent.get(a, 0.0))
            S0 += pop * (1.0 - p_ever)
            L_fast0 += pop * p_recent
            L_slow0 += pop * (p_ever - p_recent)

        # Seed I0 using: prevalence ≈ incidence * duration (duration = pre_det_months)
        duration_years = pre_det_months / 12.0
        I0 = N_target * (initial_inc / 100000.0) * duration_years

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
    annual_incidence = np.zeros(int(years), dtype=float)

    # ---------------------------
    # Intervention functions
    # ---------------------------
    def tau_at_time(time_years):
        if (not intervention) or rollout_years <= 0:
            return 0.0
        if time_years <= rollout_years:
            return annual_fraction_treated
        return 0.0

    def gamma_at_time(time_years):
        if not intervention:
            return delta_pre
        if rollout_years <= 0:
            return delta_post
        frac = min(1.0, time_years / rollout_years)
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

        # Year index for piecewise-constant beta_series
        year_idx = int(np.floor(time_prev))
        if year_idx < 0:
            year_idx = 0

        if beta_series is None:
            beta_now = beta_fallback_scalar
        else:
            year_idx = min(year_idx, int(beta_series.size) - 1)
            beta_now = float(beta_series[year_idx])

        lambda_t = beta_now * I[i - 1] / N
        tau_now = tau_at_time(time_prev)
        gamma_now = gamma_at_time(time_prev)

        # New TB cases (incidence flow) this instant
        new_cases_rate = sigma_fast_eff * L_fast[i - 1] + sigma_slow_eff * L_slow[i - 1]

        # Accumulate annual cases
        if 0 <= year_idx < int(years):
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
    annual_years = np.arange(0, int(years) + 1)
    annual_prevalence_I = np.zeros(int(years) + 1, dtype=float)
    for y in annual_years:
        idx = int(round(y / dt))
        idx = min(idx, n_steps)
        annual_prevalence_I[int(y)] = I[idx]

    final_state = {
        "S": float(S[-1]),
        "L_fast": float(L_fast[-1]),
        "L_slow": float(L_slow[-1]),
        "I": float(I[-1]),
        "R": float(R[-1]),
    }

    return {
        "time": annual_years,
        "annual_incidence_time": np.arange(0, int(years), dtype=float),
        "annual_incidence": annual_incidence,  # counts per year (annual totals)
        "annual_prevalence_I": annual_prevalence_I,  # I at year boundaries
        "S": S,
        "L_fast": L_fast,
        "L_slow": L_slow,
        "I": I,
        "R": R,
        "final_state": final_state,
    }
