import numpy as np
import pandas as pd


def run_static_mechanistic(params, L_fast0=None, L_slow0=None):
    """
    One-year static mechanistic TB model with sub-year time steps (dt = 0.1):
    - If L_fast0 and L_slow0 are None, seed from LTBI distribution (year 0).
    - Otherwise, carry forward latent states from previous year.
    - Applies pulse-style LTBI Test & Treat over the year:
        LTBI → S according to coverage * regimen efficacy.
    - Computes annual incidence as flow out of LTBI.
    - Computes prevalence via diagnosis delay (pre_det_months).
    - Returns updated L_fast1, L_slow1, plus annual incidence/prevalence.
    """


    smoker_pct   = params["smoker_pct"]
    diabetes_pct = params["diabetes_pct"]
    renal_pct    = params["renal_pct"]
    HIV_untreated_pct = params["HIV_untreated_pct"]
    HIV_treated_pct   = params["HIV_treated_pct"] 
    alcohol_pct   = params["alcohol_pct"]
    ltbi_ever   = params["ltbi_ever"]        # dict age → Prob(ever infected)
    ltbi_recent = params["ltbi_recent"]      # dict age → Prob(recent infection)
    age_counts  = params["age_counts"]       # dict age → population


    pre_det_months = params["pre_det_months"]

    treatment_method = params.get("treatment_method", "None")
    ltbi_coverage    = params.get("ltbi_coverage", 0.0)

    # Natural history
    sigma_fast = 0.01
    sigma_slow = 0.001

    # Risk multiplier
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
        risk_multiplier *= (1 - p[k]) + p[k] * RR[k]


    sigma_fast_eff = sigma_fast * risk_multiplier
    sigma_slow_eff = sigma_slow * risk_multiplier

    # Treatment efficacy table
    regimen_efficacy = {
        "None": 0.0,
        "1HP": 0.90,
        "3HP": 0.90,
        "4R": 0.80,
        "6H": 0.70,
        "9H": 0.75,
    }
    eff = regimen_efficacy.get(treatment_method, 0.0)

    # Fraction of LTBI cured over a year
    coverage_total = ltbi_coverage * eff

    # For dt stepping, approximate as a constant per-year hazard:
    # treat_rate ≈ coverage_total per year (if small/moderate coverage)
    treat_rate = coverage_total

    # -----------------------
    # Seed latent compartments
    # -----------------------
    N = sum(age_counts.values())

    if L_fast0 is None or L_slow0 is None:
        S0 = 0.0
        L_fast0 = 0.0
        L_slow0 = 0.0

        for a, pop in age_counts.items():
            p_ever   = ltbi_ever.get(a,   0.0)
            p_recent = ltbi_recent.get(a, 0.0)

            S0      += pop * (1.0 - p_ever)
            L_fast0 += pop * p_recent
            L_slow0 += pop * (p_ever - p_recent)

        N = S0 + L_fast0 + L_slow0
    else:
        S0 = N - L_fast0 - L_slow0

    if N <= 0:
        N = 1.0

    # -----------------------
    # Sub-step integration (dt = 0.1)
    # -----------------------
    dt = 0.1
    steps = int(1.0 / dt)

    S = S0
    L_fast = L_fast0
    L_slow = L_slow0

    incidence_count = 0.0

    for _ in range(steps):
        # Progression out of LTBI (per year)
        prog_fast = sigma_fast_eff * L_fast
        prog_slow = sigma_slow_eff * L_slow

        # LTBI treatment "rate" (approx)
        treat_fast = treat_rate * L_fast
        treat_slow = treat_rate * L_slow

        # Apply dt
        dL_fast = -(prog_fast + treat_fast) * dt
        dL_slow = -(prog_slow + treat_slow) * dt
        dS      = (treat_fast + treat_slow) * dt

        # Update states
        L_fast += dL_fast
        L_slow += dL_slow
        S      += dS

        # Accumulate incidence flow
        incidence_count += (prog_fast + prog_slow) * dt

    L_fast1 = max(L_fast, 0.0)
    L_slow1 = max(L_slow, 0.0)
    S1 = S  # not used further except for interpretation

    # -----------------------
    # Annual prevalence = incidence × duration
    # -----------------------
    duration = pre_det_months / 12.0
    prevalence_count = incidence_count * duration

    total_pop = N

    out = {
        "population": total_pop,
        "L_fast0": L_fast0, "L_slow0": L_slow0,
        "L_fast1": L_fast1, "L_slow1": L_slow1,
        "S0": S0, "S1": S1,
        "Incidence_count": incidence_count,
        "Incidence_per100k": 100000.0 * incidence_count / total_pop,
        "Prevalence_count": prevalence_count,
        "Prevalence_per100k": 100000.0 * prevalence_count / total_pop,
    }

    return out
