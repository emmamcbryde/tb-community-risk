import numpy as np
import pandas as pd


def run_static_mechanistic(params):
    """
    One-year static mechanistic TB model:
    - Uses LTBI seeding (L_fast, L_slow) exactly as dynamic model
    - Applies pulse-style treatment to remove LTBI → S
    - Computes incidence as progression out of LTBI
    - Computes prevalence using diagnosis delay
    """

    # Extract core parameters
    smoker_pct   = params["smoker_pct"]
    diabetes_pct = params["diabetes_pct"]
    renal_pct    = params["renal_pct"]
    immune_pct   = params["immune_pct"]

    ltbi_ever   = params["ltbi_ever"]
    ltbi_recent = params["ltbi_recent"]
    age_counts  = params["age_counts"]

    delta_pre  = params["delta_pre"]
    delta_post = params["delta_post"]
    pre_det_months = params["pre_det_months"]

    treatment_method = params.get("treatment_method", "None")
    ltbi_coverage    = params.get("ltbi_coverage", 0.0)

    # Natural history
    sigma_fast = 0.01
    sigma_slow = 0.001
    omega      = 1/5   # not used in static model

    # Risk multiplier
    risk_multiplier = (
        1
        + 0.01 * smoker_pct
        + 0.02 * diabetes_pct
        + 0.03 * renal_pct
        + 0.05 * immune_pct
    )

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

    coverage_total = ltbi_coverage * eff  # total fraction of LTBI cured this year

    # -----------------------
    # Seed latent compartments
    # -----------------------
    S0 = 0
    L_fast0 = 0
    L_slow0 = 0

    for a, pop in age_counts.items():
        p_ever   = ltbi_ever.get(a,   0.0)
        p_recent = ltbi_recent.get(a, 0.0)

        S0      += pop * (1 - p_ever)
        L_fast0 += pop * p_recent
        L_slow0 += pop * (p_ever - p_recent)

    N = S0 + L_fast0 + L_slow0

    # -----------------------
    # Apply LTBI Test & Treat
    # -----------------------
    L_fast_treated = coverage_total * L_fast0
    L_slow_treated = coverage_total * L_slow0

    L_fast1 = L_fast0 - L_fast_treated
    L_slow1 = L_slow0 - L_slow_treated
    S1 = S0 + L_fast_treated + L_slow_treated

    # -----------------------
    # Incidence = L_fast1*σ_fast + L_slow1*σ_slow
    # -----------------------
    incidence_count = (
        sigma_fast_eff * L_fast1 +
        sigma_slow_eff * L_slow1
    )

    # Annual prevalence = incidence × duration
    duration = pre_det_months / 12
    prevalence_count = incidence_count * duration

    total_pop = N

    out = {
        "population": total_pop,
        "L_fast0": L_fast0, "L_slow0": L_slow0,
        "L_fast1": L_fast1, "L_slow1": L_slow1,
        "S0": S0, "S1": S1,
        "Incidence_count": incidence_count,
        "Incidence_per100k": 100000 * incidence_count / total_pop,
        "Prevalence_count": prevalence_count,
        "Prevalence_per100k": 100000 * prevalence_count / total_pop,
    }

    return out
