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
    Dynamic age-structured LTBI → TB model with:
      - S        : Susceptible
      - L_fast   : Recent latent TB
      - L_slow   : Remote latent TB
      - A_tb     : Active TB, undiagnosed (fully infectious)
      - I        : Active TB, detected/on treatment (1/6 infectiousness)
      - R        : Recovered
    """

    # ------------------------------
    # 1. FORCE OF INFECTION (base ARI)
    # ------------------------------
    base_ari = inputs.get("base_ari", None)

    if base_ari is None:
        user_incidence = inputs.get("user_incidence", 30)
        base_ari = user_incidence / 5000.0  # fallback Stýblo-like rule

    base_ari = min(base_ari, 0.10)  # safety cap

    # ------------------------------
    # 2. CASE DETECTION INPUTS (OPTION C – instantaneous change)
    # ------------------------------
    delta_pre = inputs.get("delta_pre", 1.0)    # per year
    delta_post = inputs.get("delta_post", 2.0)  # per year

    # Option C: instantaneous shift to post-intervention detection at t=0
    delta_t = delta_post

    # ------------------------------
    # 3. Load demographic and TB parameters
    # ------------------------------
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

    # ------------------------------
    # 4. Regimen cure logic (LTBI treatment)
    # ------------------------------
    treatment = inputs.get("treatment", "None")
    coverage_treatment = inputs.get("coverage_treatment", 0.0)
    completion_override = inputs.get("completion_override", None)

    reg = REGIMENS.get(treatment, REGIMENS["None"])
    efficacy = reg["efficacy"]
    completion_used = completion_override if completion_override is not None else reg["completion"]
    regimen_cure = efficacy * completion_used  # cure probability for treated LTBI

    # ------------------------------
    # 5. Initialise compartments
    # ------------------------------
    S = np.zeros((A, T + 1))
    L_fast = np.zeros((A, T + 1))
    L_slow = np.zeros((A, T + 1))
    A_tb = np.zeros((A, T + 1))     # Active TB undiagnosed
    I = np.zeros((A, T + 1))        # Active TB detected
    R = np.zeros((A, T + 1))

    N0 = np.array([age_counts[a] for a in ages])
    LTBI0 = np.array([ltbi_prev.get(a, 0.0) for a in ages])

    S[:, 0] = (1 - LTBI0) * N0
    L_fast[:, 0] = LTBI0 * N0 * 0.05
    L_slow[:, 0] = LTBI0 * N0 * 0.95
    A_tb[:, 0] = 0.0
    I[:, 0] = 0.0

    # ------------------------------
    # 6. Transition rates
    # ------------------------------
    r_fast = 0.01
    r_slow = 0.001
    tau_fast_to_slow = 1 / 5.0
    gamma_cure = 0.7  # I → R cure

    # ------------------------------
    # 7. Simulation loop
    # ------------------------------
    total_incidence = np.zeros(T + 1)

    for t in range(T):
        # Total infectious population for FOI
        infectious_pool = np.sum(A_tb[:, t]) + (1/6) * np.sum(I[:, t])
        foi_t = base_ari * infectious_pool

        for i, a in enumerate(ages):

            mu = mort.get(a, 0.0)
            mu_tb = tbmort.get(a, 0.05)

            # New i
