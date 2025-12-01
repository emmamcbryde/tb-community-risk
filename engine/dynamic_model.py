# engine/dynamic_model.py

import numpy as np
import pandas as pd
from engine.params import extract_core_parameters
from engine.intervention import REGIMENS, TESTS


def simulate_dynamic_ltbi(
    age_counts,
    beta,
    inputs,
    file_path="data/parameters.xlsx",
):
    
    """
    Dynamic age-structured LTBI → TB model with testing and LTBI treatment.

    Compartments per age:
      - S        : Susceptible
      - L_fast   : Recent latent TB
      - L_slow   : Remote latent TB
      - A_tb     : Active TB, undiagnosed (fully infectious)
      - I        : Active TB, detected/on treatment (1/6 infectiousness)
      - R        : Recovered

    We simulate two branches:
      - Baseline:       detection = delta_pre,   LTBI treatment OFF (no test-driven cures)
      - Intervention:   detection = delta_post,  LTBI treatment ON with rollout and testing.
    """

    # ------------------------------
    # 1. Core inputs
    # ------------------------------
    T = int(inputs.get("time_horizon", 20))

    # Detection rates (per year)
    delta_pre = float(inputs.get("delta_pre", 1.0))
    delta_post = float(inputs.get("delta_post", 2.0))

    # LTBI treatment coverage (intervention only)
    coverage_treatment = float(inputs.get("coverage_treatment", 0.0))
    rollout_years = float(inputs.get("rollout_years", 3.0))

    # Testing
    testing_name = inputs.get("testing", "None")
    coverage_testing = float(inputs.get("coverage_testing", 0.0))
    test = TESTS.get(testing_name, TESTS.get("None", {"sensitivity": 0.0, "specificity": 1.0}))
    sens = float(test["sensitivity"])
    spec = float(test["specificity"])


    # ------------------------------
    # 2. Load parameters
    # ------------------------------
    params = extract_core_parameters(file_path)
    mort = params.get("mortality", {})
    tbmort = params.get("tb_mortality", {})
    react = params.get("reactivation", {})

    if not react:
        react = {a: 0.001 for a in range(0, 101)}
    if not mort:
        mort = {a: 0.01 for a in range(0, 101)}
    if not tbmort:
        tbmort = {a: 0.20 for a in range(0, 101)}

    ages = sorted(age_counts.keys())
    A_age = len(ages)

    # Initial population and LTBI
    N0_vec = np.array([age_counts[a] for a in ages])
    total_pop0 = float(np.sum(N0_vec))

    # Prefer LTBI recent/remote split if provided
    ltbi_recent_by_age = inputs.get("ltbi_recent_by_age")
    ltbi_remote_by_age = inputs.get("ltbi_remote_by_age")

    if ltbi_recent_by_age is not None and ltbi_remote_by_age is not None:
        LTBI_recent0 = np.array([ltbi_recent_by_age.get(a, 0.0) for a in ages])
        LTBI_remote0 = np.array([ltbi_remote_by_age.get(a, 0.0) for a in ages])
        LTBI_total0 = np.minimum(LTBI_recent0 + LTBI_remote0, 1.0)

        S0 = (1.0 - LTBI_total0) * N0_vec
        L_fast0 = LTBI_recent0 * N0_vec
        L_slow0 = LTBI_remote0 * N0_vec
    else:
        raise ValueError("LTBI recent/remote split is required but missing.")


    # ------------------------------
    # 3. Regimen cure logic (LTBI treatment)
    # ------------------------------
    treatment_name = inputs.get("treatment", "None")
    completion_override = inputs.get("completion_override", None)

    reg = REGIMENS.get(treatment_name, REGIMENS["None"])
    eff = float(reg["efficacy"])
    if completion_override is not None:
        comp = float(completion_override)
    else:
        comp = float(reg["completion"])
    regimen_cure = eff * comp  # cure probability per treated LTBI


    # ------------------------------------------------------
    # 4. Transition rates WITH risk-factor scaling (Option 1)
    # ------------------------------------------------------

    # Risk factors from inputs (same as static model)
    smoker = inputs.get("smoker_pct", 0) / 100.0
    diabetes = inputs.get("diabetes_pct", 0) / 100.0
    renal = inputs.get("renal_pct", 0) / 100.0
    immune = inputs.get("immune_pct", 0) / 100.0
    recent_inf = inputs.get("recent_infection", 0) / 100.0

    # Relative risks (same as static model)
    rr_smoker = 1.5
    rr_diabetes = 1.8
    rr_renal = 2.0
    rr_immune = 4.0
    
    # Combined multiplier
    baseline_risk = (
        1.0
        + smoker * (rr_smoker - 1.0)
        + diabetes * (rr_diabetes - 1.0)
        + renal * (rr_renal - 1.0)
        + immune * (rr_immune - 1.0)
    )

    # Base progression hazards
    r_fast_base = 0.01
    r_slow_base = 0.001

    # Apply multiplier
    r_fast = r_fast_base * baseline_risk
    r_slow = r_slow_base * baseline_risk

    # Other transitions remain unchanged
    tau_fast_to_slow = 1 / 5.0
    gamma_cure = 0.7

    # ------------------------------
    # 5. Helper: simulate one branch
    # ------------------------------
    def run_branch(delta_detection: float, treat_cov: float, use_testing_for_ltbi: bool):
        """
        Run the dynamic model for a given detection rate and LTBI treatment coverage.
        If use_testing_for_ltbi is False, we ignore testing effects on L_fast / L_slow.
        Returns total_incidence array of length T+1.
        """

        S = np.zeros((A_age, T + 1))
        L_fast = np.zeros((A_age, T + 1))
        L_slow = np.zeros((A_age, T + 1))
        A_tb = np.zeros((A_age, T + 1))
        I = np.zeros((A_age, T + 1))
        R = np.zeros((A_age, T + 1))

        # Initial conditions
        S[:, 0] = S0
        L_fast[:, 0] = L_fast0
        L_slow[:, 0] = L_slow0

        total_inc = np.zeros(T + 1)

        for t in range(T):
            # Current total population (for FOI scaling)
            N_t = np.sum(S[:, t] + L_fast[:, t] + L_slow[:, t] + A_tb[:, t] + I[:, t] + R[:, t])
            if N_t <= 0:
                break

            # Infectious pool: A_tb full, I at 1/6 infectiousness
            infectious_pool = np.sum(A_tb[:, t]) + (1.0 / 6.0) * np.sum(I[:, t])

            # Force of infection per susceptible person
            foi_t = beta * infectious_pool / N_t
            foi_t = min(max(foi_t, 0.0), 10)  # cap to avoid absurdly high yearly risk

            for idx, a in enumerate(ages):
                mu = mort.get(a, 0.0)
                mu_tb = tbmort.get(a, 0.05)

                S_t = S[idx, t]
                Lf_t = L_fast[idx, t]
                Ls_t = L_slow[idx, t]
                A_t = A_tb[idx, t]
                I_t = I[idx, t]
                R_t = R[idx, t]

                # New infections S → L_fast
                new_inf = S_t * foi_t

                # Progression L → A_tb
                new_TB_fast = Lf_t * r_fast
                new_TB_slow = Ls_t * r_slow
                new_TB = new_TB_fast + new_TB_slow

                # Recent LTBI → remote LTBI
                move_fs = Lf_t * tau_fast_to_slow

                # --- LTBI testing + treatment (intervention only) ---
                cure_fast = 0.0
                cure_slow = 0.0

                if use_testing_for_ltbi and treat_cov > 0.0 and testing_name != "None":
                    # Individuals tested in each compartment
                    tested_S = S_t * coverage_testing
                    tested_Lf = Lf_t * coverage_testing
                    tested_Ls = Ls_t * coverage_testing

                    # True positives (latent) and false positives (susceptible)
                    TP_fast = tested_Lf * sens
                    TP_slow = tested_Ls * sens
                    FP = tested_S * (1.0 - spec)

                    # Of true positives, fraction accept LTBI treatment
                    treated_fast = TP_fast * treat_cov
                    treated_slow = TP_slow * treat_cov
                    treated_FP = FP * treat_cov  # SAEs/costs later; no compartment change now

                    # Cures among treated latent infections
                    cure_fast = treated_fast * regimen_cure
                    cure_slow = treated_slow * regimen_cure

                    # Note: we remove ONLY the cured latent from L_fast/L_slow.
                    # Treatment failures remain latent (still in L).
                else:
                    treated_fast = 0.0
                    treated_slow = 0.0
                    treated_FP = 0.0

                # Detection A_tb → I
                detected = delta_detection * A_t

                # Deaths
                dS = S_t * mu
                dLf = Lf_t * mu
                dLs = Ls_t * mu
                dA = A_t * mu_tb
                dI = I_t * mu_tb
                dR = R_t * mu

                # Cure of detected TB
                cure_I = I_t * gamma_cure

                # ----- Updates -----
                S[idx, t + 1] = (
                    S_t
                    - dS
                    - new_inf
                    + cure_fast
                    + cure_slow
                    # FP are treated but remain in S; no compartment change here
                )

                L_fast[idx, t + 1] = (
                    Lf_t
                    + new_inf
                    - new_TB_fast
                    - move_fs
                    # only cured latent leave L_fast
                    - cure_fast
                    - dLf
                )

                L_slow[idx, t + 1] = (
                    Ls_t
                    - new_TB_slow
                    - dLs
                    # only cured latent leave L_slow
                    - cure_slow
                    + move_fs
                )

                A_tb[idx, t + 1] = (
                    A_t
                    + new_TB
                    - detected
                    - dA
                )

                I[idx, t + 1] = (
                    I_t
                    + detected
                    - cure_I
                    - dI
                )

                R[idx, t + 1] = R_t + cure_I - dR

                total_inc[t] += new_TB

        if T >= 1:
            total_inc[T] = total_inc[T - 1]
        return total_inc

    # ------------------------------
    # 6. Run baseline and intervention branches
    # ------------------------------
    # Baseline: detection = delta_pre, NO LTBI treatment effect, testing ignored
    inc_baseline = run_branch(delta_pre, treat_cov=0.0, use_testing_for_ltbi=False)

    # Intervention: detection = delta_post, full LTBI testing+treatment
    inc_interv = run_branch(delta_post, treat_cov=coverage_treatment, use_testing_for_ltbi=True)

    # ------------------------------
    # 7. Build output
    # ------------------------------
    years = np.arange(0, T + 1)
    inc_b_rate = inc_baseline / total_pop0 * 1e5
    inc_i_rate = inc_interv / total_pop0 * 1e5

    df = pd.DataFrame({
        "Year": np.concatenate([years, years]),
        "Scenario": np.repeat(["Dynamic_baseline", "Dynamic_intervention"], len(years)),
        "Incidence_per_100k": np.concatenate([inc_b_rate, inc_i_rate]),
        "Cases": np.concatenate([inc_baseline, inc_interv]),
    })

    summary = {
        "Peak_incidence_baseline": float(np.max(inc_b_rate)),
        "Peak_incidence_intervention": float(np.max(inc_i_rate)),
        "Total_cases_baseline": float(np.sum(inc_baseline)),
        "Total_cases_intervention": float(np.sum(inc_interv)),
    }

    return df, summary
