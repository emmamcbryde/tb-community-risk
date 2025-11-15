import numpy as np
import pandas as pd
from engine.params import extract_core_parameters

import numpy as np
import pandas as pd

import numpy as np
import pandas as pd

def simulate_community(inputs, file_path="data/parameters.xlsx"):
    """
    Static incidence model with proper 3-year intervention rollout.
    Always returns (df, summary).
    """

    # Inputs
    N = inputs.get("population", 10000)
    T = inputs.get("time_horizon", 20)
    base_incidence = float(inputs.get("user_incidence", 30))
    testing = inputs.get("testing", "None")
    treatment = inputs.get("treatment", "None")

    # Risk factors
    indigenous = inputs.get("indigenous_pct", 0) / 100
    smoker     = inputs.get("smoker_pct", 0) / 100
    diabetes   = inputs.get("diabetes_pct", 0) / 100
    renal      = inputs.get("renal_pct", 0) / 100
    immune     = inputs.get("immune_pct", 0) / 100
    recent_inf = inputs.get("recent_infection", 0) / 100

    # Multipliers
    rr_indigenous = 2.0
    rr_smoker     = 1.5
    rr_diabetes   = 1.8
    rr_renal      = 2.0
    rr_immune     = 4.0
    rr_recent_inf = 2.0

    # Combined baseline risk
    baseline_risk = (
        1
        + indigenous * (rr_indigenous - 1)
        + smoker     * (rr_smoker - 1)
        + diabetes   * (rr_diabetes - 1)
        + renal      * (rr_renal - 1)
        + immune     * (rr_immune - 1)
        + recent_inf * (rr_recent_inf - 1)
    )

    # --- 3-year rollout of interventions (static model) ---
    def rollout_factor(t):
        """
        Return a time-dependent intervention multiplier.
        """
        if t < 1:
            return 0.0
        elif t < 4:
            return (t - 0) / 3.0   # years 1–3 ramp up
        else:
            return 1.0

    # ============================================================
    # Intervention Cascade Calculator
    # ============================================================

    TESTS = {
        "None":      {"sensitivity": 0.0,  "specificity": 1.0},
        "TST":       {"sensitivity": 0.80, "specificity": 0.56},
        "IGRA":      {"sensitivity": 0.84, "specificity": 0.97},
        "QFT-Plus":  {"sensitivity": 0.85, "specificity": 0.98},
    }

    REGIMENS = {
        "None": {"efficacy": 0.0,  "completion": 0.0},
        "1HP":  {"efficacy": 0.90, "completion": 0.88},
        "3HP":  {"efficacy": 0.92, "completion": 0.82},
        "6H":   {"efficacy": 0.65, "completion": 0.55},
    }


    def compute_full_effect(
        testing: str,
        treatment: str,
        ltbi_prev: float,
        coverage_testing: float,
        coverage_treatment: float
    ):
        """
        Returns full_effect and a cascade breakdown.
        
        Args:
            testing: name of test (e.g. "TST", "IGRA")
            treatment: regimen name (e.g. "1HP", "3HP")
            ltbi_prev: baseline LTBI prevalence (0–1)
            coverage_testing: proportion tested (0–1)
            coverage_treatment: proportion accepting treatment among positives (0–1)
        """

        # --------------------------------------------------------
        # 1. Pull test and treatment parameters
        # --------------------------------------------------------
        test = TESTS.get(testing, TESTS["None"])
        regimen = REGIMENS.get(treatment, REGIMENS["None"])

        sens = test["sensitivity"]
        spec = test["specificity"]
        efficacy = regimen["efficacy"]
        completion = regimen["completion"]

        # --------------------------------------------------------
        # 2. Testing step
        # --------------------------------------------------------
        tested = coverage_testing

        # True LTBI among the population
        tp_rate = ltbi_prev * sens              # true positives
        fn_rate = ltbi_prev * (1 - sens)        # false negatives

        # False positives (important for cost/SAEs later)
        fp_rate = (1 - ltbi_prev) * (1 - spec)
        tn_rate = (1 - ltbi_prev) * spec

        # All expressed as *proportions of the total population*
        true_positives = tested * tp_rate
        false_positives = tested * fp_rate

        # --------------------------------------------------------
        # 3. Treatment cascade
        # --------------------------------------------------------
        treated_tp = true_positives * coverage_treatment
        treated_fp = false_positives * coverage_treatment  # not used yet

        completed_tp = treated_tp * completion
        completed_fp = treated_fp * completion             # for SAEs later

        # --------------------------------------------------------
        # 4. Protection
        # --------------------------------------------------------
        protected = completed_tp * efficacy

        hazard_reduction = protected
        full_effect = 1 - hazard_reduction

        # --------------------------------------------------------
        # 5. Return both the multiplier + debug cascade
        # --------------------------------------------------------
        cascade = {
            "tested": tested,
            "tp_rate": tp_rate,
            "fp_rate": fp_rate,
            "true_positives": true_positives,
            "false_positives": false_positives,
            "treated_tp": treated_tp,
            "treated_fp": treated_fp,
            "completed_tp": completed_tp,
            "completed_fp": completed_fp,
            "protected": protected,
            "hazard_reduction": hazard_reduction,
            "full_effect": full_effect,
        }

        return full_effect, cascade


    # Time steps
    years = np.arange(0, T + 1)

    # Apply rollout year-by-year
    intervention_multiplier = np.array([
        1 - rollout_factor(t) * (1 - full_effect)
        for t in years
    ])

    # Baseline incidence
    baseline_incidence = base_incidence * baseline_risk * np.exp(-0.02 * years)
    baseline_cases     = N * baseline_incidence / 100_000

    # Intervention incidence (with rollout)
    intervention_incidence = (
        base_incidence * baseline_risk * intervention_multiplier * np.exp(-0.02 * years)
    )
    intervention_cases = N * intervention_incidence / 100_000

    # Build output DataFrame
    df = pd.DataFrame({
        "Year": np.tile(years, 2),
        "Scenario": np.repeat(["Baseline", "Intervention"], len(years)),
        "Incidence_per_100k": np.concatenate([baseline_incidence, intervention_incidence]),
        "Cases": np.concatenate([baseline_cases, intervention_cases])
    })

    # Summary — consistent with app
    summary = {
        "Peak_incidence_baseline": float(np.max(baseline_incidence)),
        "Peak_incidence_intervention": float(np.max(intervention_incidence)),
        "Total_cases_baseline": float(np.sum(baseline_cases)),
        "Total_cases_intervention": float(np.sum(intervention_cases)),
    }

    return df, summary
