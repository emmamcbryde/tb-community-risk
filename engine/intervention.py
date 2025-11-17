
# engine/intervention.py

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
    "4R": {"efficacy": 0.90,"completion": 0.85},
    "9H": {"efficacy": 0.65,"completion": 0.50},
}


def compute_full_effect(
    testing: str,
    treatment: str,
    ltbi_prev: float,
    coverage_testing: float,
    coverage_treatment: float
):
    """
    Compute intervention hazard reduction using test+regimen cascade.
    Returns (full_effect, cascade_dict)
    """
    test = TESTS.get(testing, TESTS["None"])
    regimen = REGIMENS.get(treatment, REGIMENS["None"])

    sens = test["sensitivity"]
    spec = test["specificity"]
    efficacy = regimen["efficacy"]
    completion = regimen["completion"]

    # Testing
    tested = coverage_testing
    tp_rate = ltbi_prev * sens
    fp_rate = (1 - ltbi_prev) * (1 - spec)

    true_positives  = tested * tp_rate
    false_positives = tested * fp_rate

    # Treatment
    treated_tp = true_positives * coverage_treatment
    completed_tp = treated_tp * completion

    # Protection
    protected = completed_tp * efficacy

    hazard_reduction = protected
    full_effect = 1 - hazard_reduction

    cascade = {
        "tested": tested,
        "true_positives": true_positives,
        "false_positives": false_positives,
        "treated_tp": treated_tp,
        "completed_tp": completed_tp,
        "protected": protected,
        "hazard_reduction": hazard_reduction,
        "full_effect": full_effect,
    }

    return full_effect, cascade
