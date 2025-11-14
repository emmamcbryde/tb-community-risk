import numpy as np
import pandas as pd

def calc_ari_from_incidence(
    inc_history,     # dict: year offset (0 now, −1 last year etc) → incidence per 100k
    f_ssplus=0.6,    # fraction smear‐positive of all forms (constant or dict by year)
    K=5000.0,        # divisor in classic Stýblo rule
    adjustment=1.0   # optional multiplier for shorter infectious period etc
):
    """
    Converts incidence history to ARI history (probability) via modernised Stýblo‐style rule.
    Returns dict: year offset → ARI.
    """
    ari = {}
    for t, inc in inc_history.items():
        ari_t = (inc * f_ssplus / K) * adjustment
        ari[t] = min(max(ari_t, 0.0), 1.0)  # clamp between 0 and 1
    return ari

def infection_prob_by_age(
    ages,
    ari_history,
    clearance_rate=0.0
):
    """
    Computes P(ever infected by now) for each age in `ages`, using ARI history.
    ari_history keys: 0 (now), -1, -2, ..., −(maxAge‐1)
    clearance_rate: annual self‐clearance probability (default 0)
    Returns dict: age → probability (0-1)
    """
    prob = {}
    for a in ages:
        if a <= 0:
            prob[a] = 0.0
            continue
        prod = 1.0
        for k in range(a):
            prod *= (1.0 - ari_history.get(-k, ari_history.get(min(ari_history.keys()), 0)))
            # optionally account for clearance
            prod *= (1.0 - clearance_rate)
        prob[a] = 1.0 - prod
    return prob
