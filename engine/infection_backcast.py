import numpy as np
import pandas as pd


def calc_ari_from_incidence(
    inc_history,  # dict: year offset (0 now, −1 last year etc) → incidence per 100k
    f_ssplus=0.6,  # fraction smear‐positive of all forms (constant or dict by year)
    K=5000.0,  # divisor in classic Stýblo rule
    adjustment=1.0,  # optional multiplier for shorter infectious period etc
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


def infection_prob_by_age(ages, ari_history, clearance_rate=0.0):
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
            prod *= 1.0 - ari_history.get(
                -k, ari_history.get(min(ari_history.keys()), 0)
            )
            # optionally account for clearance
            prod *= 1.0 - clearance_rate
        prob[a] = 1.0 - prod
    return prob


def infection_prob_by_age_split(ages, ari_history, window_recent=5):
    """
    Split lifetime infection probability into recent vs remote components.

    ages: list of ages (ints)
    ari_history: dict mapping years in the past to ARI, e.g. 0, -1, -2, ...
    window_recent: number of most recent years considered "recent" infection.

    Returns three dicts:
      ever[a]   = P(ever infected by age a)
      recent[a] = P(infected at least once in the last `window_recent` years)
      remote[a] = P(infected only more than `window_recent` years ago)
    """
    ever = {}
    recent = {}
    remote = {}

    min_key = min(ari_history.keys())

    for a in ages:
        if a <= 0:
            ever[a] = 0.0
            recent[a] = 0.0
            remote[a] = 0.0
            continue

        # 1. P(no infection over whole life)
        prod_all = 1.0
        for k in range(a):
            ari_k = ari_history.get(-k, ari_history.get(min_key, 0.0))
            prod_all *= 1.0 - ari_k
        P_ever = 1.0 - prod_all

        # 2. P(no infection in last `window_recent` years)
        max_recent_k = min(window_recent, a)
        prod_last = 1.0
        for k in range(max_recent_k):
            ari_k = ari_history.get(-k, ari_history.get(min_key, 0.0))
            prod_last *= 1.0 - ari_k

        # 3. P(no infection > window_recent years ago)
        prod_old = 1.0
        if a > window_recent:
            for k in range(window_recent, a):
                ari_k = ari_history.get(-k, ari_history.get(min_key, 0.0))
                prod_old *= 1.0 - ari_k

        # Infection only > window_recent years ago AND none in recent window
        P_remote_only = (1.0 - prod_old) * prod_last

        # Infected at least once in recent window
        P_recent_any = max(P_ever - P_remote_only, 0.0)  # numeric safety

        ever[a] = P_ever
        recent[a] = P_recent_any
        remote[a] = P_remote_only

    return ever, recent, remote
