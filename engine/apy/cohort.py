from __future__ import annotations

from typing import Any

import numpy as np

from engine.apy.age_distribution import broad_age_group_from_years
from engine.apy.calibration import age_cumulative_infection_hazard
from engine.apy.timing import (
    average_survival_to_uniform_screen,
    preventable_active_risk_to_uniform_screen,
)


def make_rng(seed: int | None = None) -> np.random.Generator:
    return np.random.default_rng(seed)


def discrete_draw_values(values, probabilities, n: int, rng) -> np.ndarray:
    value_array = np.asarray(values)
    probabilities_array = np.asarray(probabilities, dtype=float)
    probabilities_array = probabilities_array / probabilities_array.sum()
    indices = rng.choice(len(value_array), size=int(n), p=probabilities_array)
    return value_array[indices]


def age_lookup(values_by_age_group, age_group) -> np.ndarray:
    values = np.asarray(values_by_age_group)
    groups = np.asarray(age_group, dtype=int)
    return values[groups - 1]


def exprnd_local(mean, rng) -> np.ndarray:
    mean_array = np.asarray(mean, dtype=float)
    return rng.exponential(scale=mean_array)


def infection_probability(
    age_years,
    pars: dict[str, Any],
    log_lambda: float,
    gamma: float,
    marijuana,
    contact,
    renal,
) -> np.ndarray:
    base_cum_hazard = age_cumulative_infection_hazard(age_years, log_lambda, gamma)
    eta = (
        np.log(base_cum_hazard)
        + np.log(pars["infOR"]["MJ"]) * np.asarray(marijuana, dtype=float)
        + np.log(pars["infOR"]["contact"]) * np.asarray(contact, dtype=float)
        + np.log(pars["infOR"]["renal"]) * np.asarray(renal, dtype=float)
    )
    return 1.0 - np.exp(-np.exp(eta))


def disease_multiplier(
    pars: dict[str, Any],
    marijuana,
    contact,
    renal,
    diabetes,
    smoking,
    cld,
    alcohol,
) -> np.ndarray:
    dis_or = pars["disOR"]
    return (
        (dis_or["MJ"] ** np.asarray(marijuana, dtype=float))
        * (dis_or["contact"] ** np.asarray(contact, dtype=float))
        * (dis_or["renal"] ** np.asarray(renal, dtype=float))
        * (dis_or["diabetes"] ** np.asarray(diabetes, dtype=float))
        * (dis_or["smoking"] ** np.asarray(smoking, dtype=float))
        * (dis_or["cld"] ** np.asarray(cld, dtype=float))
        * (dis_or["alcohol"] ** np.asarray(alcohol, dtype=float))
    )


def average_survival_to_random_screen(
    mult_disease,
    lambda_early: float,
    screen_window: float,
    lambda_late: float | None = None,
    early_progression_period_years: float | None = None,
) -> np.ndarray:
    if lambda_late is not None and early_progression_period_years is not None:
        return average_survival_to_uniform_screen(
            mult_disease,
            lambda_early,
            lambda_late,
            screen_window,
            early_progression_period_years,
        )
    rate = float(lambda_early) * np.asarray(mult_disease, dtype=float)
    avg_survival = np.ones_like(rate, dtype=float)
    idx = rate > 0
    denominator = np.maximum(rate[idx] * screen_window, np.finfo(float).eps)
    avg_survival[idx] = (1.0 - np.exp(-rate[idx] * screen_window)) / denominator
    return avg_survival


def preventable_active_risk(
    mult_disease,
    lambda_early: float,
    lambda_late: float,
    screen_window: float,
    follow_horizon: float,
    early_progression_period_years: float | None = None,
) -> np.ndarray:
    if early_progression_period_years is not None:
        return preventable_active_risk_to_uniform_screen(
            mult_disease,
            lambda_early,
            lambda_late,
            screen_window,
            follow_horizon,
            early_progression_period_years,
        )
    mult = np.asarray(mult_disease, dtype=float)
    avg_surv_to_screen = average_survival_to_random_screen(
        mult, lambda_early, screen_window
    )
    surv_to_horizon = np.exp(
        -(lambda_early * mult * screen_window)
        - (lambda_late * mult * max(follow_horizon - screen_window, 0))
    )
    return np.maximum(avg_surv_to_screen - surv_to_horizon, 0)


def select_screened_people(
    ltbi_risk_score,
    cure_target_score,
    prevent_target_score,
    screen_coverage: float,
    strategy: str,
    rng,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    strategy_key = str(strategy).strip().lower()
    if strategy_key == "random":
        priority_score = rng.random(len(ltbi_risk_score))
    elif strategy_key == "ltbi":
        priority_score = np.asarray(ltbi_risk_score, dtype=float)
    elif strategy_key == "cure":
        priority_score = np.asarray(cure_target_score, dtype=float)
    elif strategy_key == "prevent":
        priority_score = np.asarray(prevent_target_score, dtype=float)
    else:
        raise ValueError(f"Unknown screeningStrategy: {strategy}")

    n = len(priority_score)
    n_to_screen = min(max(round(float(screen_coverage) * n), 0), n)
    tie_breaker = rng.random(n) * 1e-12
    order = np.argsort(-(priority_score + tie_breaker))
    priority_rank = np.empty(n, dtype=int)
    priority_rank[order] = np.arange(1, n + 1)
    screened = np.zeros(n, dtype=bool)
    if n_to_screen > 0:
        screened[order[:n_to_screen]] = True
    return screened, priority_score, priority_rank


def get_test_performance(bcg, opts: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    bcg_array = np.asarray(bcg, dtype=bool)
    sensitivity = np.zeros_like(bcg_array, dtype=float)
    specificity = np.zeros_like(bcg_array, dtype=float)
    if str(opts.get("testType", "IGRA")).upper() == "IGRA":
        sensitivity[:] = float(opts.get("testSensitivity", 0.95))
        specificity[:] = float(opts.get("testSpecificity", 0.98))
    else:
        sensitivity[:] = float(opts.get("tstSensitivity", 0.80))
        specificity[~bcg_array] = float(opts.get("tstSpecificityNoBCG", 0.97))
        specificity[bcg_array] = float(opts.get("tstSpecificityBCG", 0.55))
    return sensitivity, specificity


def get_counterfactual_no_bcg_specificity(bcg, opts: dict[str, Any]) -> np.ndarray:
    bcg_array = np.asarray(bcg, dtype=bool)
    specificity = np.zeros_like(bcg_array, dtype=float)
    if str(opts.get("testType", "IGRA")).upper() == "IGRA":
        specificity[:] = float(opts.get("testSpecificity", 0.98))
    else:
        specificity[:] = float(opts.get("tstSpecificityNoBCG", 0.97))
    return specificity


def draw_untreated_active_times(
    infected,
    mult_disease,
    lambda_early: float,
    lambda_late: float,
    screen_window: float,
    rng,
    early_progression_period_years: float | None = None,
) -> np.ndarray:
    infected_array = np.asarray(infected, dtype=bool)
    mult = np.asarray(mult_disease, dtype=float)
    t_active = np.full(infected_array.shape, np.inf, dtype=float)
    idx_inf = np.flatnonzero(infected_array)
    if len(idx_inf) == 0:
        return t_active

    early_period = float(early_progression_period_years or screen_window)
    rate1 = lambda_early * mult[idx_inf]
    t1 = exprnd_local(1.0 / rate1, rng)
    early = t1 <= early_period
    t_active[idx_inf[early]] = t1[early]
    idx_late = idx_inf[~early]
    if len(idx_late) > 0:
        rate2 = lambda_late * mult[idx_late]
        t2 = exprnd_local(1.0 / rate2, rng)
        t_active[idx_late] = early_period + t2
    return t_active


def draw_base_population(
    n: int,
    pars: dict[str, Any],
    log_lambda: float,
    gamma: float,
    rng,
) -> dict[str, np.ndarray]:
    age_years = discrete_draw_values(pars["exactAgeValues"], pars["exactAgeProb"], n, rng)
    age_group = np.asarray(broad_age_group_from_years(age_years), dtype=int)

    female = rng.random(n) < pars["totalFemalePrev"]
    bcg = rng.random(n) < pars["totalBCGPrev"]

    p_mj = age_lookup(pars["mjPrevByAge"], age_group)
    p_contact = age_lookup(pars["contactPrevByAge"], age_group)
    p_renal = age_lookup(pars["renalPrevByAge"], age_group)
    p_diabetes = age_lookup(pars["diabetesPrevByAge"], age_group)
    p_smoking = age_lookup(pars["smokingPrevByAge"], age_group)
    p_cld = age_lookup(pars["cldPrevByAge"], age_group)
    p_alcohol = age_lookup(pars["alcoholPrevByAge"], age_group)

    marijuana = rng.random(n) < p_mj
    contact = rng.random(n) < p_contact
    renal = rng.random(n) < p_renal
    diabetes = rng.random(n) < p_diabetes
    smoking = rng.random(n) < p_smoking
    cld = rng.random(n) < p_cld
    alcohol = rng.random(n) < p_alcohol

    p_inf = infection_probability(
        age_years, pars, log_lambda, gamma, marijuana, contact, renal
    )
    infected = rng.random(n) < p_inf
    mult_disease = disease_multiplier(
        pars, marijuana, contact, renal, diabetes, smoking, cld, alcohol
    )

    return {
        "ageYears": age_years,
        "ageGroup": age_group,
        "female": female,
        "BCG": bcg,
        "MJ": marijuana,
        "contact": contact,
        "renal": renal,
        "diabetes": diabetes,
        "smoking": smoking,
        "cld": cld,
        "alcohol": alcohol,
        "pInfection": p_inf,
        "infected": infected,
        "diseaseMultiplier": mult_disease,
    }


def add_targeting_scores(
    population: dict[str, np.ndarray],
    lambda_early: float,
    lambda_late: float,
    screen_window: float,
    follow_horizon: float,
    early_progression_period_years: float | None = None,
) -> dict[str, np.ndarray]:
    out = dict(population)
    avg_latent_at_screen = average_survival_to_random_screen(
        out["diseaseMultiplier"],
        lambda_early,
        screen_window,
        lambda_late,
        early_progression_period_years,
    )
    preventable = preventable_active_risk(
        out["diseaseMultiplier"],
        lambda_early,
        lambda_late,
        screen_window,
        follow_horizon,
        early_progression_period_years,
    )
    out["ltbiRiskScore"] = out["pInfection"]
    out["cureTargetScore"] = out["pInfection"] * avg_latent_at_screen
    out["preventTargetScore"] = out["pInfection"] * preventable
    return out
