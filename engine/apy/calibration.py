from __future__ import annotations

import math
from itertools import product
from typing import Any, Iterable

import numpy as np

from engine.apy.age_distribution import broad_age_group_from_years
from engine.apy.config import normalise_config
from engine.apy.data import load_parameters_from_config
from engine.apy.ltbi_state import (
    RECENT_TO_REMOTE_RATE_PER_YEAR,
    mixed_baseline_event_between,
    require_numeric_ltbi_state_assumptions,
)
from engine.apy.timing import resolve_time_settings


EPS = np.finfo(float).eps
MATLAB_V9_IMPLICIT_EARLY_LATE = "matlab_v9_implicit_early_late"


def bern_prob(x: int | bool, p: float) -> float:
    return float(p) if bool(x) else 1.0 - float(p)


def odds_from_prob(p: float) -> float:
    return float(p) / max(1.0 - float(p), EPS)


def age_cumulative_infection_hazard(age, log_lambda: float, gamma: float):
    hazard = np.exp(log_lambda) * ((np.asarray(age, dtype=float) + 0.5) ** gamma)
    return np.maximum(hazard, EPS)


def infection_probability_from_eta(eta):
    return 1.0 - np.exp(-np.exp(eta))


def expected_infection_prevalence_exact(
    log_lambda: float,
    gamma: float,
    pars: dict[str, Any],
) -> float:
    return _expected_infection_for_ages(
        log_lambda,
        gamma,
        pars,
        pars["exactAgeValues"],
        pars["exactAgeProb"],
    )


def expected_infection_in_age_subset(
    log_lambda: float,
    gamma: float,
    pars: dict[str, Any],
    age_mask,
) -> float:
    ages = list(pars["exactAgeValues"])
    probs = list(pars["exactAgeProb"])
    mask = list(age_mask)
    subset_ages = [age for age, keep in zip(ages, mask) if keep]
    subset_probs = [prob for prob, keep in zip(probs, mask) if keep]
    total = sum(subset_probs)
    if total <= 0:
        return math.nan
    subset_probs = [p / total for p in subset_probs]
    return _expected_infection_for_ages(log_lambda, gamma, pars, subset_ages, subset_probs)


def expected_infection_by_age_threshold(
    log_lambda: float,
    gamma: float,
    pars: dict[str, Any],
    age_cut: float = 25,
) -> tuple[float, float]:
    ages = list(pars["exactAgeValues"])
    older_mask = [age >= age_cut for age in ages]
    younger_mask = [age < age_cut for age in ages]
    return (
        expected_infection_in_age_subset(log_lambda, gamma, pars, older_mask),
        expected_infection_in_age_subset(log_lambda, gamma, pars, younger_mask),
    )


def expected_age_or_exact(
    log_lambda: float,
    gamma: float,
    pars: dict[str, Any],
    age_cut: float = 25,
) -> float:
    p_older, p_younger = expected_infection_by_age_threshold(
        log_lambda, gamma, pars, age_cut
    )
    return odds_from_prob(p_older) / max(odds_from_prob(p_younger), EPS)


def solve_loglambda_for_gamma(
    gamma: float,
    pars: dict[str, Any],
    target_inf_prev: float,
) -> float:
    def objective(log_lambda: float) -> float:
        return (
            expected_infection_prevalence_exact(log_lambda, gamma, pars)
            - target_inf_prev
        )

    lo = -30.0
    hi = 10.0
    f_lo = objective(lo)
    f_hi = objective(hi)
    while f_lo > 0:
        lo -= 5.0
        f_lo = objective(lo)
    while f_hi < 0:
        hi += 5.0
        f_hi = objective(hi)
    return _bisect_root(objective, lo, hi)


def calibrate_age_infection_model(
    pars: dict[str, Any],
    target_inf_prev: float,
    target_age_or: float,
) -> dict[str, float]:
    def objective(log_gamma: float) -> float:
        gamma = math.exp(log_gamma)
        log_lambda = solve_loglambda_for_gamma(gamma, pars, target_inf_prev)
        return math.log(expected_age_or_exact(log_lambda, gamma, pars) / target_age_or)

    grid = np.linspace(-3.0, 3.0, 49)
    values = [_safe_eval(objective, x) for x in grid]
    log_gamma = None
    for left, right, f_left, f_right in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if math.isfinite(f_left) and math.isfinite(f_right) and np.sign(f_left) != np.sign(f_right):
            log_gamma = _bisect_root(objective, float(left), float(right))
            break
    if log_gamma is None:
        log_gamma = _minimise_scalar_grid(lambda x: objective(x) ** 2, -3.0, 3.0)

    gamma = math.exp(log_gamma)
    log_lambda = solve_loglambda_for_gamma(gamma, pars, target_inf_prev)
    return {"logLambda": log_lambda, "gamma": gamma}


def disease_multiplier_from_flags(
    pars: dict[str, Any],
    mj=False,
    contact=False,
    renal=False,
    diabetes=False,
    smoking=False,
    cld=False,
    alcohol=False,
):
    dis_or = pars["disOR"]
    return (
        (dis_or["MJ"] ** int(bool(mj)))
        * (dis_or["contact"] ** int(bool(contact)))
        * (dis_or["renal"] ** int(bool(renal)))
        * (dis_or["diabetes"] ** int(bool(diabetes)))
        * (dis_or["smoking"] ** int(bool(smoking)))
        * (dis_or["cld"] ** int(bool(cld)))
        * (dis_or["alcohol"] ** int(bool(alcohol)))
    )


def expected_active_within_window(
    lambda_early: float,
    pars: dict[str, Any],
    log_lambda: float,
    gamma: float,
    screen_window: float,
    early_progression_period_years: float | None = None,
    lambda_late: float | None = None,
    baseline_recent_proportion: float = 1.0,
    recent_to_remote_rate: float | None = None,
) -> float:
    early_period = float(early_progression_period_years or screen_window)
    late = float(lambda_late if lambda_late is not None else lambda_early)
    prev = 0.0
    for age, age_prob in zip(pars["exactAgeValues"], pars["exactAgeProb"]):
        age_group = broad_age_group_from_years([age])[0] - 1
        p_mj = pars["mjPrevByAge"][age_group]
        p_contact = pars["contactPrevByAge"][age_group]
        p_renal = pars["renalPrevByAge"][age_group]
        p_diabetes = pars["diabetesPrevByAge"][age_group]
        p_smoking = pars["smokingPrevByAge"][age_group]
        p_cld = pars["cldPrevByAge"][age_group]
        p_alcohol = pars["alcoholPrevByAge"][age_group]
        p_inf_by_flags = _infection_probability_for_age_and_primary_flags(
            age, log_lambda, gamma, pars
        )
        for mj, contact, renal, diabetes, smoking, cld, alcohol in product([0, 1], repeat=7):
            p_inf = p_inf_by_flags[(mj, contact, renal)]
            mult = disease_multiplier_from_flags(
                pars, mj, contact, renal, diabetes, smoking, cld, alcohol
            )
            p_active = float(
                mixed_baseline_event_between(
                    0.0,
                    screen_window,
                    mult,
                    lambda_early,
                    late,
                    RECENT_TO_REMOTE_RATE_PER_YEAR
                    if recent_to_remote_rate is None
                    else recent_to_remote_rate,
                    baseline_recent_proportion,
                )
            )
            prev += (
                age_prob
                * bern_prob(mj, p_mj)
                * bern_prob(contact, p_contact)
                * bern_prob(renal, p_renal)
                * bern_prob(diabetes, p_diabetes)
                * bern_prob(smoking, p_smoking)
                * bern_prob(cld, p_cld)
                * bern_prob(alcohol, p_alcohol)
                * p_inf
                * p_active
            )
    return prev


def expected_active_within_window_matlab_v9(
    lambda_early: float,
    pars: dict[str, Any],
    log_lambda: float,
    gamma: float,
    screen_window: float,
) -> float:
    """MATLAB v9 compatibility target: all infections begin in the early phase.

    This preserves the reference software semantics for stored APY v9
    validation scenarios. It is not a scientific recent-LTBI estimate.
    """
    weights, multipliers = _progression_weighted_multipliers(pars, log_lambda, gamma)
    return float((weights * (1.0 - np.exp(-float(screen_window) * lambda_early * multipliers))).sum())


def calibrate_early_hazard(
    pars: dict[str, Any],
    log_lambda: float,
    gamma: float,
    target_active_2y: float,
    screen_window: float,
    early_late_ratio: float,
    early_progression_period_years: float | None = None,
    baseline_recent_proportion: float = 1.0,
    recent_to_remote_rate: float | None = None,
    natural_history_semantics: str | None = None,
) -> dict[str, float]:
    if early_late_ratio < 1:
        raise ValueError("earlyLateRatio must be >= 1")
    weights, multipliers = _progression_weighted_multipliers(pars, log_lambda, gamma)

    def objective(lambda_early: float) -> float:
        lambda_late = lambda_early / early_late_ratio
        if natural_history_semantics == MATLAB_V9_IMPLICIT_EARLY_LATE:
            return (
                expected_active_within_window_matlab_v9(
                    lambda_early,
                    pars,
                    log_lambda,
                    gamma,
                    screen_window,
                )
                - target_active_2y
            )
        return (
            float(
                (
                    weights
                    * mixed_baseline_event_between(
                        0.0,
                        screen_window,
                        multipliers,
                        lambda_early,
                        lambda_late,
                        RECENT_TO_REMOTE_RATE_PER_YEAR
                        if recent_to_remote_rate is None
                        else recent_to_remote_rate,
                        baseline_recent_proportion,
                    )
                ).sum()
            )
            - target_active_2y
        )

    lambda_early = _bisect_root(objective, 1e-8, 10.0)
    return {
        "lambdaEarly": lambda_early,
        "lambdaLate": lambda_early / early_late_ratio,
        "expectedActiveAtCalibrationHorizon": target_active_2y + objective(lambda_early),
    }


def _progression_weighted_multipliers(
    pars: dict[str, Any],
    log_lambda: float,
    gamma: float,
) -> tuple[np.ndarray, np.ndarray]:
    weights = []
    multipliers = []
    for age, age_prob in zip(pars["exactAgeValues"], pars["exactAgeProb"]):
        age_group = broad_age_group_from_years([age])[0] - 1
        p_mj = pars["mjPrevByAge"][age_group]
        p_contact = pars["contactPrevByAge"][age_group]
        p_renal = pars["renalPrevByAge"][age_group]
        p_diabetes = pars["diabetesPrevByAge"][age_group]
        p_smoking = pars["smokingPrevByAge"][age_group]
        p_cld = pars["cldPrevByAge"][age_group]
        p_alcohol = pars["alcoholPrevByAge"][age_group]
        p_inf_by_flags = _infection_probability_for_age_and_primary_flags(
            age, log_lambda, gamma, pars
        )
        for mj, contact, renal, diabetes, smoking, cld, alcohol in product([0, 1], repeat=7):
            weights.append(
                age_prob
                * bern_prob(mj, p_mj)
                * bern_prob(contact, p_contact)
                * bern_prob(renal, p_renal)
                * bern_prob(diabetes, p_diabetes)
                * bern_prob(smoking, p_smoking)
                * bern_prob(cld, p_cld)
                * bern_prob(alcohol, p_alcohol)
                * p_inf_by_flags[(mj, contact, renal)]
            )
            multipliers.append(
                disease_multiplier_from_flags(
                    pars, mj, contact, renal, diabetes, smoking, cld, alcohol
                )
            )
    return np.asarray(weights, dtype=float), np.asarray(multipliers, dtype=float)


def calibrate_from_config(config: dict[str, Any]) -> dict[str, Any]:
    cfg = normalise_config(config)
    pars = load_parameters_from_config(cfg)
    target_inf_prev = _default_if_empty(cfg["ltbiPrevalence"], 47 / 624)
    target_active_2y = (
        _default_if_empty(cfg["activeTBPrevalence"], 10 / 770)
    )
    early_late_ratio = _default_if_empty(cfg["earlyLateRatio"], 5)
    timing = resolve_time_settings(cfg)
    ltbi_state = require_numeric_ltbi_state_assumptions(cfg)

    age_calibration = calibrate_age_infection_model(
        pars, target_inf_prev, cfg["targetAgeOR"]
    )
    hazard = calibrate_early_hazard(
        pars,
        age_calibration["logLambda"],
        age_calibration["gamma"],
        target_active_2y,
        timing["activeTBCalibrationHorizonYears"],
        early_late_ratio,
        timing["earlyProgressionPeriodYears"],
        ltbi_state["baselineRecentLTBIProportion"],
        ltbi_state["recentToRemoteTransitionRatePerYear"],
        cfg.get("naturalHistorySemantics"),
    )
    return {
        "parameters": pars,
        "ageInfLogLambda": age_calibration["logLambda"],
        "ageInfGamma": age_calibration["gamma"],
        "expectedInfPrev": expected_infection_prevalence_exact(
            age_calibration["logLambda"], age_calibration["gamma"], pars
        ),
        "expectedAgeOR": expected_age_or_exact(
            age_calibration["logLambda"], age_calibration["gamma"], pars
        ),
        "lambdaEarly": hazard["lambdaEarly"],
        "lambdaLate": hazard["lambdaLate"],
        "expectedActiveAtCalibrationHorizon": hazard["expectedActiveAtCalibrationHorizon"],
        "expectedActive2y": hazard["expectedActiveAtCalibrationHorizon"],
        "targetInfPrev": target_inf_prev,
        "targetAgeOR": cfg["targetAgeOR"],
        "targetActiveAtCalibrationHorizon": target_active_2y,
        "targetActive2y": target_active_2y,
        "earlyProgressionPeriodYears": timing["earlyProgressionPeriodYears"],
        "activeTBCalibrationHorizonYears": timing["activeTBCalibrationHorizonYears"],
        "baselineRecentLTBIProportion": ltbi_state["baselineRecentLTBIProportion"],
        "recentToRemoteTransitionRatePerYear": ltbi_state[
            "recentToRemoteTransitionRatePerYear"
        ],
        "ltbiStateAssumptionStatus": ltbi_state["status"],
    }


def _expected_infection_for_ages(
    log_lambda: float,
    gamma: float,
    pars: dict[str, Any],
    ages: Iterable[float],
    age_probs: Iterable[float],
) -> float:
    prev = 0.0
    for age, age_prob in zip(ages, age_probs):
        p_inf_by_flags = _infection_probability_for_age_and_primary_flags(
            age, log_lambda, gamma, pars
        )
        age_group = broad_age_group_from_years([age])[0] - 1
        p_mj = pars["mjPrevByAge"][age_group]
        p_contact = pars["contactPrevByAge"][age_group]
        p_renal = pars["renalPrevByAge"][age_group]
        for mj, contact, renal in product([0, 1], repeat=3):
            prev += (
                age_prob
                * bern_prob(mj, p_mj)
                * bern_prob(contact, p_contact)
                * bern_prob(renal, p_renal)
                * p_inf_by_flags[(mj, contact, renal)]
            )
    return prev


def _infection_probability_for_age_and_primary_flags(
    age: float, log_lambda: float, gamma: float, pars: dict[str, Any]
) -> dict[tuple[int, int, int], float]:
    base_cum_hazard = float(age_cumulative_infection_hazard(age, log_lambda, gamma))
    out = {}
    for mj, contact, renal in product([0, 1], repeat=3):
        eta = (
            math.log(base_cum_hazard)
            + math.log(pars["infOR"]["MJ"]) * mj
            + math.log(pars["infOR"]["contact"]) * contact
            + math.log(pars["infOR"]["renal"]) * renal
        )
        out[(mj, contact, renal)] = float(infection_probability_from_eta(eta))
    return out


def _bisect_root(func, lo: float, hi: float, *, tol: float = 1e-12, max_iter: int = 200) -> float:
    f_lo = func(lo)
    f_hi = func(hi)
    if f_lo == 0:
        return lo
    if f_hi == 0:
        return hi
    if np.sign(f_lo) == np.sign(f_hi):
        raise ValueError("Root is not bracketed.")
    for _ in range(max_iter):
        mid = (lo + hi) / 2.0
        f_mid = func(mid)
        if abs(f_mid) < tol or (hi - lo) / 2.0 < tol:
            return mid
        if np.sign(f_mid) == np.sign(f_lo):
            lo = mid
            f_lo = f_mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def _minimise_scalar_grid(func, lo: float, hi: float) -> float:
    grid = np.linspace(lo, hi, 801)
    values = [_safe_eval(func, x) for x in grid]
    best_index = int(np.nanargmin(values))
    return float(grid[best_index])


def _safe_eval(func, x: float) -> float:
    try:
        return float(func(float(x)))
    except (ValueError, OverflowError, ZeroDivisionError):
        return math.nan


def _default_if_empty(value: Any, default: float) -> float:
    if value is None or value == []:
        return float(default)
    return float(value)
