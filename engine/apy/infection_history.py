from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import product
from typing import Any

import numpy as np

from engine.apy.age_distribution import broad_age_group_from_years


_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(16)
_UNIT_NODES = (_GL_NODES + 1.0) / 2.0
_UNIT_WEIGHTS = _GL_WEIGHTS / 2.0
_AGE_RISK_GRID_CACHE: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
_CALIBRATION_CACHE: dict[tuple[Any, ...], dict[str, Any]] = {}

TRAJECTORY_MODEL = "calibrated_historical_infection_pressure"
DERIVATION_METHOD = "infection_history_trajectory"
COMPATIBILITY_REFERENCE_BASIS = "sa_health_matlab_v9_compatibility_reference"
MATLAB_V9_COMPATIBILITY_SEMANTICS = "matlab_v9_implicit_early_late"
EXPLICIT_EXPERIMENTAL_BASIS = "explicit_recent_remote_scientific_scenario"
DEFAULT_RECENT_WINDOW_YEARS = 2.0
DEFAULT_TREND_ABS_PER_YEAR = 0.01
EXPERIMENTAL_STATUS_LABEL = "Experimental - not the report reference"
EXPERIMENTAL_ECONOMICS_GUARD_MESSAGE = (
    "Health-economic conclusions are not available for this experimental "
    "infection-history scenario because future TB progression calibration is not "
    "yet validated. Restore the SA Health reference to reproduce the report analysis."
)
EXPERIMENTAL_READINESS_ITEMS = [
    {
        "item": "Provenance and interpretation of 10/770",
        "status": "unresolved",
        "reason": "The inherited target is not yet established as future incident progression from LTBI.",
    },
    {
        "item": "Baseline active TB separated from future incident progression",
        "status": "unresolved",
        "reason": "Prevalent or near-baseline TB is not represented as a distinct compartment.",
    },
    {
        "item": "Recent and remote progression hazards",
        "status": "unresolved",
        "reason": "Hazards are software-calibrated scenario values, not reviewed natural-history estimates.",
    },
    {
        "item": "Two-year infection timing versus early higher-risk state duration",
        "status": "unresolved",
        "reason": "Infection acquired within two years and a five-year mean early-risk state are related but not equivalent definitions.",
    },
    {
        "item": "Age odds-ratio interpretation",
        "status": "unresolved",
        "reason": "The implemented age association is an odds ratio for LTBI prevalence in >=25 versus <25 years.",
    },
    {
        "item": "Joint risk-factor multiplier handling",
        "status": "unresolved",
        "reason": "Disease-risk odds ratios are multiplied as hazard multipliers and can produce extreme joint risks.",
    },
    {
        "item": "Deterministic-stochastic reconciliation",
        "status": "unresolved",
        "reason": "Deterministic compatibility expected values do not reproduce the stochastic reference mean.",
    },
]
TRAJECTORY_PRESETS = {
    "rising": {
        "label": "Rising",
        "trendRatePerYear": DEFAULT_TREND_ABS_PER_YEAR,
        "description": "infection pressure has increased toward the present",
    },
    "steady": {
        "label": "Steady",
        "trendRatePerYear": 0.0,
        "description": "infection pressure has remained approximately constant",
    },
    "falling": {
        "label": "Falling",
        "trendRatePerYear": -DEFAULT_TREND_ABS_PER_YEAR,
        "description": "infection pressure has declined toward the present",
    },
}


@dataclass(frozen=True)
class InfectionHistoryCalibration:
    trajectory: str
    trend_rate_per_year: float
    recent_window_years: float
    log_scale: float
    age_shape_gamma: float
    expected_prevalence: float
    expected_age_or: float
    recent_fraction: float
    remote_fraction: float
    age_rows: tuple[dict[str, float | str], ...]

    def as_dict(self) -> dict[str, Any]:
        return {
            "model": TRAJECTORY_MODEL,
            "trajectory": self.trajectory,
            "trendRatePerYear": self.trend_rate_per_year,
            "recentWindowYears": self.recent_window_years,
            "logScale": self.log_scale,
            "ageShapeGamma": self.age_shape_gamma,
            "expectedPrevalence": self.expected_prevalence,
            "expectedAgeOR": self.expected_age_or,
            "recentFraction": self.recent_fraction,
            "remoteFraction": self.remote_fraction,
            "ageRows": list(self.age_rows),
            "trendStatus": "scenario_assumption",
            "identifiability": (
                "Overall LTBI prevalence and an age odds ratio calibrate the scale "
                "and age-shape for a fixed trajectory, but do not identify the "
                "calendar-time trend. Rising, steady and falling slopes are scenario "
                "assumptions."
            ),
        }


def trajectory_label(value: str | None) -> str:
    key = normalise_trajectory(value)
    return TRAJECTORY_PRESETS[key]["label"]


def normalise_trajectory(value: str | None) -> str:
    key = str(value or "steady").strip().lower()
    if key not in TRAJECTORY_PRESETS:
        raise ValueError(f"Unknown infection-pressure trajectory: {value}")
    return key


def configure_infection_history_assumptions(config: dict[str, Any], trajectory: str) -> dict[str, Any]:
    key = normalise_trajectory(trajectory)
    preset = TRAJECTORY_PRESETS[key]
    out = dict(config)
    nested = dict(out.get("ltbiStateAssumptions") or {})
    nested.update(
        {
            "baselineRecentLTBIProportion": None,
            "baselineRecentLTBIDerivationMethod": DERIVATION_METHOD,
            "infectionPressureTrajectory": key,
            "infectionPressureTrajectoryLabel": preset["label"],
            "infectionPressureTrendRatePerYear": preset["trendRatePerYear"],
            "recentDefinitionYears": DEFAULT_RECENT_WINDOW_YEARS,
            "stateDefinition": (
                "Infection acquired within the preceding two years is derived "
                "from the historical infection-pressure model. The early "
                "higher-progression-risk state is separate and currently uses a "
                "five-year mean residence time before lower-risk progression."
            ),
            "baselineRecentLTBIProportionSource": (
                "Model-derived from calibrated historical infection-pressure "
                f"trajectory ({preset['label'].lower()} scenario)."
            ),
            "baselineRecentLTBIProportionStatus": "model_derived_reviewed",
            "source": (
                "Scenario assumption using repository LTBI prevalence, APY age "
                "distribution and LTBI age odds ratio; trajectory slope is not "
                "directly observed."
            ),
            "status": "model_derived_reviewed",
            "provisional": True,
            "developmentCompatibilityMode": False,
            "warnings": [
                "Infection acquired within the preceding two years is model-derived from the selected transmission-history assumption; it is not directly observed.",
                "The two-year infection-timing definition is not equivalent to the early higher-progression-risk state duration.",
            ],
        }
    )
    out["ltbiStateAssumptions"] = nested
    out["analysisBasis"] = EXPLICIT_EXPERIMENTAL_BASIS
    out["naturalHistorySemantics"] = EXPLICIT_EXPERIMENTAL_BASIS
    return out


def configure_compatibility_reference_assumptions(config: dict[str, Any]) -> dict[str, Any]:
    out = dict(config)
    nested = dict(out.get("ltbiStateAssumptions") or {})
    nested.update(
        {
            "baselineRecentLTBIProportion": None,
            "baselineRecentLTBIDerivationMethod": "",
            "infectionPressureTrajectory": "",
            "infectionPressureTrajectoryLabel": "",
            "infectionPressureTrendRatePerYear": None,
            "recentDefinitionYears": None,
            "baselineRecentLTBIProportionSource": "",
            "baselineRecentLTBIProportionStatus": "unresolved",
            "source": (
                "Transition structure from older static/transmission-dynamic "
                "architecture; APY-specific baseline recent fraction unresolved."
            ),
            "status": "unresolved_development_compatibility",
            "provisional": True,
            "developmentCompatibilityMode": True,
            "warnings": [],
        }
    )
    out["ltbiStateAssumptions"] = nested
    out.pop("baselineRecentLTBIProportion", None)
    out.pop("recentToRemoteTransitionRatePerYear", None)
    out["analysisBasis"] = COMPATIBILITY_REFERENCE_BASIS
    out["naturalHistorySemantics"] = MATLAB_V9_COMPATIBILITY_SEMANTICS
    return out


def has_explicit_infection_history(config: dict[str, Any]) -> bool:
    nested = config.get("ltbiStateAssumptions") if isinstance(config, dict) else {}
    return (
        isinstance(nested, dict)
        and nested.get("baselineRecentLTBIDerivationMethod") == DERIVATION_METHOD
    )


def is_experimental_infection_history_config(config: dict[str, Any] | None) -> bool:
    if not isinstance(config, dict):
        return False
    return (
        config.get("analysisBasis") == EXPLICIT_EXPERIMENTAL_BASIS
        or config.get("naturalHistorySemantics") == EXPLICIT_EXPERIMENTAL_BASIS
        or has_explicit_infection_history(config)
    )


def infection_history_readiness_status() -> dict[str, Any]:
    return {
        "status": "not_validated_for_decision_use",
        "label": EXPERIMENTAL_STATUS_LABEL,
        "items": [dict(item) for item in EXPERIMENTAL_READINESS_ITEMS],
    }


def infection_history_basis_from_results(results_bundle: dict[str, Any] | None) -> str:
    if not isinstance(results_bundle, dict):
        return ""
    metadata = results_bundle.get("metadata") if isinstance(results_bundle.get("metadata"), dict) else {}
    if metadata.get("analysisBasis"):
        return str(metadata.get("analysisBasis"))
    technical = results_bundle.get("technical") if isinstance(results_bundle.get("technical"), dict) else {}
    config = technical.get("interfaceConfig") if isinstance(technical.get("interfaceConfig"), dict) else {}
    return str(config.get("analysisBasis") or "")


def is_experimental_infection_history_results(results_bundle: dict[str, Any] | None) -> bool:
    if not isinstance(results_bundle, dict):
        return False
    if infection_history_basis_from_results(results_bundle) == EXPLICIT_EXPERIMENTAL_BASIS:
        return True
    technical = results_bundle.get("technical") if isinstance(results_bundle.get("technical"), dict) else {}
    config = technical.get("interfaceConfig") if isinstance(technical.get("interfaceConfig"), dict) else {}
    return is_experimental_infection_history_config(config)


def calibrate_infection_history(
    pars: dict[str, Any],
    *,
    target_prevalence: float,
    target_age_or: float,
    trajectory: str = "steady",
    trend_rate_per_year: float | None = None,
    recent_window_years: float = DEFAULT_RECENT_WINDOW_YEARS,
    age_shape_gamma: float | None = None,
) -> dict[str, Any]:
    key = normalise_trajectory(trajectory)
    trend = (
        float(TRAJECTORY_PRESETS[key]["trendRatePerYear"])
        if trend_rate_per_year is None
        else float(trend_rate_per_year)
    )
    recent_window = float(recent_window_years)
    if recent_window <= 0:
        raise ValueError("recent_window_years must be positive.")
    cache_key = (
        id(pars),
        round(float(target_prevalence), 12),
        round(float(target_age_or), 12),
        key,
        round(float(trend), 12),
        round(float(recent_window), 12),
        None if age_shape_gamma is None else round(float(age_shape_gamma), 12),
    )
    if cache_key in _CALIBRATION_CACHE:
        return dict(_CALIBRATION_CACHE[cache_key])

    if age_shape_gamma is not None:
        gamma = float(age_shape_gamma)
    else:
        def objective(log_gamma: float) -> float:
            candidate_gamma = math.exp(log_gamma)
            candidate_scale = _solve_log_scale_for_gamma(
                pars,
                gamma=candidate_gamma,
                trend_rate_per_year=trend,
                target_prevalence=target_prevalence,
            )
            expected = expected_age_or(
                pars,
                log_scale=candidate_scale,
                gamma=candidate_gamma,
                trend_rate_per_year=trend,
            )
            return math.log(expected / float(target_age_or))

        grid = np.linspace(-8.0, 3.0, 89)
        values = [_safe_eval(objective, x) for x in grid]
        log_gamma = None
        for left, right, f_left, f_right in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if math.isfinite(f_left) and math.isfinite(f_right) and np.sign(f_left) != np.sign(f_right):
                log_gamma = _bisect_root(objective, float(left), float(right))
                break
        if log_gamma is None:
            log_gamma = _minimise_scalar_grid(lambda x: objective(x) ** 2, -8.0, 3.0)
        gamma = math.exp(log_gamma)
    log_scale = _solve_log_scale_for_gamma(
        pars,
        gamma=gamma,
        trend_rate_per_year=trend,
        target_prevalence=target_prevalence,
    )
    expected_prev = expected_prevalence(
        pars,
        log_scale=log_scale,
        gamma=gamma,
        trend_rate_per_year=trend,
    )
    expected_or = expected_age_or(
        pars,
        log_scale=log_scale,
        gamma=gamma,
        trend_rate_per_year=trend,
    )
    recent = expected_recent_prevalent_fraction(
        pars,
        log_scale=log_scale,
        gamma=gamma,
        trend_rate_per_year=trend,
        recent_window_years=recent_window,
    )
    age_rows = tuple(
        _age_group_rows(
            pars,
            log_scale=log_scale,
            gamma=gamma,
            trend_rate_per_year=trend,
            recent_window_years=recent_window,
        )
    )
    result = InfectionHistoryCalibration(
        trajectory=key,
        trend_rate_per_year=trend,
        recent_window_years=recent_window,
        log_scale=log_scale,
        age_shape_gamma=gamma,
        expected_prevalence=expected_prev,
        expected_age_or=expected_or,
        recent_fraction=recent,
        remote_fraction=max(1.0 - recent, 0.0),
        age_rows=age_rows,
    ).as_dict()
    _CALIBRATION_CACHE[cache_key] = dict(result)
    return result


def _default_age_shape_gamma(target_age_or: float) -> float:
    # Used only when the caller has not supplied the APY age-shape calibration.
    # It gives a monotone age pattern without pretending that prevalence and a
    # trajectory slope identify both age shape and timing.
    return max(0.05, min(math.log(max(float(target_age_or), 1.01)) / math.log(25.5 / 5.5), 4.0))


def prevalent_infection_probabilities(
    age_years,
    pars: dict[str, Any],
    infection_history: dict[str, Any],
    marijuana,
    contact,
    renal,
) -> np.ndarray:
    log_scale = float(infection_history["logScale"])
    gamma = float(infection_history["ageShapeGamma"])
    trend = float(infection_history["trendRatePerYear"])
    cum_hazard = _cumulative_hazard(np.asarray(age_years, dtype=float), log_scale, gamma, trend)
    multiplier = _infection_multiplier(pars, marijuana, contact, renal)
    return 1.0 - np.exp(-cum_hazard * multiplier)


def conditional_recent_probability(
    age_years,
    pars: dict[str, Any],
    infection_history: dict[str, Any],
    marijuana,
    contact,
    renal,
) -> np.ndarray:
    log_scale = float(infection_history["logScale"])
    gamma = float(infection_history["ageShapeGamma"])
    trend = float(infection_history["trendRatePerYear"])
    window = float(infection_history["recentWindowYears"])
    ages = np.asarray(age_years, dtype=float)
    multiplier = _infection_multiplier(pars, marijuana, contact, renal)
    total_h = _cumulative_hazard(ages, log_scale, gamma, trend) * multiplier
    recent_h = _cumulative_hazard(ages, log_scale, gamma, trend, upper=np.minimum(ages, window)) * multiplier
    p_recent_any = 1.0 - np.exp(-recent_h)
    p_ever = 1.0 - np.exp(-total_h)
    with np.errstate(divide="ignore", invalid="ignore"):
        conditional = np.divide(p_recent_any, p_ever, out=np.zeros_like(p_ever), where=p_ever > 0)
    return np.clip(conditional, 0.0, 1.0)


def expected_prevalence(
    pars: dict[str, Any],
    *,
    log_scale: float,
    gamma: float,
    trend_rate_per_year: float,
) -> float:
    return _expected_component(pars, log_scale, gamma, trend_rate_per_year, component="ever")


def expected_age_or(
    pars: dict[str, Any],
    *,
    log_scale: float,
    gamma: float,
    trend_rate_per_year: float,
    age_cut: float = 25.0,
) -> float:
    old = _expected_component(pars, log_scale, gamma, trend_rate_per_year, component="ever", age_min=age_cut)
    young = _expected_component(pars, log_scale, gamma, trend_rate_per_year, component="ever", age_max=age_cut)
    return _odds(old) / max(_odds(young), np.finfo(float).eps)


def expected_recent_prevalent_fraction(
    pars: dict[str, Any],
    *,
    log_scale: float,
    gamma: float,
    trend_rate_per_year: float,
    recent_window_years: float,
) -> float:
    recent = _expected_component(
        pars,
        log_scale,
        gamma,
        trend_rate_per_year,
        component="recent",
        recent_window_years=recent_window_years,
    )
    ever = expected_prevalence(
        pars,
        log_scale=log_scale,
        gamma=gamma,
        trend_rate_per_year=trend_rate_per_year,
    )
    return 0.0 if ever <= 0 else float(recent / ever)


def _age_group_rows(
    pars: dict[str, Any],
    *,
    log_scale: float,
    gamma: float,
    trend_rate_per_year: float,
    recent_window_years: float,
) -> list[dict[str, float | str]]:
    rows = []
    groups = broad_age_group_from_years(pars["exactAgeValues"])
    labels = {1: "0-4 years", 2: "5-14 years", 3: "15 years and older"}
    for group in (1, 2, 3):
        mask = [g == group for g in groups]
        prevalence = _expected_component(
            pars,
            log_scale,
            gamma,
            trend_rate_per_year,
            component="ever",
            mask=mask,
        )
        recent = _expected_component(
            pars,
            log_scale,
            gamma,
            trend_rate_per_year,
            component="recent",
            recent_window_years=recent_window_years,
            mask=mask,
        )
        rows.append(
            {
                "ageGroup": labels[group],
                "ltbiPrevalence": prevalence,
                "recentFractionAmongPrevalent": 0.0 if prevalence <= 0 else recent / prevalence,
                "remoteFractionAmongPrevalent": 1.0 if prevalence <= 0 else max(1.0 - recent / prevalence, 0.0),
            }
        )
    return rows


def _expected_component(
    pars: dict[str, Any],
    log_scale: float,
    gamma: float,
    trend_rate_per_year: float,
    *,
    component: str,
    recent_window_years: float = DEFAULT_RECENT_WINDOW_YEARS,
    age_min: float | None = None,
    age_max: float | None = None,
    mask: list[bool] | None = None,
) -> float:
    ages, weights, multipliers, age_indices = _age_risk_grid(pars)
    keep = np.ones_like(ages, dtype=bool)
    if mask is not None:
        keep &= np.asarray([bool(mask[idx]) for idx in age_indices], dtype=bool)
    if age_min is not None:
        keep &= ages >= float(age_min)
    if age_max is not None:
        keep &= ages < float(age_max)
    if not keep.any():
        return 0.0
    h_total = _cumulative_hazard(ages[keep], log_scale, gamma, trend_rate_per_year)
    if component == "ever":
        values = 1.0 - np.exp(-h_total * multipliers[keep])
    else:
        h_recent = _cumulative_hazard(
            ages[keep],
            log_scale,
            gamma,
            trend_rate_per_year,
            upper=np.minimum(ages[keep], float(recent_window_years)),
        )
        values = 1.0 - np.exp(-h_recent * multipliers[keep])
    numerator = float(np.sum(weights[keep] * values))
    if mask is not None or age_min is not None or age_max is not None:
        denominator = float(np.sum(weights[keep]))
        return 0.0 if denominator <= 0 else numerator / denominator
    return numerator


def _age_risk_grid(pars: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    cache_key = id(pars)
    if cache_key in _AGE_RISK_GRID_CACHE:
        return _AGE_RISK_GRID_CACHE[cache_key]
    ages = []
    weights = []
    multipliers = []
    age_indices = []
    for age_idx, (age, age_prob) in enumerate(zip(pars["exactAgeValues"], pars["exactAgeProb"])):
        age_group = broad_age_group_from_years([age])[0] - 1
        p_mj = pars["mjPrevByAge"][age_group]
        p_contact = pars["contactPrevByAge"][age_group]
        p_renal = pars["renalPrevByAge"][age_group]
        for mj, contact, renal in product([0, 1], repeat=3):
            ages.append(float(age))
            weights.append(
                float(age_prob)
                * _bern_prob(mj, p_mj)
                * _bern_prob(contact, p_contact)
                * _bern_prob(renal, p_renal)
            )
            multipliers.append(float(_infection_multiplier(pars, [mj], [contact], [renal])[0]))
            age_indices.append(age_idx)
    grid = (
        np.asarray(ages, dtype=float),
        np.asarray(weights, dtype=float),
        np.asarray(multipliers, dtype=float),
        np.asarray(age_indices, dtype=int),
    )
    _AGE_RISK_GRID_CACHE[cache_key] = grid
    return grid


def _cumulative_hazard(
    age,
    log_scale: float,
    gamma: float,
    trend_rate_per_year: float,
    *,
    upper=None,
):
    ages = np.asarray(age, dtype=float)
    upper_arr = ages if upper is None else np.asarray(upper, dtype=float)
    upper_arr = np.maximum(np.minimum(upper_arr, ages), 0.0)
    tau = np.expand_dims(upper_arr, axis=-1) * _UNIT_NODES
    age_at_infection = np.maximum(np.expand_dims(ages, axis=-1) - tau, 0.0)
    calendar = np.exp(-float(trend_rate_per_year) * tau)
    age_effect = (age_at_infection + 0.5) ** float(gamma)
    integral = np.sum(calendar * age_effect * _UNIT_WEIGHTS, axis=-1) * upper_arr
    return math.exp(float(log_scale)) * integral


def _infection_multiplier(pars: dict[str, Any], marijuana, contact, renal) -> np.ndarray:
    inf_or = pars["infOR"]
    return (
        (float(inf_or["MJ"]) ** np.asarray(marijuana, dtype=float))
        * (float(inf_or["contact"]) ** np.asarray(contact, dtype=float))
        * (float(inf_or["renal"]) ** np.asarray(renal, dtype=float))
    )


def _solve_log_scale_for_gamma(
    pars: dict[str, Any],
    *,
    gamma: float,
    trend_rate_per_year: float,
    target_prevalence: float,
) -> float:
    def objective(log_scale: float) -> float:
        return expected_prevalence(
            pars,
            log_scale=log_scale,
            gamma=gamma,
            trend_rate_per_year=trend_rate_per_year,
        ) - float(target_prevalence)

    return _bisect_root(objective, -40.0, 10.0)


def _bisect_root(func, lo: float, hi: float, *, tol: float = 1e-11, max_iter: int = 200) -> float:
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
    return float(grid[int(np.nanargmin(values))])


def _safe_eval(func, x: float) -> float:
    try:
        return float(func(float(x)))
    except (ValueError, OverflowError, ZeroDivisionError):
        return math.nan


def _odds(p: float) -> float:
    return float(p) / max(1.0 - float(p), np.finfo(float).eps)


def _bern_prob(x: int | bool, p: float) -> float:
    return float(p) if bool(x) else 1.0 - float(p)
