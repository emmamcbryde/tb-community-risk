"""Descriptive trend estimation for estimated TB disease incidence.

The trend of estimated TB disease incidence summarises WHO or local estimates. It
is *not* the trend in infection pressure or force of infection: disease incidence
also reflects progression from past infection, risk-factor prevalence and case
detection. Nothing here feeds the transmission, infection-history or dynamic
model.

Methods
-------
Log-linear (primary)
    log(I_t) = a + b (t - t_mean) + e_t, ordinary least squares on WHO point
    estimates in the fitting window. Annual percentage change (APC) =
    100 (exp(b) - 1). The regression interval uses the t distribution of b.
Segmented disruption (log-linear option)
    log(I_t) = a + b (t - t_mean) + c D_t + e_t with D_t = 1 in nominated
    disruption years: a temporary level shift, not a causal effect estimate. APC
    is derived from b.
Penalised smooth (sensitivity)
    Whittaker-Henderson smoother on the annual grid (a discrete P-spline):
    minimise sum w_t (y_t - s_t)^2 + lambda sum (second difference of s)^2 on
    y = log(I). lambda is chosen by generalised cross-validation on a fixed grid,
    never below a conservative floor. The recent slope is the mean annual change
    of s over the last three years of the window.

Uncertainty
-----------
WHO bounds are treated as the 2.5% and 97.5% points of a split-normal
distribution on the log scale (sigma_low = (log I - log lo)/1.96, sigma_high =
(log hi - log I)/1.96). Draws are independent between years by default, which is
conservative for slope uncertainty; a common-shift option treats years as fully
correlated. Each draw is refitted with the same method and settings, with a fixed
recorded seed. The result is a *propagated uncertainty interval* for the APC. It
covers only WHO estimate uncertainty as expressed by the published bounds; it
does not include WHO methodological uncertainty, structural change, model
structural uncertainty, infection-pressure uncertainty or other parameters.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from enum import Enum
import math
from typing import Any, Sequence

import numpy as np

from engine.profiles.population_profile import IncidencePoint


class TrendMethod(str, Enum):
    NOT_ESTIMATED = "not_estimated"
    LOG_LINEAR_RECENT = "log_linear_recent"
    PENALISED_SPLINE = "penalised_spline"
    STATE_SPACE = "state_space"
    USER_OVERRIDE = "user_override"


class CovidHandling(str, Enum):
    INCLUDE = "include"
    EXCLUDE = "exclude"
    SEGMENTED = "segmented"


IMPLEMENTED_METHODS = frozenset({TrendMethod.NOT_ESTIMATED, TrendMethod.LOG_LINEAR_RECENT, TrendMethod.PENALISED_SPLINE, TrendMethod.USER_OVERRIDE})
DEFAULT_COVID_DISRUPTION_YEARS = (2020, 2021, 2022)
DEFAULT_WINDOW_YEARS = 10
WINDOW_CHOICES = {"5 years": 5, "10 years": 10, "15 years": 15, "Full series": None}
DEFAULT_SEED = 20250101
DEFAULT_DRAWS = 1000
MIN_POINTS_LOG_LINEAR = 5
MIN_POINTS_SPLINE = 8
MIN_POINTS_ABSOLUTE = 3
SPLINE_LAMBDA_GRID = tuple(10.0 ** power for power in np.arange(0.0, 6.01, 0.25))
SPLINE_LAMBDA_FLOOR = 10.0
SPLINE_SLOPE_YEARS = 3
NEAR_ZERO_RATE = 0.1
CURVATURE_P_THRESHOLD = 0.05
RESIDUAL_SD_THRESHOLD = 0.10
DISRUPTION_APC_DIFFERENCE = 1.0
DISCONTINUITY_LOG_JUMP = math.log(1.5)
NO_CLEAR_CHANGE_BAND = 1.0
Z_975 = 1.959963984540054

CLASSIFICATION_LABELS = {
    "adequate_log_linear": "Adequate log-linear fit",
    "non_linear": "Non-linear trend - review smooth estimate",
    "short_series": "Unstable because of short series",
    "disruption_sensitive": "Disruption-sensitive",
    "bounds_incomplete": "Bounds incomplete",
    "series_discontinuity": "Possible series discontinuity",
    "unsuitable": "Unsuitable for automated trend inference",
    "smooth_boundary": "Smooth estimate - boundary slope",
    "user_override": "User-specified trend",
    "not_estimated": "Trend not estimated",
}
SUMMARY_LABELS = {
    "increasing": "Estimated disease incidence increasing",
    "decreasing": "Estimated disease incidence decreasing",
    "no_clear_change": "No clear recent change",
    "uncertain": "Trend uncertain",
    "not_estimated": "Trend not estimated",
}
NO_CLEAR_CHANGE_NOTE = "No clear recent change is not evidence that the epidemic is constant."

EPIDEMIOLOGICAL_QUANTITIES = {
    "estimated_tb_disease_incidence": "Estimated new and relapse TB disease episodes per 100,000 population per year (WHO or local estimate).",
    "notification_rate": "Diagnosed and reported TB cases per 100,000; depends on case detection as well as disease occurrence.",
    "active_tb_prevalence": "People with TB disease at a point in time; depends on incidence and disease duration.",
    "force_of_infection": "Per-capita hazard of new M. tuberculosis infection; a model quantity, not directly observed.",
    "annual_risk_of_infection": "Probability of infection within a year; related to force of infection.",
    "recent_infection": "Infection acquired within the last few years, with a higher progression hazard.",
    "remote_infection": "Longer-standing infection with a lower progression hazard.",
    "progression_to_disease": "Rate at which infected people develop TB disease, modified by time since infection and risk factors.",
    "detection_and_treatment": "Diagnosis, notification and treatment; changes alter notifications and infectious duration, not incidence directly.",
}
INCIDENCE_TO_INFECTION_POLICY = (
    "The slope of estimated TB disease incidence is not used as the slope of infection "
    "pressure without an explicit, documented linkage model."
)
DESCRIPTIVE_STATEMENT = (
    "Country incidence data currently describe the TB disease burden and trend. "
    "They are not yet used to infer infection pressure or transmission."
)


@dataclass(frozen=True)
class TrendSettings:
    method: TrendMethod = TrendMethod.LOG_LINEAR_RECENT
    window_years: int | None = DEFAULT_WINDOW_YEARS
    start_year: int | None = None
    end_year: int | None = None
    covid_handling: CovidHandling = CovidHandling.INCLUDE
    disruption_years: tuple[int, ...] = DEFAULT_COVID_DISRUPTION_YEARS
    excluded_years: tuple[int, ...] = ()
    propagate_uncertainty: bool = True
    draws: int = DEFAULT_DRAWS
    seed: int = DEFAULT_SEED
    year_correlation: str = "independent"
    user_annual_percent_change: float | None = None

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["method"] = self.method.value
        payload["covid_handling"] = self.covid_handling.value
        payload["disruption_years"] = list(self.disruption_years)
        payload["excluded_years"] = list(self.excluded_years)
        return payload

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "TrendSettings":
        data = dict(payload)
        data["method"] = TrendMethod(data.get("method", TrendMethod.LOG_LINEAR_RECENT.value))
        data["covid_handling"] = CovidHandling(data.get("covid_handling", CovidHandling.INCLUDE.value))
        data["disruption_years"] = tuple(int(y) for y in data.get("disruption_years", DEFAULT_COVID_DISRUPTION_YEARS))
        data["excluded_years"] = tuple(int(y) for y in data.get("excluded_years", ()))
        allowed = {f for f in cls.__dataclass_fields__}
        return cls(**{key: value for key, value in data.items() if key in allowed})


@dataclass(frozen=True)
class TrendPoint:
    year: int
    observed: float | None
    observed_lower: float | None
    observed_upper: float | None
    fitted: float | None = None
    used_in_fit: bool = False
    in_window: bool = False
    exclusion_reason: str = ""
    covid_disruption: bool = False


@dataclass(frozen=True)
class TrendResult:
    method: TrendMethod
    settings: TrendSettings
    points: tuple[TrendPoint, ...]
    period: tuple[int, int] | None = None
    years_used: int = 0
    annual_percent_change: float | None = None
    fit_interval: tuple[float, float] | None = None
    propagated_interval: tuple[float, float] | None = None
    classification: str = "not_estimated"
    flags: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    summary: str = "not_estimated"
    diagnostics: dict[str, Any] = field(default_factory=dict)
    quantity: str = "estimated_tb_disease_incidence"
    status: str = "not_estimated"

    @property
    def classification_label(self) -> str:
        return CLASSIFICATION_LABELS[self.classification]

    @property
    def summary_label(self) -> str:
        return SUMMARY_LABELS[self.summary]

    def to_dict(self) -> dict[str, Any]:
        return {
            "quantity": self.quantity,
            "method": self.method.value,
            "settings": self.settings.to_dict(),
            "period": list(self.period) if self.period else None,
            "yearsUsed": self.years_used,
            "annualPercentChange": self.annual_percent_change,
            "fitInterval": list(self.fit_interval) if self.fit_interval else None,
            "propagatedInterval": list(self.propagated_interval) if self.propagated_interval else None,
            "classification": self.classification,
            "classificationLabel": self.classification_label,
            "flags": list(self.flags),
            "warnings": list(self.warnings),
            "summary": self.summary,
            "summaryLabel": self.summary_label,
            "diagnostics": _jsonable(self.diagnostics),
            "status": self.status,
            "points": [asdict(point) for point in self.points],
        }


def covid_disruption_flags(years: Sequence[int], disruption_years: Sequence[int] = DEFAULT_COVID_DISRUPTION_YEARS) -> list[bool]:
    disrupted = set(disruption_years)
    return [year in disrupted for year in years]


def fitting_period(series: Sequence[IncidencePoint], settings: TrendSettings) -> tuple[int, int] | None:
    years = sorted(point.year for point in series if point.estimate is not None)
    if not years:
        return None
    end = settings.end_year if settings.end_year is not None else years[-1]
    if settings.start_year is not None:
        start = settings.start_year
    elif settings.window_years is None:
        start = years[0]
    else:
        start = end - settings.window_years + 1
    return (max(start, years[0]), min(end, years[-1]))


def estimate_trend(series: Sequence[IncidencePoint], settings: TrendSettings | None = None) -> TrendResult:
    settings = settings or TrendSettings()
    ordered = sorted(series, key=lambda point: point.year)
    if settings.method is TrendMethod.STATE_SPACE:
        raise NotImplementedError("The state-space trend method is designed but not implemented.")
    if settings.method is TrendMethod.NOT_ESTIMATED:
        return TrendResult(method=settings.method, settings=settings, points=_points(ordered, settings, None, set()))
    if settings.method is TrendMethod.USER_OVERRIDE:
        if settings.user_annual_percent_change is None:
            raise ValueError("A user-specified trend requires an annual percentage change.")
        return TrendResult(
            method=settings.method,
            settings=settings,
            points=_points(ordered, settings, None, set()),
            annual_percent_change=float(settings.user_annual_percent_change),
            classification="user_override",
            summary=_summary(float(settings.user_annual_percent_change), None),
            status="user_override",
            warnings=("User-specified annual percentage change; not estimated from data.",),
        )
    if settings.method is TrendMethod.PENALISED_SPLINE and settings.covid_handling is CovidHandling.SEGMENTED:
        raise ValueError("Segmented disruption adjustment is available for the log-linear method only.")

    period = fitting_period(ordered, settings)
    if period is None:
        return _unsuitable(ordered, settings, None, ["No incidence estimates are available."])
    start, end = period
    window = [point for point in ordered if start <= point.year <= end and point.estimate is not None]
    excluded = _excluded_years(settings)
    used = [point for point in window if point.year not in excluded]
    missing_years = sorted(set(range(start, end + 1)) - {point.year for point in window})
    flags: list[str] = []
    warnings: list[str] = []
    if missing_years:
        warnings.append(f"No estimates for {missing_years} inside the fitting period; these years are not imputed.")
    if any(point.estimate <= 0 for point in used):
        return _unsuitable(ordered, settings, period, ["Zero incidence estimates cannot be analysed on the log scale."], used, excluded)
    minimum = MIN_POINTS_SPLINE if settings.method is TrendMethod.PENALISED_SPLINE else MIN_POINTS_LOG_LINEAR
    if len(used) < MIN_POINTS_ABSOLUTE:
        return _unsuitable(ordered, settings, period, [f"Only {len(used)} year(s) available for fitting."], used, excluded)
    if len(used) < minimum:
        flags.append("short_series")
        warnings.append(f"Only {len(used)} years used; at least {minimum} are recommended for this method.")
    if any(point.estimate < NEAR_ZERO_RATE for point in used):
        warnings.append("Some estimates are near zero; log-scale changes are unstable.")

    years = np.array([point.year for point in used], dtype=float)
    y = np.log(np.array([point.estimate for point in used], dtype=float))
    fit = _fit(years, y, settings, start, end)
    if fit is None:
        return _unsuitable(ordered, settings, period, ["The trend model could not be fitted to these data."], used, excluded)
    apc = _apc(fit["slope"])
    diagnostics = dict(fit["diagnostics"])
    fit_interval = None
    if fit.get("slope_se") is not None and fit["df"] > 0:
        t_crit = _t_quantile(0.975, fit["df"])
        fit_interval = (_apc(fit["slope"] - t_crit * fit["slope_se"]), _apc(fit["slope"] + t_crit * fit["slope_se"]))

    propagated = None
    bounds_ok = all(point.lower is not None and point.upper is not None and point.lower > 0 for point in used)
    if not bounds_ok:
        flags.append("bounds_incomplete")
        warnings.append("Uncertainty bounds are missing or zero for some years; propagated uncertainty is unavailable.")
    elif any(point.lower > point.estimate or point.upper < point.estimate for point in used):
        return _unsuitable(ordered, settings, period, ["Uncertainty bounds are inconsistent with the estimates."], used, excluded)
    elif settings.propagate_uncertainty and settings.draws > 0:
        propagated, draw_diag = _propagate(used, years, settings, start, end)
        diagnostics.update(draw_diag)

    if settings.method is TrendMethod.LOG_LINEAR_RECENT:
        if fit["diagnostics"].get("curvature_p") is not None and fit["diagnostics"]["curvature_p"] < CURVATURE_P_THRESHOLD and len(used) >= 6:
            flags.append("non_linear")
            warnings.append("The log-linear description is poor (significant curvature); review the smooth estimate.")
        elif fit["diagnostics"].get("residual_sd", 0.0) > RESIDUAL_SD_THRESHOLD:
            flags.append("non_linear")
            warnings.append("Residual scatter around the log-linear trend is large.")
    else:
        flags.append("smooth_boundary")
        warnings.append("The smooth slope is estimated at the end of the series, where smoothing is least stable.")

    disrupted_in_window = [year for year in settings.disruption_years if start <= year <= end]
    if disrupted_in_window and settings.method is TrendMethod.LOG_LINEAR_RECENT:
        comparison = _disruption_sensitivity([p for p in window if p.year not in set(settings.excluded_years)], settings)
        diagnostics["disruptionSensitivity"] = comparison
        if comparison is not None and abs(comparison["apcDifference"]) > DISRUPTION_APC_DIFFERENCE:
            flags.append("disruption_sensitive")
            warnings.append(
                f"The fitting period includes {disrupted_in_window} and the estimated change differs by "
                f"{comparison['apcDifference']:.1f} percentage points depending on how those years are handled."
            )
    jumps = _discontinuities(window)
    if jumps:
        flags.append("series_discontinuity")
        warnings.append(f"Large year-to-year jumps at {jumps}; a method or data change cannot be ruled out.")
    diagnostics["missingYears"] = missing_years
    diagnostics["excludedYears"] = sorted(year for year in excluded if start <= year <= end)

    classification = _classify(flags)
    interval = propagated or fit_interval
    fitted_by_year = {int(year): float(np.exp(value)) for year, value in zip(fit["grid_years"], fit["grid_fitted"])}
    return TrendResult(
        method=settings.method,
        settings=settings,
        points=_points(ordered, settings, period, excluded, fitted_by_year, {point.year for point in used}),
        period=period,
        years_used=len(used),
        annual_percent_change=apc,
        fit_interval=fit_interval,
        propagated_interval=propagated,
        classification=classification,
        flags=tuple(dict.fromkeys(flags)),
        warnings=tuple(warnings),
        summary=_summary(apc, interval) if classification != "unsuitable" else "not_estimated",
        diagnostics=diagnostics,
        status="estimated",
    )


def _fit(years: np.ndarray, y: np.ndarray, settings: TrendSettings, start: int, end: int) -> dict[str, Any] | None:
    if settings.method is TrendMethod.PENALISED_SPLINE:
        return _fit_spline(years, y, start, end, settings)
    segmented = settings.covid_handling is CovidHandling.SEGMENTED
    return _fit_log_linear(years, y, settings.disruption_years if segmented else ())


def _fit_log_linear(years: np.ndarray, y: np.ndarray, disruption_years: Sequence[int]) -> dict[str, Any] | None:
    centre = years.mean()
    columns = [np.ones_like(years), years - centre]
    dummy = np.array([1.0 if int(year) in set(disruption_years) else 0.0 for year in years])
    has_dummy = bool(disruption_years) and 0 < dummy.sum() < len(years)
    if has_dummy:
        columns.append(dummy)
    X = np.column_stack(columns)
    n, k = X.shape
    if n <= k - 1 or np.linalg.matrix_rank(X) < k:
        return None
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    residuals = y - X @ beta
    df = n - k
    sigma2 = float(residuals @ residuals / df) if df > 0 else float("nan")
    cov = sigma2 * np.linalg.inv(X.T @ X) if df > 0 else None
    total = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float(residuals @ residuals) / total if total > 0 else 1.0
    diagnostics: dict[str, Any] = {
        "model": "segmented log-linear" if has_dummy else "log-linear",
        "rSquared": r2,
        "residual_sd": math.sqrt(sigma2) if df > 0 else None,
        "residuals": [float(value) for value in residuals],
        "durbinWatson": float(np.sum(np.diff(residuals) ** 2) / (residuals @ residuals)) if float(residuals @ residuals) > 0 else None,
        "maxAbsResidual": float(np.max(np.abs(residuals))),
        "curvature_p": _curvature_p(years - centre, y) if not has_dummy else None,
    }
    if has_dummy:
        diagnostics["disruptionLevelShift"] = float(beta[2])
        diagnostics["disruptionLevelShiftPercent"] = 100.0 * (math.exp(float(beta[2])) - 1.0)
        diagnostics["disruptionNote"] = "Temporary level shift during nominated years; not a causal effect estimate."
    grid = np.arange(years.min(), years.max() + 1)
    grid_fitted = beta[0] + beta[1] * (grid - centre)
    return {
        "slope": float(beta[1]),
        "slope_se": math.sqrt(float(cov[1, 1])) if cov is not None else None,
        "df": df,
        "diagnostics": diagnostics,
        "grid_years": grid,
        "grid_fitted": grid_fitted,
    }


def _curvature_p(centred: np.ndarray, y: np.ndarray) -> float | None:
    n = len(y)
    if n < 5:
        return None
    X1 = np.column_stack([np.ones(n), centred])
    X2 = np.column_stack([np.ones(n), centred, centred**2])
    r1 = y - X1 @ np.linalg.lstsq(X1, y, rcond=None)[0]
    r2 = y - X2 @ np.linalg.lstsq(X2, y, rcond=None)[0]
    rss1, rss2 = float(r1 @ r1), float(r2 @ r2)
    df2 = n - 3
    if rss2 <= 1e-15:
        return 0.0 if rss1 > 1e-12 else 1.0
    f_stat = (rss1 - rss2) / (rss2 / df2)
    return _f_sf(f_stat, 1, df2)


def _whittaker(y_grid: np.ndarray, weights: np.ndarray, lam: float) -> tuple[np.ndarray, float]:
    n = len(y_grid)
    D = np.diff(np.eye(n), n=2, axis=0)
    W = np.diag(weights)
    A = W + lam * D.T @ D
    smooth = np.linalg.solve(A, W @ y_grid)
    hat_trace = float(np.trace(np.linalg.solve(A, W)))
    return smooth, hat_trace


def _fit_spline(years: np.ndarray, y: np.ndarray, start: int, end: int, settings: TrendSettings, fixed_lambda: float | None = None) -> dict[str, Any] | None:
    grid = np.arange(int(years.min()), int(years.max()) + 1)
    if len(grid) < 3:
        return None
    y_grid = np.zeros(len(grid))
    weights = np.zeros(len(grid))
    for year, value in zip(years, y):
        index = int(year) - int(grid[0])
        y_grid[index] = value
        weights[index] = 1.0
    n_obs = weights.sum()
    if fixed_lambda is None:
        best = None
        for lam in SPLINE_LAMBDA_GRID:
            if lam < SPLINE_LAMBDA_FLOOR:
                continue
            smooth, edf = _whittaker(y_grid, weights, lam)
            rss = float(np.sum(weights * (y_grid - smooth) ** 2))
            denominator = (1.0 - edf / n_obs) ** 2
            gcv = (rss / n_obs) / denominator if denominator > 1e-12 else math.inf
            if best is None or gcv < best[0] - 1e-15:
                best = (gcv, lam)
        lam = best[1]
    else:
        lam = fixed_lambda
    smooth, edf = _whittaker(y_grid, weights, lam)
    span = min(SPLINE_SLOPE_YEARS, len(grid) - 1)
    slope = float((smooth[-1] - smooth[-1 - span]) / span)
    residuals = [float(value) for value, weight in zip(y_grid - smooth, weights) if weight > 0]
    return {
        "slope": slope,
        "slope_se": None,
        "df": 0,
        "lambda": lam,
        "diagnostics": {
            "model": "penalised smooth (Whittaker-Henderson, second differences)",
            "smoothingParameter": lam,
            "smoothingSelection": f"generalised cross-validation on 10^0..10^6, floor {SPLINE_LAMBDA_FLOOR:g}",
            "effectiveDegreesOfFreedom": edf,
            "slopeDefinition": f"mean annual change of the smooth over the last {span} year(s)",
            "residuals": residuals,
            "residual_sd": float(np.std(residuals, ddof=1)) if len(residuals) > 1 else None,
        },
        "grid_years": grid,
        "grid_fitted": smooth,
    }


def _propagate(used: Sequence[IncidencePoint], years: np.ndarray, settings: TrendSettings, start: int, end: int) -> tuple[tuple[float, float] | None, dict[str, Any]]:
    rng = np.random.default_rng(settings.seed)
    m = np.log(np.array([point.estimate for point in used]))
    sigma_low = (m - np.log(np.array([point.lower for point in used]))) / Z_975
    sigma_high = (np.log(np.array([point.upper for point in used])) - m) / Z_975
    if settings.year_correlation == "common_shift":
        z = np.repeat(rng.standard_normal((settings.draws, 1)), len(used), axis=1)
    else:
        z = rng.standard_normal((settings.draws, len(used)))
    draws = m + np.where(z > 0, sigma_high, sigma_low) * z
    slopes = _vectorised_slopes(years, m, draws, settings, start, end)
    if slopes is None or len(slopes) == 0:
        return None, {}
    apcs = 100.0 * (np.exp(slopes) - 1.0)
    lower, upper = np.percentile(apcs, [2.5, 97.5])
    return (float(lower), float(upper)), {
        "propagation": {
            "distribution": "split-normal on log incidence; bounds as 2.5% and 97.5% points",
            "yearCorrelation": settings.year_correlation,
            "draws": settings.draws,
            "seed": settings.seed,
            "asymmetry": [float(high / low) if low > 0 else None for low, high in zip(sigma_low, sigma_high)],
            "interval": "propagated uncertainty interval (2.5th-97.5th percentile of refitted APC)",
        }
    }


def _vectorised_slopes(years: np.ndarray, m: np.ndarray, draws: np.ndarray, settings: TrendSettings, start: int, end: int) -> np.ndarray | None:
    """Slopes for all draws at once; both estimators are linear in log incidence."""
    if settings.method is TrendMethod.PENALISED_SPLINE:
        base = _fit_spline(years, m, start, end, settings)
        if base is None:
            return None
        grid = np.arange(int(years.min()), int(years.max()) + 1)
        weights = np.zeros(len(grid))
        index = (years - grid[0]).astype(int)
        weights[index] = 1.0
        n = len(grid)
        D = np.diff(np.eye(n), n=2, axis=0)
        A = np.diag(weights) + base["lambda"] * D.T @ D
        y_grid = np.zeros((draws.shape[0], n))
        y_grid[:, index] = draws
        smooth = np.linalg.solve(A, (weights * y_grid).T).T
        span = min(SPLINE_SLOPE_YEARS, n - 1)
        return (smooth[:, -1] - smooth[:, -1 - span]) / span
    centre = years.mean()
    columns = [np.ones_like(years), years - centre]
    if settings.covid_handling is CovidHandling.SEGMENTED:
        dummy = np.array([1.0 if int(year) in set(settings.disruption_years) else 0.0 for year in years])
        if 0 < dummy.sum() < len(years):
            columns.append(dummy)
    X = np.column_stack(columns)
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return None
    return (np.linalg.pinv(X) @ draws.T)[1]


def _disruption_sensitivity(window: Sequence[IncidencePoint], settings: TrendSettings) -> dict[str, Any] | None:
    """Log-linear APC with and without the disruption years (diagnostic only)."""
    disrupted = set(settings.disruption_years)
    results = {}
    for label, points in (("apcIncluded", list(window)), ("apcExcluded", [p for p in window if p.year not in disrupted])):
        if len(points) < MIN_POINTS_ABSOLUTE or any(p.estimate <= 0 for p in points):
            return None
        fit = _fit_log_linear(np.array([p.year for p in points], dtype=float), np.log(np.array([p.estimate for p in points])), ())
        if fit is None:
            return None
        results[label] = _apc(fit["slope"])
    results["apcDifference"] = results["apcIncluded"] - results["apcExcluded"]
    return results


def _discontinuities(window: Sequence[IncidencePoint]) -> list[int]:
    jumps = []
    for previous, current in zip(window, window[1:]):
        if current.year == previous.year + 1 and previous.estimate and current.estimate and previous.estimate > 0 and current.estimate > 0:
            if abs(math.log(current.estimate / previous.estimate)) > DISCONTINUITY_LOG_JUMP:
                jumps.append(current.year)
    return jumps


def _classify(flags: Sequence[str]) -> str:
    for code in ("unsuitable", "short_series", "series_discontinuity", "disruption_sensitive", "non_linear", "bounds_incomplete", "smooth_boundary"):
        if code in flags:
            return code
    return "adequate_log_linear"


def _summary(apc: float, interval: tuple[float, float] | None) -> str:
    if interval is None:
        return "uncertain"
    lower, upper = interval
    if lower > 0:
        return "increasing"
    if upper < 0:
        return "decreasing"
    if -NO_CLEAR_CHANGE_BAND <= lower and upper <= NO_CLEAR_CHANGE_BAND:
        return "no_clear_change"
    return "uncertain"


def _excluded_years(settings: TrendSettings) -> set[int]:
    excluded = set(settings.excluded_years)
    if settings.covid_handling is CovidHandling.EXCLUDE:
        excluded |= set(settings.disruption_years)
    return excluded


def _points(
    ordered: Sequence[IncidencePoint],
    settings: TrendSettings,
    period: tuple[int, int] | None,
    excluded: set[int],
    fitted: dict[int, float] | None = None,
    used_years: set[int] | None = None,
) -> tuple[TrendPoint, ...]:
    fitted = fitted or {}
    used_years = used_years or set()
    disrupted = set(settings.disruption_years)
    out = []
    for point in ordered:
        in_window = period is not None and period[0] <= point.year <= period[1]
        reason = ""
        if in_window and point.year in excluded:
            reason = "excluded by trend settings"
        out.append(
            TrendPoint(
                year=point.year,
                observed=point.estimate,
                observed_lower=point.lower,
                observed_upper=point.upper,
                fitted=fitted.get(point.year),
                used_in_fit=point.year in used_years,
                in_window=in_window,
                exclusion_reason=reason,
                covid_disruption=point.year in disrupted,
            )
        )
    return tuple(out)


def _unsuitable(ordered, settings, period, warnings, used=(), excluded=None) -> TrendResult:
    return TrendResult(
        method=settings.method,
        settings=settings,
        points=_points(ordered, settings, period, excluded or set(), None, {p.year for p in used}),
        period=period,
        years_used=len(used),
        classification="unsuitable",
        flags=("unsuitable",),
        warnings=tuple(warnings),
        summary="not_estimated",
        status="unsuitable",
    )


def _apc(slope: float) -> float:
    return 100.0 * (math.exp(slope) - 1.0)


def _t_quantile(p: float, df: int) -> float:
    from scipy import stats

    return float(stats.t.ppf(p, df))


def _f_sf(value: float, d1: int, d2: int) -> float:
    from scipy import stats

    return float(stats.f.sf(value, d1, d2))


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def settings_with(settings: TrendSettings, **changes: Any) -> TrendSettings:
    return replace(settings, **changes)
