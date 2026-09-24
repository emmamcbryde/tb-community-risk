"""Incidence-trend estimation interface (design only in this milestone).

The trend of *estimated TB disease incidence* is a descriptive summary of WHO or
local estimates. It is not, by itself, the trend in infection pressure or force
of infection: disease incidence also reflects progression from past infection,
risk-factor prevalence, and changes in case detection and notification. No
estimator here may be wired into the transmission or infection-history model
without an explicit, documented linkage model.

Planned methods (not yet implemented): log-linear recent trend, penalised spline
and state-space trend, each preserving observed values (with WHO uncertainty
bounds) alongside fitted values, with a selectable fitting period, COVID-era
disruption flags, diagnostic plots, annual percentage change and user override.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Protocol, Sequence

from engine.profiles.population_profile import IncidencePoint


class TrendMethod(str, Enum):
    NOT_ESTIMATED = "not_estimated"
    LOG_LINEAR_RECENT = "log_linear_recent"
    PENALISED_SPLINE = "penalised_spline"
    STATE_SPACE = "state_space"
    USER_OVERRIDE = "user_override"


IMPLEMENTED_METHODS = frozenset({TrendMethod.NOT_ESTIMATED, TrendMethod.USER_OVERRIDE})
DEFAULT_COVID_DISRUPTION_YEARS = (2020, 2021, 2022)

EPIDEMIOLOGICAL_QUANTITIES = {
    "estimated_tb_disease_incidence": (
        "Estimated new and relapse TB disease episodes per 100,000 population per year "
        "(WHO or local estimate, with uncertainty bounds)."
    ),
    "infection_pressure": "The overall rate at which susceptible people acquire new M. tuberculosis infection.",
    "force_of_infection": "Per-capita hazard of new infection; a model quantity, not directly observed.",
    "progression_to_disease": "Rate at which infected people develop TB disease, modified by time since infection and risk factors.",
    "case_detection_and_notification": "The fraction of incident disease diagnosed and reported; changes in it alter notifications, not incidence.",
}
INCIDENCE_TO_INFECTION_POLICY = (
    "The slope of estimated TB disease incidence is not used as the slope of infection "
    "pressure without an explicit, documented linkage model."
)


@dataclass(frozen=True)
class TrendSettings:
    method: TrendMethod = TrendMethod.NOT_ESTIMATED
    fit_start_year: int | None = None
    fit_end_year: int | None = None
    use_uncertainty_bounds: bool = True
    covid_disruption_years: tuple[int, ...] = DEFAULT_COVID_DISRUPTION_YEARS
    covid_handling: str = "flag"  # "flag", "exclude" or "include"
    user_annual_percent_change: float | None = None
    options: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class TrendPoint:
    year: int
    observed: float | None
    observed_lower: float | None
    observed_upper: float | None
    fitted: float | None = None
    fitted_lower: float | None = None
    fitted_upper: float | None = None
    covid_disruption: bool = False
    used_in_fit: bool = False


@dataclass(frozen=True)
class TrendResult:
    method: TrendMethod
    settings: TrendSettings
    points: tuple[TrendPoint, ...]
    annual_percent_change: float | None = None
    annual_percent_change_lower: float | None = None
    annual_percent_change_upper: float | None = None
    quantity: str = "estimated_tb_disease_incidence"
    diagnostics: dict[str, Any] = field(default_factory=dict)
    status: str = "not_estimated"


class TrendEstimator(Protocol):
    method: TrendMethod

    def fit(self, series: Sequence[IncidencePoint], settings: TrendSettings) -> TrendResult:
        """Fit a trend while preserving observed values and bounds."""


def covid_disruption_flags(years: Sequence[int], disruption_years: Sequence[int] = DEFAULT_COVID_DISRUPTION_YEARS) -> list[bool]:
    disrupted = set(disruption_years)
    return [year in disrupted for year in years]


def observed_points(series: Sequence[IncidencePoint], settings: TrendSettings) -> tuple[TrendPoint, ...]:
    flags = covid_disruption_flags([point.year for point in series], settings.covid_disruption_years)
    return tuple(
        TrendPoint(
            year=point.year,
            observed=point.estimate,
            observed_lower=point.lower,
            observed_upper=point.upper,
            covid_disruption=flag,
        )
        for point, flag in zip(sorted(series, key=lambda item: item.year), flags)
    )


def estimate_trend(series: Sequence[IncidencePoint], settings: TrendSettings) -> TrendResult:
    """Dispatch to a trend method; fitted methods are not yet implemented."""
    points = observed_points(series, settings)
    if settings.method is TrendMethod.NOT_ESTIMATED:
        return TrendResult(method=settings.method, settings=settings, points=points)
    if settings.method is TrendMethod.USER_OVERRIDE:
        if settings.user_annual_percent_change is None:
            raise ValueError("A user-override trend requires an annual percentage change.")
        return TrendResult(
            method=settings.method,
            settings=settings,
            points=points,
            annual_percent_change=float(settings.user_annual_percent_change),
            status="user_override",
        )
    raise NotImplementedError(
        f"Trend method {settings.method.value!r} is designed but not implemented in this milestone."
    )
