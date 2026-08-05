from __future__ import annotations

from typing import Any

import numpy as np


def resolve_time_settings(config: dict[str, Any]) -> dict[str, float]:
    screening = _float_alias(config, "screeningWindowYears", "screenWindow", default=3.0)
    early = _float_alias(config, "earlyProgressionPeriodYears", default=2.0)
    calibration = _float_alias(
        config,
        "activeTBCalibrationHorizonYears",
        default=2.0,
    )
    follow = _float_alias(config, "followUpHorizonYears", "followHorizon", default=20.0)
    return {
        "screeningWindowYears": screening,
        "earlyProgressionPeriodYears": early,
        "activeTBCalibrationHorizonYears": calibration,
        "followUpHorizonYears": follow,
        "legacyScreenWindow": screening,
        "legacyFollowHorizon": follow,
    }


def survival_probability(
    time,
    mult_disease,
    lambda_early: float,
    lambda_late: float,
    early_period_years: float,
):
    time_array = np.asarray(time, dtype=float)
    mult = np.asarray(mult_disease, dtype=float)
    early_time = np.minimum(np.maximum(time_array, 0.0), float(early_period_years))
    late_time = np.maximum(time_array - float(early_period_years), 0.0)
    hazard = lambda_early * mult * early_time + lambda_late * mult * late_time
    return np.exp(-hazard)


def event_probability_between(
    start: float,
    end: float,
    mult_disease,
    lambda_early: float,
    lambda_late: float,
    early_period_years: float,
):
    return np.maximum(
        survival_probability(
            start,
            mult_disease,
            lambda_early,
            lambda_late,
            early_period_years,
        )
        - survival_probability(
            end,
            mult_disease,
            lambda_early,
            lambda_late,
            early_period_years,
        ),
        0.0,
    )


def average_survival_to_uniform_screen(
    mult_disease,
    lambda_early: float,
    lambda_late: float,
    screening_window_years: float,
    early_period_years: float,
    *,
    quadrature_points: int = 32,
):
    window = float(screening_window_years)
    mult = np.asarray(mult_disease, dtype=float)
    if window <= 0:
        return np.ones_like(mult, dtype=float)
    times = (np.arange(quadrature_points, dtype=float) + 0.5) / quadrature_points * window
    values = [
        survival_probability(t, mult, lambda_early, lambda_late, early_period_years)
        for t in times
    ]
    return np.mean(values, axis=0)


def preventable_active_risk_to_uniform_screen(
    mult_disease,
    lambda_early: float,
    lambda_late: float,
    screening_window_years: float,
    follow_horizon_years: float,
    early_period_years: float,
    *,
    quadrature_points: int = 32,
):
    avg_survival = average_survival_to_uniform_screen(
        mult_disease,
        lambda_early,
        lambda_late,
        screening_window_years,
        early_period_years,
        quadrature_points=quadrature_points,
    )
    horizon_survival = survival_probability(
        follow_horizon_years,
        mult_disease,
        lambda_early,
        lambda_late,
        early_period_years,
    )
    return np.maximum(avg_survival - horizon_survival, 0.0)


def _float_alias(config: dict[str, Any], primary: str, legacy: str | None = None, *, default: float) -> float:
    value = config.get(primary)
    if value is None and legacy is not None:
        value = config.get(legacy)
    if value is None:
        value = default
    return float(value)
