"""Pure mathematics for exogenous background TB infection pressure (catalytic hazard).

Status: isolated and validated, **not connected** to cohort generation, the
epidemiological engine, the event ledger, economics or any page. Only tests import
it. The specification is ``docs/catalytic_infection_pressure_spec.md``.

Conventions (section 3a of the specification):

* The hazard λ(t) is in infections per susceptible person-year. It is a rate, not a
  probability; values above 1 are mathematically valid.
* Time is calendar time in decimal years. An annual series value for year ``y``
  applies on the half-open interval ``[y, y + 1)``.
* Intervals are half-open ``[t0, t1)`` with ``t0 <= t1``. A zero-duration interval
  accumulates no hazard and needs no coverage.
* A time series must cover every year that the interval touches. Nothing is carried
  forward, backfilled, interpolated or replaced with zero.
* ``H(t0, t1) = ∫ λ(s) ds`` and ``P = 1 - exp(-H)``, computed as ``-expm1(-H)``.
* First infection: the caller supplies ``U`` in the open interval (0, 1), or the
  exponential threshold ``E = -log(U)`` directly (finite, ``E >= 0``). Infection
  occurs in ``[t0, t1)`` if and only if ``E < H(t0, t1)``, at
  ``τ = inf{s : H(t0, s) > E}``. Otherwise the result carries the residual threshold
  ``E - H(t0, t1)`` for the next interval. No function here draws random numbers.
* Only first infection of people susceptible at ``t0`` is covered. Reinfection,
  progression after a new infection and baseline LTBI prevalence are out of scope:
  baseline infection status remains a separate input.
* Age-specific hazards are specified in ``background_exposure_v1`` but are not yet
  executable; building a schedule from such a configuration raises
  :class:`AgeSpecificHazardNotExecutableError`.
* Mode ``none`` has λ = 0 exactly: H = 0, P = 0, no infection, and
  :func:`requires_random_draw` is False so callers draw nothing.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
import math
from typing import Union

import numpy as np

from engine.profiles.background_exposure import (
    PLAUSIBILITY_WARNING_HAZARD,
    BackgroundExposure,
    BackgroundExposureValidationError,
    ExposureMode,
    hazard_plausibility_warning,
)
from engine.profiles.population_profile import ValueState


class HazardInputError(BackgroundExposureValidationError):
    """An invalid hazard, time, interval or random variate."""


class HazardCoverageError(BackgroundExposureValidationError):
    """An annual series does not cover the requested interval."""


class AgeSpecificHazardNotExecutableError(BackgroundExposureValidationError):
    """Age-specific hazards are specified in the schema but not yet executable."""


@dataclass(frozen=True)
class HazardSchedule:
    """Executable hazard λ(t). Build with :meth:`none`, :meth:`constant` or :meth:`annual_series`."""

    mode: ExposureMode
    constant_hazard: float = 0.0
    series: tuple[tuple[int, float], ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.mode, ExposureMode):
            raise HazardInputError("mode must be an ExposureMode.")
        if self.mode is ExposureMode.NONE:
            if self.constant_hazard != 0.0 or self.series:
                raise HazardInputError("Mode 'none' takes no hazard values.")
        elif self.mode is ExposureMode.CONSTANT:
            _check_hazard(self.constant_hazard, "constant hazard")
            if self.series:
                raise HazardInputError("Mode 'constant' takes no time series.")
        else:
            if self.constant_hazard != 0.0:
                raise HazardInputError("Mode 'time_series' takes no constant hazard.")
            if not self.series:
                raise HazardInputError("Mode 'time_series' needs at least one year.")
            years = [year for year, _ in self.series]
            for year, hazard in self.series:
                _check_year(year)
                _check_hazard(hazard, f"year {year}")
            duplicates = sorted({year for year in years if years.count(year) > 1})
            if duplicates:
                raise HazardInputError(f"Duplicate years in time series: {duplicates}.")
            if years != sorted(years):
                raise HazardInputError("Time-series years must be in increasing order.")

    @classmethod
    def none(cls) -> "HazardSchedule":
        return cls(ExposureMode.NONE)

    @classmethod
    def constant(cls, hazard: float) -> "HazardSchedule":
        return cls(ExposureMode.CONSTANT, constant_hazard=_as_float(hazard, "constant hazard"))

    @classmethod
    def annual_series(cls, values: Mapping[int, float] | Iterable[tuple[int, float]]) -> "HazardSchedule":
        """Annual values; year ``y`` applies on ``[y, y + 1)``. Pairs may be unsorted; duplicates are rejected."""
        pairs = list(values.items()) if isinstance(values, Mapping) else list(values)
        rows = []
        for pair in pairs:
            if not isinstance(pair, tuple) or len(pair) != 2:
                raise HazardInputError("Time-series entries must be (year, hazard) pairs.")
            year, hazard = pair
            _check_year(year)
            rows.append((year, _as_float(hazard, f"year {year}")))
        years = [year for year, _ in rows]
        duplicates = sorted({year for year in years if years.count(year) > 1})
        if duplicates:
            raise HazardInputError(f"Duplicate years in time series: {duplicates}.")
        return cls(ExposureMode.TIME_SERIES, series=tuple(sorted(rows)))

    @property
    def years(self) -> tuple[int, ...]:
        return tuple(year for year, _ in self.series)


ExponentialThreshold = float


@dataclass(frozen=True)
class InfectionOccurs:
    """First infection at ``time`` within ``[t0, t1)``."""

    time: float
    infected: bool = True


@dataclass(frozen=True)
class NoInfectionInInterval:
    """No infection within ``[t0, t1)``.

    ``residual_threshold`` is ``E - H(t0, t1)``: the exponential threshold to carry into
    the next interval for the same person. It is ``None`` for mode ``none``, where no
    variate is used.
    """

    cumulative_hazard: float
    residual_threshold: float | None
    infected: bool = False


FirstInfectionResult = Union[InfectionOccurs, NoInfectionInInterval]


# --- validation helpers ---------------------------------------------------------------


def _as_float(value: object, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise HazardInputError(f"{context}: hazard must be a number.")
    return float(value)


def _check_hazard(value: float, context: str) -> None:
    if not isinstance(value, float) or not math.isfinite(value):
        raise HazardInputError(f"{context}: hazard must be a finite float.")
    if value < 0.0:
        raise HazardInputError(f"{context}: hazard cannot be negative.")


def _check_year(year: object) -> None:
    if isinstance(year, bool) or not isinstance(year, (int, np.integer)):
        raise HazardInputError(f"Time-series year {year!r} must be an integer.")


def _check_time(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise HazardInputError(f"{name} must be a number.")
    value = float(value)
    if not math.isfinite(value):
        raise HazardInputError(f"{name} must be finite.")
    return value


def _check_interval(t0: object, t1: object) -> tuple[float, float]:
    start, stop = _check_time(t0, "t0"), _check_time(t1, "t1")
    if stop < start:
        raise HazardInputError(f"Reversed interval: t1 ({stop:g}) is before t0 ({start:g}).")
    return start, stop


def _check_threshold(threshold: object) -> float:
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float, np.integer, np.floating)):
        raise HazardInputError("The exponential threshold must be a number.")
    value = float(threshold)
    if not math.isfinite(value) or value < 0.0:
        raise HazardInputError("The exponential threshold must be finite and non-negative.")
    return value


def exponential_threshold(uniform: float) -> ExponentialThreshold:
    """``E = -log(U)`` for a caller-supplied ``U`` in the open interval (0, 1)."""
    if isinstance(uniform, bool) or not isinstance(uniform, (int, float, np.integer, np.floating)):
        raise HazardInputError("The uniform variate must be a number.")
    value = float(uniform)
    if not (0.0 < value < 1.0):
        raise HazardInputError(f"The uniform variate must lie in the open interval (0, 1); got {value!r}.")
    return -math.log(value)


# --- configuration bridge -------------------------------------------------------------


def schedule_from_config(config: BackgroundExposure) -> HazardSchedule:
    """Executable schedule for a validated ``background_exposure_v1`` configuration.

    Missing time-series values are left out, so any interval that touches them fails
    the coverage check; they are never read as zero. Age-specific hazards raise
    :class:`AgeSpecificHazardNotExecutableError` rather than being ignored.
    """
    if config.mode is ExposureMode.NONE:
        return HazardSchedule.none()
    if config.age_specific_hazards:
        raise AgeSpecificHazardNotExecutableError(
            "Age-specific background hazards are specified in background_exposure_v1 but are not yet "
            "executable. Ageing across band boundaries has not been approved; use one constant hazard."
        )
    if config.mode is ExposureMode.CONSTANT:
        return HazardSchedule.constant(config.constant_hazard.value)
    return HazardSchedule.annual_series(
        (row.year, row.hazard.value) for row in config.time_series if row.hazard.state is ValueState.VALUE
    )


def requires_random_draw(schedule: HazardSchedule) -> bool:
    """False for mode ``none``: callers must not draw a variate, so random streams are untouched."""
    return schedule.mode is not ExposureMode.NONE


def plausibility_warnings(schedule: HazardSchedule) -> tuple[str, ...]:
    """Non-blocking warnings for hazards above the provisional threshold. Not a validity check."""
    if schedule.mode is ExposureMode.CONSTANT:
        values = [("constant hazard", schedule.constant_hazard)]
    else:
        values = [(f"year {year}", hazard) for year, hazard in schedule.series]
    return tuple(hazard_plausibility_warning(h, c) for c, h in values if h > PLAUSIBILITY_WARNING_HAZARD)


# --- core mathematics -----------------------------------------------------------------


def hazard_at(schedule: HazardSchedule, t: float) -> float:
    """λ(t) in infections per susceptible person-year. Year ``y`` covers ``[y, y + 1)``."""
    t = _check_time(t, "t")
    if schedule.mode is ExposureMode.NONE:
        return 0.0
    if schedule.mode is ExposureMode.CONSTANT:
        return schedule.constant_hazard
    year = math.floor(t)
    for row_year, hazard in schedule.series:
        if row_year == year:
            return hazard
    raise HazardCoverageError(_coverage_message(schedule, [year], t, t))


def _coverage_message(schedule: HazardSchedule, missing: list[int], t0: float, t1: float) -> str:
    return (
        f"The background-exposure time series does not cover [{t0:g}, {t1:g}): no value for "
        f"year(s) {missing}. Supplied years: {list(schedule.years)}. Values are not carried "
        "forward, backfilled, interpolated or set to zero."
    )


def _segments(schedule: HazardSchedule, t0: float, t1: float) -> list[tuple[float, float, float]]:
    """Constant-hazard pieces ``(start, stop, λ)`` that partition ``[t0, t1)``; empty if ``t0 == t1``."""
    if t1 == t0:
        return []
    if schedule.mode is ExposureMode.NONE:
        return [(t0, t1, 0.0)]
    if schedule.mode is ExposureMode.CONSTANT:
        return [(t0, t1, schedule.constant_hazard)]
    first, last = math.floor(t0), math.ceil(t1) - 1
    lookup = dict(schedule.series)
    missing = [year for year in range(first, last + 1) if year not in lookup]
    if missing:
        raise HazardCoverageError(_coverage_message(schedule, missing, t0, t1))
    return [(max(t0, float(year)), min(t1, float(year + 1)), lookup[year]) for year in range(first, last + 1)]


def _total(segments: list[tuple[float, float, float]]) -> float:
    total = math.fsum(hazard * (stop - start) for start, stop, hazard in segments)
    if not math.isfinite(total):
        raise HazardInputError("Cumulative hazard overflowed; the hazard or interval is too large to represent.")
    return total


def cumulative_hazard(schedule: HazardSchedule, t0: float, t1: float) -> float:
    """``H(t0, t1) = ∫_{t0}^{t1} λ(s) ds`` over the half-open interval ``[t0, t1)``."""
    start, stop = _check_interval(t0, t1)
    if schedule.mode is ExposureMode.NONE:
        return 0.0
    return _total(_segments(schedule, start, stop))


def probability_from_cumulative_hazard(cumulative: float) -> float:
    """``1 - exp(-H)`` computed as ``-expm1(-H)`` (accurate for small H, exactly 1.0 for large H)."""
    value = _check_threshold(cumulative)
    return -math.expm1(-value)


def infection_probability(schedule: HazardSchedule, t0: float, t1: float) -> float:
    """Probability that a person susceptible at ``t0`` is first infected within ``[t0, t1)``."""
    return probability_from_cumulative_hazard(cumulative_hazard(schedule, t0, t1))


def first_infection_time(
    schedule: HazardSchedule,
    t0: float,
    t1: float,
    *,
    uniform: float | None = None,
    threshold: float | None = None,
) -> FirstInfectionResult:
    """Invert the cumulative hazard for one person susceptible at ``t0``.

    Supply exactly one of ``uniform`` (U in (0, 1); ``E = -log U``) or ``threshold``
    (E). Mode ``none`` takes neither and always returns no infection.
    """
    start, stop = _check_interval(t0, t1)
    if schedule.mode is ExposureMode.NONE:
        if uniform is not None or threshold is not None:
            raise HazardInputError("Mode 'none' uses no random variate; do not draw one.")
        return NoInfectionInInterval(cumulative_hazard=0.0, residual_threshold=None)
    if (uniform is None) == (threshold is None):
        raise HazardInputError("Supply exactly one of uniform or threshold.")
    target = exponential_threshold(uniform) if uniform is not None else _check_threshold(threshold)
    segments = _segments(schedule, start, stop)
    total = _total(segments)
    if not target < total:
        return NoInfectionInInterval(cumulative_hazard=total, residual_threshold=target - total)
    # Zero-hazard pieces add nothing to H, so the infection falls in a positive piece.
    positive = [segment for segment in segments if segment[2] > 0.0]
    accumulated = 0.0
    for index, (seg_start, seg_stop, hazard) in enumerate(positive):
        piece = hazard * (seg_stop - seg_start)
        if accumulated + piece > target or index == len(positive) - 1:
            # The last positive piece also absorbs any rounding difference from fsum.
            return InfectionOccurs(time=_place(seg_start + (target - accumulated) / hazard, seg_start, seg_stop))
        accumulated += piece
    raise AssertionError("unreachable: E < H implies a positive-hazard piece")


def _place(time: float, seg_start: float, seg_stop: float) -> float:
    """Keep an inverted time inside its half-open piece ``[seg_start, seg_stop)``."""
    return float(min(max(time, seg_start), np.nextafter(seg_stop, seg_start)))


# --- vectorised forms (same semantics as the scalar functions) ------------------------


def first_infection_times(
    schedule: HazardSchedule, t0: float, t1: float, thresholds: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised :func:`first_infection_time` for many people sharing ``[t0, t1)``.

    ``thresholds`` are exponential thresholds E (finite, non-negative). Returns
    ``(infected, times)``; ``times`` is NaN where no infection occurs. Mode ``none``
    rejects any thresholds array except an empty one: draw nothing.
    """
    start, stop = _check_interval(t0, t1)
    e = np.asarray(thresholds, dtype=float)
    if schedule.mode is ExposureMode.NONE:
        if e.size:
            raise HazardInputError("Mode 'none' uses no random variates; do not draw them.")
        return np.zeros(e.shape, dtype=bool), np.full(e.shape, np.nan)
    if not np.all(np.isfinite(e)) or np.any(e < 0.0):
        raise HazardInputError("Exponential thresholds must be finite and non-negative.")
    segments = _segments(schedule, start, stop)
    total = _total(segments)
    infected = e < total
    times = np.full(e.shape, np.nan)
    if not infected.any():
        return infected, times
    positive = [segment for segment in segments if segment[2] > 0.0]
    seg_start = np.array([s[0] for s in positive])
    seg_stop = np.array([s[1] for s in positive])
    hazard = np.array([s[2] for s in positive])
    pieces = hazard * (seg_stop - seg_start)
    cum_start = np.concatenate(([0.0], np.cumsum(pieces)[:-1]))
    # First positive piece whose cumulative end exceeds E; the last absorbs rounding.
    index = np.minimum(np.searchsorted(cum_start + pieces, e[infected], side="right"), len(positive) - 1)
    tau = seg_start[index] + (e[infected] - cum_start[index]) / hazard[index]
    times[infected] = np.clip(tau, seg_start[index], np.nextafter(seg_stop[index], seg_start[index]))
    return infected, times


def infection_probabilities(schedule: HazardSchedule, intervals: Iterable[tuple[float, float]]) -> np.ndarray:
    """:func:`infection_probability` for several intervals."""
    return np.array([infection_probability(schedule, t0, t1) for t0, t1 in intervals], dtype=float)
