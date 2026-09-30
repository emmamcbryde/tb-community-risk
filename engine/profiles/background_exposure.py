"""Versioned configuration for background exposure to TB infection (schema only).

Status: preparatory scaffolding. Nothing in the epidemiological engine, the engine
mapping or the ordinary interface reads this module; the current model has no
infection after baseline, which is mode ``none``. The scientific specification is
``docs/catalytic_infection_pressure_spec.md``.

Background exposure is an *exogenous force of infection* (a catalytic model input):
an annual infection hazard supplied from outside the model. It never depends on
infectious people generated within the model, and interventions never change it.

The schema keeps these distinctions:

* mode ``none`` means the hazard after baseline is exactly zero;
* a constant or time-series hazard of 0 is a value, not a missing value;
* a missing time-series value stays missing (its handling is an open decision);
* the quantity is an infection hazard in infections per person-year, never
  estimated TB disease incidence, and WHO incidence snapshots are not an accepted
  source.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import hashlib
import json
import math
from typing import Any

from engine.profiles.population_profile import (
    INCIDENCE_MEASURE,
    ProfileValidationError,
    ProfileValue,
    Provenance,
    ReviewStatus,
    ValueState,
)


BACKGROUND_EXPOSURE_SCHEMA_VERSION = "background_exposure_v1"
SUPPORTED_BACKGROUND_EXPOSURE_VERSIONS = (BACKGROUND_EXPOSURE_SCHEMA_VERSION,)
EXPOSURE_QUANTITY = "exogenous_infection_hazard"
HAZARD_UNIT = "infections per person-year"
# Provisional plausibility bound: a hazard of 1 per person-year is an annual infection
# probability of 1 - exp(-1), about 63%. Larger values are almost always unit errors
# (for example an incidence per 100,000 entered as a hazard).
MAX_PLAUSIBLE_HAZARD = 1.0
PROFILE_PAYLOAD_KEY = "backgroundExposure"
ACCEPTED_PROVENANCE = (Provenance.BUNDLED, Provenance.USER_DEFINED, Provenance.LOCAL_UPLOAD)


class ExposureMode(str, Enum):
    NONE = "none"
    CONSTANT = "constant"
    TIME_SERIES = "time_series"


class BackgroundExposureValidationError(ProfileValidationError):
    """Raised when a background-exposure payload is invalid."""


def _fail(message: str) -> None:
    raise BackgroundExposureValidationError(message)


def _hazard(payload: Any, context: str, *, allow_missing: bool) -> ProfileValue:
    if not isinstance(payload, dict):
        _fail(f"{context}: hazard must be an object with value, state, provenance, reviewStatus and unit.")
    try:
        value = ProfileValue.from_dict(payload)
    except ProfileValidationError as exc:
        raise BackgroundExposureValidationError(f"{context}: {exc}") from exc
    return _check_hazard(value, context, allow_missing=allow_missing)


def _check_hazard(value: ProfileValue, context: str, *, allow_missing: bool) -> ProfileValue:
    if not isinstance(value, ProfileValue):
        _fail(f"{context}: hazard must be a ProfileValue.")
    if value.unit != HAZARD_UNIT:
        _fail(f"{context}: unit must be '{HAZARD_UNIT}', not '{value.unit}'. TB disease incidence is not an infection hazard.")
    if value.state is ValueState.VALUE:
        if value.value is None or not math.isfinite(value.value):
            _fail(f"{context}: hazard must be a finite number.")
        if value.value < 0:
            _fail(f"{context}: hazard cannot be negative.")
        if value.value > MAX_PLAUSIBLE_HAZARD:
            _fail(f"{context}: hazard {value.value:g} per person-year exceeds {MAX_PLAUSIBLE_HAZARD:g}; check the units.")
    elif value.state is ValueState.MISSING:
        if not allow_missing:
            _fail(f"{context}: a value is required.")
        if value.value is not None:
            _fail(f"{context}: a missing hazard cannot carry a value.")
    else:
        _fail(f"{context}: state '{value.state.value}' is not accepted for a hazard.")
    if value.provenance not in ACCEPTED_PROVENANCE:
        _fail(f"{context}: provenance '{value.provenance.value}' is not accepted; WHO incidence snapshots do not supply infection hazards.")
    return value


@dataclass(frozen=True)
class AgeBandHazard:
    """Hazard for ages in [age_lower, age_upper); ``age_upper=None`` is open-ended."""

    age_lower: int
    age_upper: int | None
    hazard: ProfileValue

    def to_dict(self) -> dict[str, Any]:
        return {"ageLower": self.age_lower, "ageUpper": self.age_upper, "hazard": self.hazard.to_dict()}

    @classmethod
    def from_dict(cls, payload: dict[str, Any], index: int) -> "AgeBandHazard":
        context = f"ageSpecificHazards[{index}]"
        lower, upper = payload.get("ageLower"), payload.get("ageUpper")
        if isinstance(lower, bool) or not isinstance(lower, int) or lower < 0:
            _fail(f"{context}: ageLower must be a non-negative integer.")
        if upper is not None and (isinstance(upper, bool) or not isinstance(upper, int) or upper <= lower):
            _fail(f"{context}: ageUpper must be an integer greater than ageLower, or null.")
        return cls(lower, upper, _hazard(payload.get("hazard"), context, allow_missing=False))


@dataclass(frozen=True)
class HazardYear:
    """Annual hazard for one calendar year (applied over [year, year + 1))."""

    year: int
    hazard: ProfileValue

    def to_dict(self) -> dict[str, Any]:
        return {"year": self.year, "hazard": self.hazard.to_dict()}

    @classmethod
    def from_dict(cls, payload: dict[str, Any], index: int) -> "HazardYear":
        context = f"timeSeries[{index}]"
        year = payload.get("year")
        if isinstance(year, bool) or not isinstance(year, int) or not 1900 <= year <= 2200:
            _fail(f"{context}: year must be an integer between 1900 and 2200.")
        return cls(year, _hazard(payload.get("hazard"), f"{context} (year {year})", allow_missing=True))


@dataclass(frozen=True)
class BackgroundExposure:
    mode: ExposureMode = ExposureMode.NONE
    constant_hazard: ProfileValue | None = None
    age_specific_hazards: tuple[AgeBandHazard, ...] = ()
    time_series: tuple[HazardYear, ...] = ()
    source: str = ""
    citation: str = ""
    review_status: ReviewStatus = ReviewStatus.NOT_REQUIRED
    applicability_notes: str = ""
    schema_version: str = BACKGROUND_EXPOSURE_SCHEMA_VERSION
    quantity: str = EXPOSURE_QUANTITY
    unit: str = HAZARD_UNIT

    def __post_init__(self) -> None:
        _validate_structure(self)

    @property
    def hazard_is_identically_zero(self) -> bool:
        """True only for mode ``none``: no infection after baseline, as in the current model."""
        return self.mode is ExposureMode.NONE

    def to_dict(self) -> dict[str, Any]:
        return {
            "schemaVersion": self.schema_version,
            "mode": self.mode.value,
            "quantity": self.quantity,
            "unit": self.unit,
            "constantHazard": None if self.constant_hazard is None else self.constant_hazard.to_dict(),
            "ageSpecificHazards": [band.to_dict() for band in self.age_specific_hazards],
            "timeSeries": [row.to_dict() for row in self.time_series],
            "source": self.source,
            "citation": self.citation,
            "reviewStatus": self.review_status.value,
            "applicabilityNotes": self.applicability_notes,
        }

    def to_json(self, *, indent: int | None = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True, ensure_ascii=False, allow_nan=False)

    def exposure_hash(self) -> str:
        canonical = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    @classmethod
    def from_dict(cls, payload: dict[str, Any] | None) -> "BackgroundExposure":
        if payload is None:
            return no_background_exposure()
        if not isinstance(payload, dict):
            _fail("Background exposure must be an object.")
        version = payload.get("schemaVersion")
        if version not in SUPPORTED_BACKGROUND_EXPOSURE_VERSIONS:
            _fail(f"Unsupported background-exposure schema version: {version!r}.")
        quantity = payload.get("quantity", EXPOSURE_QUANTITY)
        if quantity == INCIDENCE_MEASURE:
            _fail("Estimated TB disease incidence cannot be used as background exposure to infection.")
        if quantity != EXPOSURE_QUANTITY:
            _fail(f"quantity must be '{EXPOSURE_QUANTITY}'.")
        unit = payload.get("unit", HAZARD_UNIT)
        if unit != HAZARD_UNIT:
            _fail(f"unit must be '{HAZARD_UNIT}', not '{unit}'. TB disease incidence is not an infection hazard.")
        try:
            mode = ExposureMode(payload.get("mode"))
        except ValueError:
            _fail(f"mode must be one of {[m.value for m in ExposureMode]}.")
        try:
            review_status = ReviewStatus(payload.get("reviewStatus", ReviewStatus.NOT_REQUIRED.value))
        except ValueError:
            _fail("reviewStatus is not recognised.")
        constant = payload.get("constantHazard")
        bands = payload.get("ageSpecificHazards") or []
        rows = payload.get("timeSeries") or []
        if not isinstance(bands, list) or not isinstance(rows, list):
            _fail("ageSpecificHazards and timeSeries must be lists.")
        return cls(
            mode=mode,
            constant_hazard=None if constant is None else _hazard(constant, "constantHazard", allow_missing=False),
            age_specific_hazards=tuple(sorted((AgeBandHazard.from_dict(b, i) for i, b in enumerate(bands)), key=lambda b: b.age_lower)),
            time_series=tuple(sorted((HazardYear.from_dict(r, i) for i, r in enumerate(rows)), key=lambda r: r.year)),
            source=str(payload.get("source") or ""),
            citation=str(payload.get("citation") or ""),
            review_status=review_status,
            applicability_notes=str(payload.get("applicabilityNotes") or ""),
        )

    @classmethod
    def from_json(cls, text: str) -> "BackgroundExposure":
        return cls.from_dict(json.loads(text))


def _validate_structure(config: BackgroundExposure) -> None:
    if config.schema_version not in SUPPORTED_BACKGROUND_EXPOSURE_VERSIONS:
        _fail(f"Unsupported background-exposure schema version: {config.schema_version!r}.")
    if config.quantity != EXPOSURE_QUANTITY or config.unit != HAZARD_UNIT:
        _fail("Background exposure must be an infection hazard in infections per person-year.")
    if not isinstance(config.mode, ExposureMode) or not isinstance(config.review_status, ReviewStatus):
        _fail("mode and review_status must be enum members.")
    if config.constant_hazard is not None:
        _check_hazard(config.constant_hazard, "constantHazard", allow_missing=False)
    for index, band in enumerate(config.age_specific_hazards):
        _check_hazard(band.hazard, f"ageSpecificHazards[{index}]", allow_missing=False)
    for index, row in enumerate(config.time_series):
        _check_hazard(row.hazard, f"timeSeries[{index}] (year {row.year})", allow_missing=True)
    has_constant = config.constant_hazard is not None
    has_bands = bool(config.age_specific_hazards)
    has_series = bool(config.time_series)
    if config.mode is ExposureMode.NONE:
        if has_constant or has_bands or has_series:
            _fail("Mode 'none' has an exactly zero hazard after baseline and takes no hazard values.")
        return
    if config.mode is ExposureMode.CONSTANT:
        if has_series:
            _fail("Mode 'constant' does not take a time series.")
        if has_constant == has_bands:
            _fail("Mode 'constant' needs either one constant hazard or a complete set of age-specific hazards.")
        if has_bands:
            _validate_age_bands(config.age_specific_hazards)
        return
    if has_constant or has_bands:
        _fail("Mode 'time_series' takes only time-series rows; age-by-year hazards are not yet specified.")
    if not has_series:
        _fail("Mode 'time_series' needs at least one row.")
    years = [row.year for row in config.time_series]
    if len(set(years)) != len(years):
        duplicates = sorted({year for year in years if years.count(year) > 1})
        _fail(f"Duplicate years in time series: {duplicates}.")
    if list(years) != sorted(years):
        _fail("Time-series rows must be ordered by year.")
    if all(row.hazard.state is ValueState.MISSING for row in config.time_series):
        _fail("A time series needs at least one non-missing hazard.")


def _validate_age_bands(bands: tuple[AgeBandHazard, ...]) -> None:
    if [band.age_lower for band in bands] != sorted(band.age_lower for band in bands):
        _fail("Age-specific hazards must be ordered by ageLower.")
    if bands[0].age_lower != 0:
        _fail("Age-specific hazards must start at age 0.")
    for previous, current in zip(bands, bands[1:]):
        if previous.age_upper is None or previous.age_upper != current.age_lower:
            _fail("Age-specific hazard bands must be contiguous and non-overlapping.")
    if bands[-1].age_upper is not None:
        _fail("The last age-specific hazard band must be open-ended (ageUpper null).")


def no_background_exposure() -> BackgroundExposure:
    """The current model: no infection after baseline."""
    return BackgroundExposure(mode=ExposureMode.NONE, review_status=ReviewStatus.NOT_REQUIRED)


def background_exposure_from_profile_payload(profile_payload: dict[str, Any]) -> BackgroundExposure:
    """Read an optional background-exposure block from a profile payload.

    Existing and legacy profiles have no such block and map to mode ``none``. The
    population-profile contract itself is unchanged.
    """
    return BackgroundExposure.from_dict(profile_payload.get(PROFILE_PAYLOAD_KEY))


def hazard_value(value: float | None, *, provenance: Provenance = Provenance.USER_DEFINED, source: str = "", notes: str = "") -> ProfileValue:
    """Build a hazard value; ``None`` is recorded as missing, never as zero."""
    return ProfileValue(
        value=None if value is None else float(value),
        state=ValueState.MISSING if value is None else ValueState.VALUE,
        provenance=provenance,
        review_status=ReviewStatus.NOT_REVIEWED,
        unit=HAZARD_UNIT,
        source=source or ("User-defined" if provenance is Provenance.USER_DEFINED else ""),
        notes=notes,
    )
