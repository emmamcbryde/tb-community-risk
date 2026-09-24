"""Serializable, versioned population-profile contract.

A population profile describes *who* is being modelled and *what evidence* the
inputs came from. It deliberately keeps these states distinct:

* zero            - a numeric value of 0 (``state=value``, ``value=0``);
* missing         - a value is expected but unavailable (``state=missing``);
* not reviewed    - a value exists but its evidence has not been reviewed;
* not applicable  - the quantity does not apply to this profile;
* excluded        - deliberately excluded from the analysis;
* bundled         - demonstration values supplied with the application;
* WHO snapshot    - taken from the installed, checksummed WHO data snapshot;
* local upload    - taken from a file supplied by the user (never labelled WHO);
* user-defined    - entered or changed by the user.

Profiles round-trip through JSON and hash deterministically so that saved
analyses can identify the exact profile and incidence snapshot used.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum
import hashlib
import json
import math
import re
from typing import Any


PROFILE_SCHEMA_VERSION = "population_profile_v2"
SUPPORTED_SCHEMA_VERSIONS = ("population_profile_v1", "population_profile_v2")
ISO3_PATTERN = re.compile(r"^[A-Z]{3}$")
INCIDENCE_MEASURE = "estimated_tb_disease_incidence"
INCIDENCE_UNIT = "per 100,000 population per year"


class ValueState(str, Enum):
    VALUE = "value"
    MISSING = "missing"
    NOT_APPLICABLE = "not_applicable"
    EXCLUDED = "excluded"


class Provenance(str, Enum):
    BUNDLED = "bundled"
    USER_DEFINED = "user_defined"
    WHO_SNAPSHOT = "who_snapshot"
    LOCAL_UPLOAD = "local_upload"


class ReviewStatus(str, Enum):
    NOT_REVIEWED = "not_reviewed"
    REVIEWED = "reviewed"
    NOT_REQUIRED = "not_required"


class EffectMeasure(str, Enum):
    """Effect-measure types are retained exactly; they are never converted."""

    RR = "RR"
    HR = "HR"
    OR = "OR"


class LocationKind(str, Enum):
    DEMONSTRATION = "demonstration"
    COUNTRY = "country"
    TERRITORY = "territory"
    SUBNATIONAL = "subnational"


STATUS_LABELS = {
    "zero": "Zero",
    "missing": "Missing",
    "not_reviewed": "Not reviewed",
    "not_applicable": "Not applicable",
    "excluded": "Excluded",
    "bundled": "Bundled",
    "user_defined": "User-defined",
    "who_snapshot": "WHO snapshot",
    "local_upload": "Local upload",
    "reviewed": "Reviewed",
}


class ProfileValidationError(ValueError):
    """Raised when a profile payload cannot be parsed into the contract."""


@dataclass(frozen=True)
class ProfileValue:
    """One numeric input with explicit value state, provenance and review status."""

    value: float | None
    state: ValueState
    provenance: Provenance
    review_status: ReviewStatus
    unit: str
    source: str = ""
    notes: str = ""

    @property
    def is_user_override(self) -> bool:
        return self.provenance is Provenance.USER_DEFINED

    @property
    def is_zero(self) -> bool:
        return self.state is ValueState.VALUE and self.value == 0

    def status_codes(self) -> tuple[str, ...]:
        """Return the distinct status codes describing this value."""
        codes: list[str] = []
        if self.state is ValueState.VALUE:
            if self.value == 0:
                codes.append("zero")
        else:
            codes.append(self.state.value)
        codes.append(self.provenance.value)
        if self.review_status is ReviewStatus.NOT_REVIEWED:
            codes.append("not_reviewed")
        elif self.review_status is ReviewStatus.REVIEWED:
            codes.append("reviewed")
        return tuple(codes)

    def status_label(self) -> str:
        return ", ".join(STATUS_LABELS[code] for code in self.status_codes())

    def to_dict(self) -> dict[str, Any]:
        return {
            "value": self.value,
            "state": self.state.value,
            "provenance": self.provenance.value,
            "reviewStatus": self.review_status.value,
            "unit": self.unit,
            "source": self.source,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "ProfileValue":
        _require_keys(payload, ("value", "state", "provenance", "reviewStatus", "unit"), "value")
        value = payload["value"]
        if value is not None:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ProfileValidationError(f"Non-numeric value: {value!r}")
            value = float(value)
        return cls(
            value=value,
            state=_enum(ValueState, payload["state"]),
            provenance=_enum(Provenance, payload["provenance"]),
            review_status=_enum(ReviewStatus, payload["reviewStatus"]),
            unit=str(payload["unit"]),
            source=str(payload.get("source") or ""),
            notes=str(payload.get("notes") or ""),
        )


def bundled_value(value: float | None, unit: str, source: str, *, notes: str = "") -> ProfileValue:
    state = ValueState.MISSING if value is None else ValueState.VALUE
    return ProfileValue(
        value=None if value is None else float(value),
        state=state,
        provenance=Provenance.BUNDLED,
        review_status=ReviewStatus.NOT_REVIEWED,
        unit=unit,
        source=source,
        notes=notes,
    )


def user_value(value: float | None, unit: str, *, notes: str = "") -> ProfileValue:
    """Return a user-defined value; ``None`` is recorded as missing, never as zero."""
    state = ValueState.MISSING if value is None else ValueState.VALUE
    return ProfileValue(
        value=None if value is None else float(value),
        state=state,
        provenance=Provenance.USER_DEFINED,
        review_status=ReviewStatus.NOT_REVIEWED,
        unit=unit,
        source="User-defined",
        notes=notes,
    )


def not_applicable_value(unit: str, *, notes: str = "") -> ProfileValue:
    return ProfileValue(
        value=None,
        state=ValueState.NOT_APPLICABLE,
        provenance=Provenance.BUNDLED,
        review_status=ReviewStatus.NOT_REQUIRED,
        unit=unit,
        notes=notes,
    )


@dataclass(frozen=True)
class Location:
    """Identity of the modelled place; the national population is context only.

    ``national_population`` is the WHO/UN population of the whole country or area;
    it is not the simulated population size.
    """

    name: str
    kind: LocationKind
    iso3: str | None = None
    who_region: str | None = None
    national_population: float | None = None
    national_population_year: int | None = None
    national_population_source: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "kind": self.kind.value,
            "iso3": self.iso3,
            "whoRegion": self.who_region,
            "nationalPopulation": self.national_population,
            "nationalPopulationYear": self.national_population_year,
            "nationalPopulationSource": self.national_population_source,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "Location":
        _require_keys(payload, ("name", "kind"), "location")
        iso3 = payload.get("iso3")
        year = payload.get("nationalPopulationYear")
        return cls(
            name=str(payload["name"]),
            kind=_enum(LocationKind, payload["kind"]),
            iso3=None if iso3 in (None, "") else str(iso3),
            who_region=payload.get("whoRegion") or None,
            national_population=_optional_float(payload.get("nationalPopulation")),
            national_population_year=None if year is None else int(year),
            national_population_source=str(payload.get("nationalPopulationSource") or ""),
        )


@dataclass(frozen=True)
class AgeBand:
    label: str
    lower_age: int
    upper_age: int | None
    proportion: ProfileValue

    def to_dict(self) -> dict[str, Any]:
        return {
            "label": self.label,
            "lowerAge": self.lower_age,
            "upperAge": self.upper_age,
            "proportion": self.proportion.to_dict(),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "AgeBand":
        _require_keys(payload, ("label", "lowerAge", "upperAge", "proportion"), "age band")
        upper = payload["upperAge"]
        return cls(
            label=str(payload["label"]),
            lower_age=int(payload["lowerAge"]),
            upper_age=None if upper is None else int(upper),
            proportion=ProfileValue.from_dict(payload["proportion"]),
        )


@dataclass(frozen=True)
class IncidencePoint:
    """Estimated TB disease incidence for one year (per 100,000 per year)."""

    year: int
    estimate: float | None
    lower: float | None
    upper: float | None

    def to_dict(self) -> dict[str, Any]:
        return {"year": self.year, "estimate": self.estimate, "lower": self.lower, "upper": self.upper}

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "IncidencePoint":
        _require_keys(payload, ("year", "estimate", "lower", "upper"), "incidence point")
        return cls(
            year=int(payload["year"]),
            estimate=_optional_float(payload["estimate"]),
            lower=_optional_float(payload["lower"]),
            upper=_optional_float(payload["upper"]),
        )


@dataclass(frozen=True)
class IncidenceData:
    """Observed incidence series attached to a profile.

    This is *estimated TB disease incidence*. It is not infection pressure, force
    of infection, progression from infection to disease, or case notification.
    """

    source: str
    snapshot_id: str | None
    provenance: Provenance
    series: tuple[IncidencePoint, ...] = ()
    measure: str = INCIDENCE_MEASURE
    unit: str = INCIDENCE_UNIT
    notes: str = ""
    source_detail: tuple[tuple[str, Any], ...] = ()

    @property
    def data_year_range(self) -> tuple[int, int] | None:
        years = [point.year for point in self.series]
        return (min(years), max(years)) if years else None

    @property
    def bounds_available(self) -> bool:
        return bool(self.series) and all(p.lower is not None and p.upper is not None for p in self.series)

    @property
    def data_hash(self) -> str | None:
        """SHA-256 of the canonical series (year, estimate, lower, upper)."""
        if not self.series:
            return None
        canonical = json.dumps([point.to_dict() for point in self.series], sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        year_range = self.data_year_range
        return {
            "source": self.source,
            "snapshotId": self.snapshot_id,
            "provenance": self.provenance.value,
            "measure": self.measure,
            "unit": self.unit,
            "dataYearRange": list(year_range) if year_range else None,
            "series": [point.to_dict() for point in self.series],
            "notes": self.notes,
            "sourceDetail": {key: value for key, value in self.source_detail},
            "dataHash": self.data_hash,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "IncidenceData":
        _require_keys(payload, ("source", "snapshotId", "provenance", "series"), "incidence")
        series = tuple(IncidencePoint.from_dict(point) for point in payload["series"] or [])
        data = cls(
            source=str(payload["source"]),
            snapshot_id=payload["snapshotId"],
            provenance=_enum(Provenance, payload["provenance"]),
            series=series,
            measure=str(payload.get("measure") or INCIDENCE_MEASURE),
            unit=str(payload.get("unit") or INCIDENCE_UNIT),
            notes=str(payload.get("notes") or ""),
            source_detail=tuple(sorted((str(k), _freeze(v)) for k, v in (payload.get("sourceDetail") or {}).items())),
        )
        declared_hash = payload.get("dataHash")
        if declared_hash is not None and declared_hash != data.data_hash:
            raise ProfileValidationError("Incidence dataHash does not match the series.")
        declared = payload.get("dataYearRange")
        derived = data.data_year_range
        if declared is not None and (derived is None or tuple(declared) != derived):
            raise ProfileValidationError("dataYearRange does not match the incidence series.")
        return data


@dataclass(frozen=True)
class TrendSpec:
    """Incidence-trend method and settings; see ``engine.who_incidence.trend``."""

    method: str = "not_estimated"
    settings: tuple[tuple[str, Any], ...] = ()
    user_annual_percent_change: ProfileValue | None = None

    def settings_dict(self) -> dict[str, Any]:
        return dict(self.settings)

    def to_dict(self) -> dict[str, Any]:
        return {
            "method": self.method,
            "settings": {key: value for key, value in self.settings},
            "userAnnualPercentChange": (
                None if self.user_annual_percent_change is None else self.user_annual_percent_change.to_dict()
            ),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "TrendSpec":
        _require_keys(payload, ("method",), "trend")
        override = payload.get("userAnnualPercentChange")
        settings = payload.get("settings") or {}
        return cls(
            method=str(payload["method"]),
            settings=tuple(sorted((str(key), _freeze(value)) for key, value in settings.items())),
            user_annual_percent_change=None if override is None else ProfileValue.from_dict(override),
        )


@dataclass(frozen=True)
class RiskFactor:
    """Optional risk-factor stratum.

    ``effect_measure`` is retained exactly as entered (RR, HR or OR). The current
    engine applies the effect estimate as a multiplier on the progression hazard
    (``engine_application``); that is an internal engine assumption and no
    conversion between measure types is performed.
    """

    risk_factor_id: str
    label: str
    enabled: bool
    prevalence: ProfileValue
    effect_estimate: ProfileValue
    effect_measure: EffectMeasure | None
    engine_key: str | None = None
    engine_application: str | None = None
    evidence_source: str = ""
    review_status: ReviewStatus = ReviewStatus.NOT_REVIEWED
    notes: str = ""
    user_modified_fields: tuple[str, ...] = ()
    prevalence_bounds: tuple[float, float] | None = None
    effect_bounds: tuple[float, float] | None = None
    affected_transition: str = "progression_to_disease"
    evidence_year: int | None = None
    population_applicability: str = ""

    @property
    def is_user_override(self) -> bool:
        return bool(self.user_modified_fields) or self.prevalence.is_user_override or self.effect_estimate.is_user_override

    def to_dict(self) -> dict[str, Any]:
        return {
            "riskFactorId": self.risk_factor_id,
            "label": self.label,
            "enabled": self.enabled,
            "prevalence": self.prevalence.to_dict(),
            "effectEstimate": self.effect_estimate.to_dict(),
            "effectMeasure": None if self.effect_measure is None else self.effect_measure.value,
            "engineKey": self.engine_key,
            "engineApplication": self.engine_application,
            "evidenceSource": self.evidence_source,
            "reviewStatus": self.review_status.value,
            "notes": self.notes,
            "userModifiedFields": list(self.user_modified_fields),
            "userOverride": self.is_user_override,
            "prevalenceBounds": None if self.prevalence_bounds is None else list(self.prevalence_bounds),
            "effectBounds": None if self.effect_bounds is None else list(self.effect_bounds),
            "affectedTransition": self.affected_transition,
            "evidenceYear": self.evidence_year,
            "populationApplicability": self.population_applicability,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "RiskFactor":
        _require_keys(
            payload,
            ("riskFactorId", "label", "enabled", "prevalence", "effectEstimate", "effectMeasure"),
            "risk factor",
        )
        measure = payload["effectMeasure"]
        factor = cls(
            risk_factor_id=str(payload["riskFactorId"]),
            label=str(payload["label"]),
            enabled=_strict_bool(payload["enabled"], "enabled"),
            prevalence=ProfileValue.from_dict(payload["prevalence"]),
            effect_estimate=ProfileValue.from_dict(payload["effectEstimate"]),
            effect_measure=None if measure is None else _enum(EffectMeasure, measure),
            engine_key=payload.get("engineKey"),
            engine_application=payload.get("engineApplication"),
            evidence_source=str(payload.get("evidenceSource") or ""),
            review_status=_enum(ReviewStatus, payload.get("reviewStatus") or ReviewStatus.NOT_REVIEWED.value),
            notes=str(payload.get("notes") or ""),
            user_modified_fields=tuple(str(item) for item in payload.get("userModifiedFields") or ()),
            prevalence_bounds=_optional_pair(payload.get("prevalenceBounds")),
            effect_bounds=_optional_pair(payload.get("effectBounds")),
            affected_transition=str(payload.get("affectedTransition") or "progression_to_disease"),
            evidence_year=None if payload.get("evidenceYear") in (None, "") else int(payload["evidenceYear"]),
            population_applicability=str(payload.get("populationApplicability") or ""),
        )
        declared = payload.get("userOverride")
        if declared is not None and bool(declared) != factor.is_user_override:
            raise ProfileValidationError(f"userOverride is inconsistent for risk factor {factor.risk_factor_id}.")
        return factor


@dataclass(frozen=True)
class PopulationProfile:
    profile_id: str
    name: str
    location: Location
    population_size: ProfileValue
    age_distribution: tuple[AgeBand, ...]
    ltbi_prevalence: ProfileValue
    incidence: IncidenceData
    trend: TrendSpec = field(default_factory=TrendSpec)
    risk_factors: tuple[RiskFactor, ...] = ()
    data_vintage: str = ""
    notes: str = ""
    demonstration: bool = False
    age_distribution_source: str = ""
    schema_version: str = PROFILE_SCHEMA_VERSION

    @property
    def has_user_overrides(self) -> bool:
        return bool(user_override_fields(self))

    def to_dict(self) -> dict[str, Any]:
        return {
            "schemaVersion": self.schema_version,
            "profileId": self.profile_id,
            "name": self.name,
            "location": self.location.to_dict(),
            "populationSize": self.population_size.to_dict(),
            "ageDistribution": [band.to_dict() for band in self.age_distribution],
            "ageDistributionSource": self.age_distribution_source,
            "ltbiPrevalence": self.ltbi_prevalence.to_dict(),
            "incidence": self.incidence.to_dict(),
            "trend": self.trend.to_dict(),
            "riskFactors": [factor.to_dict() for factor in self.risk_factors],
            "dataVintage": self.data_vintage,
            "notes": self.notes,
            "demonstration": self.demonstration,
        }

    def to_json(self, *, indent: int | None = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True, ensure_ascii=False, allow_nan=False)

    def profile_hash(self) -> str:
        return profile_hash(self)

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "PopulationProfile":
        if not isinstance(payload, dict):
            raise ProfileValidationError("Profile payload must be a JSON object.")
        version = payload.get("schemaVersion")
        if version not in SUPPORTED_SCHEMA_VERSIONS:
            raise ProfileValidationError(
                f"Unsupported profile schema version {version!r}; expected one of {', '.join(SUPPORTED_SCHEMA_VERSIONS)}."
            )
        if version == "population_profile_v1":
            payload = _migrate_v1(payload)
        _require_keys(
            payload,
            ("profileId", "name", "location", "populationSize", "ageDistribution", "ltbiPrevalence", "incidence"),
            "profile",
        )
        return cls(
            schema_version=PROFILE_SCHEMA_VERSION,
            profile_id=str(payload["profileId"]),
            name=str(payload["name"]),
            location=Location.from_dict(payload["location"]),
            population_size=ProfileValue.from_dict(payload["populationSize"]),
            age_distribution=tuple(AgeBand.from_dict(band) for band in payload["ageDistribution"] or []),
            age_distribution_source=str(payload.get("ageDistributionSource") or ""),
            ltbi_prevalence=ProfileValue.from_dict(payload["ltbiPrevalence"]),
            incidence=IncidenceData.from_dict(payload["incidence"]),
            trend=TrendSpec.from_dict(payload.get("trend") or {"method": "not_estimated"}),
            risk_factors=tuple(RiskFactor.from_dict(item) for item in payload.get("riskFactors") or []),
            data_vintage=str(payload.get("dataVintage") or ""),
            notes=str(payload.get("notes") or ""),
            demonstration=_strict_bool(payload.get("demonstration", False), "demonstration"),
        )

    @classmethod
    def from_json(cls, text: str) -> "PopulationProfile":
        try:
            payload = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ProfileValidationError(f"Profile is not valid JSON: {exc}") from exc
        return cls.from_dict(payload)


def profile_hash(profile: PopulationProfile) -> str:
    """Deterministic SHA-256 of the canonical JSON form of a profile."""
    canonical = json.dumps(
        profile.to_dict(),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def user_override_fields(profile: PopulationProfile) -> list[str]:
    """Return human-readable names of all user-defined inputs in a profile."""
    fields: list[str] = []
    if profile.population_size.is_user_override:
        fields.append("Population size")
    if profile.ltbi_prevalence.is_user_override:
        fields.append("LTBI prevalence")
    for band in profile.age_distribution:
        if band.proportion.is_user_override:
            fields.append(f"Age distribution: {band.label}")
    if profile.incidence.provenance in {Provenance.USER_DEFINED, Provenance.LOCAL_UPLOAD}:
        fields.append("Incidence series (local file)")
    if profile.trend.user_annual_percent_change is not None and profile.trend.user_annual_percent_change.is_user_override:
        fields.append("Incidence trend")
    for factor in profile.risk_factors:
        if factor.is_user_override:
            fields.append(f"Risk factor: {factor.label}")
    return fields


def validate_profile(profile: PopulationProfile) -> list[dict[str, str]]:
    """Return validation issues; an empty list means the profile is valid."""
    issues: list[dict[str, str]] = []

    def add(field_name: str, message: str, severity: str = "error") -> None:
        issues.append({"field": field_name, "severity": severity, "message": message})

    for name, value in _all_values(profile):
        for message in _value_consistency_messages(value):
            add(name, message)

    if profile.location.iso3 is not None and not ISO3_PATTERN.match(profile.location.iso3):
        add("location.iso3", f"Invalid ISO3 code {profile.location.iso3!r}.")

    size = profile.population_size
    if size.state is not ValueState.VALUE:
        add("populationSize", "Population size is required.")
    elif size.value is not None and (size.value < 1 or size.value != int(size.value)):
        add("populationSize", "Population size must be a positive whole number.")

    _check_proportion(profile.ltbi_prevalence, "ltbiPrevalence", add)

    proportions = []
    for band in profile.age_distribution:
        _check_proportion(band.proportion, f"ageDistribution.{band.label}", add)
        if band.proportion.state is ValueState.VALUE and band.proportion.value is not None:
            proportions.append(band.proportion.value)
    if proportions and len(proportions) == len(profile.age_distribution):
        if abs(sum(proportions) - 1.0) > 1e-6:
            add("ageDistribution", f"Age proportions sum to {sum(proportions):.6f}; expected 1.")

    for point in profile.incidence.series:
        prefix = f"incidence.{point.year}"
        for label, value in (("estimate", point.estimate), ("lower", point.lower), ("upper", point.upper)):
            if value is not None and value < 0:
                add(prefix, f"Negative incidence {label}.")
        if point.estimate is not None and point.lower is not None and point.lower > point.estimate:
            add(prefix, "Lower bound exceeds the estimate.")
        if point.estimate is not None and point.upper is not None and point.upper < point.estimate:
            add(prefix, "Upper bound is below the estimate.")
    years = [point.year for point in profile.incidence.series]
    if len(years) != len(set(years)):
        add("incidence", "Duplicate incidence years.")

    seen_ids: set[str] = set()
    for factor in profile.risk_factors:
        prefix = f"riskFactors.{factor.risk_factor_id}"
        if factor.risk_factor_id in seen_ids:
            add(prefix, "Duplicate risk-factor identifier.")
        seen_ids.add(factor.risk_factor_id)
        _check_proportion(factor.prevalence, f"{prefix}.prevalence", add)
        for label, bounds, value in (
            ("prevalenceBounds", factor.prevalence_bounds, factor.prevalence.value),
            ("effectBounds", factor.effect_bounds, factor.effect_estimate.value),
        ):
            if bounds is None:
                continue
            low, high = bounds
            if low > high:
                add(f"{prefix}.{label}", "Lower bound exceeds upper bound.")
            elif value is not None and not low <= value <= high:
                add(f"{prefix}.{label}", "Estimate lies outside its bounds.")
            if label == "prevalenceBounds" and not (0 <= low <= 1 and 0 <= high <= 1):
                add(f"{prefix}.{label}", "Prevalence bounds must be between 0 and 1.")
            if label == "effectBounds" and low <= 0:
                add(f"{prefix}.{label}", "Effect bounds must be greater than zero.")
        if factor.evidence_year is not None and not 1900 <= factor.evidence_year <= 2100:
            add(f"{prefix}.evidenceYear", "Evidence year is implausible.")
        effect = factor.effect_estimate
        if effect.state is ValueState.VALUE:
            if effect.value is not None and effect.value <= 0:
                add(f"{prefix}.effectEstimate", "Effect estimate must be greater than zero.")
            if factor.effect_measure is None:
                add(f"{prefix}.effectMeasure", "An effect estimate requires an effect-measure type (RR, HR or OR).")
        if factor.enabled:
            if factor.prevalence.state is ValueState.MISSING:
                add(f"{prefix}.prevalence", "Enabled risk factor has no prevalence.", "blocking")
            if effect.state is ValueState.MISSING:
                add(f"{prefix}.effectEstimate", "Enabled risk factor has no effect estimate.", "blocking")
    return issues


def unresolved_inputs(profile: PopulationProfile) -> list[dict[str, str]]:
    """List inputs that still need locally applicable evidence."""
    rows: list[dict[str, str]] = []

    def add(item: str, reason: str) -> None:
        rows.append({"Input": item, "Why it needs attention": reason})

    def check(item: str, value: ProfileValue) -> None:
        if value.state is ValueState.MISSING:
            add(item, "Missing")
        elif value.state is ValueState.VALUE and value.review_status is ReviewStatus.NOT_REVIEWED:
            if value.provenance is Provenance.BUNDLED:
                add(item, "Demonstration value; not reviewed for this setting")
            else:
                add(item, "User-defined; evidence not yet reviewed")

    check("LTBI prevalence", profile.ltbi_prevalence)
    if profile.age_distribution and any(
        band.proportion.review_status is ReviewStatus.NOT_REVIEWED for band in profile.age_distribution
    ):
        if any(band.proportion.is_user_override for band in profile.age_distribution):
            add("Age distribution", "User-defined; evidence not yet reviewed")
        else:
            add("Age distribution", "Demonstration value; not reviewed for this setting")
    if not profile.incidence.series:
        add("TB incidence", "No incidence series linked to this profile")
    for factor in profile.risk_factors:
        if not factor.enabled:
            continue
        check(f"{factor.label}: prevalence", factor.prevalence)
        check(f"{factor.label}: effect estimate", factor.effect_estimate)
    return rows


def with_population_size(profile: PopulationProfile, size: int) -> PopulationProfile:
    if profile.population_size.state is ValueState.VALUE and profile.population_size.value == float(size):
        return profile
    return replace(profile, population_size=user_value(float(size), profile.population_size.unit))


def _all_values(profile: PopulationProfile) -> list[tuple[str, ProfileValue]]:
    values = [("populationSize", profile.population_size), ("ltbiPrevalence", profile.ltbi_prevalence)]
    values.extend((f"ageDistribution.{band.label}", band.proportion) for band in profile.age_distribution)
    for factor in profile.risk_factors:
        values.append((f"riskFactors.{factor.risk_factor_id}.prevalence", factor.prevalence))
        values.append((f"riskFactors.{factor.risk_factor_id}.effectEstimate", factor.effect_estimate))
    if profile.trend.user_annual_percent_change is not None:
        values.append(("trend.userAnnualPercentChange", profile.trend.user_annual_percent_change))
    return values


def _value_consistency_messages(value: ProfileValue) -> list[str]:
    if value.state is ValueState.VALUE:
        if value.value is None or not math.isfinite(value.value):
            return ["A value marked as present must be a finite number."]
        return []
    if value.value is not None:
        return [f"A value marked {value.state.value} must not carry a number."]
    return []


def _check_proportion(value: ProfileValue, name: str, add) -> None:
    if value.state is ValueState.VALUE and value.value is not None and not 0 <= value.value <= 1:
        add(name, "Proportion must be between 0 and 1.")


def _require_keys(payload: Any, keys: tuple[str, ...], label: str) -> None:
    if not isinstance(payload, dict):
        raise ProfileValidationError(f"{label} must be an object.")
    missing = [key for key in keys if key not in payload]
    if missing:
        raise ProfileValidationError(f"{label} is missing required field(s): {', '.join(missing)}.")


def _enum(enum_type, value: Any):
    try:
        return enum_type(value)
    except ValueError as exc:
        allowed = ", ".join(member.value for member in enum_type)
        raise ProfileValidationError(f"Invalid {enum_type.__name__} {value!r}; expected one of: {allowed}.") from exc


def _strict_bool(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise ProfileValidationError(f"{label} must be true or false.")
    return value


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ProfileValidationError(f"Non-numeric value: {value!r}")
    return float(value)


def _optional_pair(value: Any) -> tuple[float, float] | None:
    if value in (None, [], ()):
        return None
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ProfileValidationError("Bounds must be a pair [lower, upper].")
    low, high = (_optional_float(item) for item in value)
    if low is None or high is None:
        raise ProfileValidationError("Bounds must contain two numbers.")
    return (low, high)


def _migrate_v1(payload: dict[str, Any]) -> dict[str, Any]:
    """Upgrade a v1 profile: WHO-snapshot and uploaded incidence get explicit provenance."""
    migrated = json.loads(json.dumps(payload))
    incidence = migrated.get("incidence") or {}
    if incidence.get("provenance") == "bundled" and incidence.get("snapshotId"):
        incidence["provenance"] = Provenance.WHO_SNAPSHOT.value
    elif incidence.get("provenance") == "user_defined":
        incidence["provenance"] = Provenance.LOCAL_UPLOAD.value
    incidence.pop("dataHash", None)
    migrated["schemaVersion"] = PROFILE_SCHEMA_VERSION
    return migrated


def _freeze(value: Any) -> Any:
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, dict):
        return tuple(sorted((str(k), _freeze(v)) for k, v in value.items()))
    return value
