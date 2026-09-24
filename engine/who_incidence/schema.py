"""Application schema and validation for country-level TB incidence snapshots."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import re
from typing import Any, Iterable


SNAPSHOT_SCHEMA_VERSION = "who_incidence_snapshot_v1"
TRANSFORMATION_VERSION = "who_incidence_transform_v1"
MEASURE = "estimated_tb_disease_incidence"
RATE_UNIT = "per 100,000 population per year"

RECORD_COLUMNS: tuple[str, ...] = (
    "iso3",
    "country",
    "year",
    "incidence_per_100k",
    "incidence_per_100k_lo",
    "incidence_per_100k_hi",
    "incident_cases",
    "incident_cases_lo",
    "incident_cases_hi",
    "population",
    "estimation_method",
)
REQUIRED_COLUMNS: tuple[str, ...] = (
    "iso3",
    "year",
    "incidence_per_100k",
    "incidence_per_100k_lo",
    "incidence_per_100k_hi",
)
NUMERIC_COLUMNS: tuple[str, ...] = (
    "incidence_per_100k",
    "incidence_per_100k_lo",
    "incidence_per_100k_hi",
    "incident_cases",
    "incident_cases_lo",
    "incident_cases_hi",
    "population",
)
BOUND_TRIPLES: tuple[tuple[str, str, str], ...] = (
    ("incidence_per_100k", "incidence_per_100k_lo", "incidence_per_100k_hi"),
    ("incident_cases", "incident_cases_lo", "incident_cases_hi"),
)
MANIFEST_REQUIRED_FIELDS: tuple[str, ...] = (
    "schemaVersion",
    "snapshotId",
    "snapshotKind",
    "isCompleteDataset",
    "sourceDataset",
    "sourceReportYear",
    "upstream",
    "extractionDate",
    "transformationVersion",
    "dataFile",
    "columns",
    "licence",
    "citation",
    "validationStatus",
)
ISO3_PATTERN = re.compile(r"^[A-Z]{3}$")
MIN_YEAR = 1900
MAX_YEAR = 2100
# Year-on-year ratio beyond which a change is flagged for review (not rejected).
DISCONTINUITY_RATIO = 2.0
DISCONTINUITY_MIN_RATE = 1.0


@dataclass(frozen=True)
class ValidationIssue:
    code: str
    severity: str
    message: str
    iso3: str | None = None
    year: int | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "severity": self.severity,
            "message": self.message,
            "iso3": self.iso3,
            "year": self.year,
        }


@dataclass
class ValidationReport:
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def errors(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.severity == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.severity == "warning"]

    @property
    def is_valid(self) -> bool:
        return not self.errors

    def codes(self) -> set[str]:
        return {issue.code for issue in self.issues}

    def add(self, code: str, severity: str, message: str, *, iso3: str | None = None, year: int | None = None) -> None:
        self.issues.append(ValidationIssue(code, severity, message, iso3, year))

    def extend(self, other: "ValidationReport") -> None:
        self.issues.extend(other.issues)

    def to_dict(self) -> dict[str, Any]:
        return {
            "isValid": self.is_valid,
            "errorCount": len(self.errors),
            "warningCount": len(self.warnings),
            "issues": [issue.to_dict() for issue in self.issues],
        }


@dataclass(frozen=True)
class IncidenceRecord:
    iso3: str
    country: str
    year: int
    incidence_per_100k: float
    incidence_per_100k_lo: float
    incidence_per_100k_hi: float
    incident_cases: float | None = None
    incident_cases_lo: float | None = None
    incident_cases_hi: float | None = None
    population: float | None = None
    estimation_method: str = ""

    def to_row(self) -> dict[str, Any]:
        return {column: getattr(self, column) for column in RECORD_COLUMNS}


def validate_columns(
    columns: Iterable[str],
    *,
    expected: tuple[str, ...] = RECORD_COLUMNS,
    required: tuple[str, ...] = REQUIRED_COLUMNS,
    warn_missing_optional: bool = True,
) -> ValidationReport:
    """Detect schema changes: missing required columns are errors, extras warnings."""
    report = ValidationReport()
    present = list(columns)
    missing_required = [column for column in required if column not in present]
    if missing_required:
        report.add("schema_changed", "error", f"Missing required column(s): {', '.join(missing_required)}.")
    missing_optional = [column for column in expected if column not in present and column not in required]
    if missing_optional and warn_missing_optional:
        report.add("schema_changed", "warning", f"Missing optional column(s): {', '.join(missing_optional)}.")
    unexpected = [column for column in present if column not in expected]
    if unexpected:
        report.add("schema_changed", "warning", f"Unexpected column(s): {', '.join(unexpected)}.")
    return report


def parse_rows(
    rows: Iterable[dict[str, Any]],
    *,
    require_iso3: bool = True,
) -> tuple[list[IncidenceRecord], ValidationReport]:
    """Parse raw rows (e.g. from ``csv.DictReader``) into validated records."""
    rows = list(rows)
    report = ValidationReport()
    columns: list[str] = list(rows[0].keys()) if rows else list(RECORD_COLUMNS)
    required = REQUIRED_COLUMNS if require_iso3 else tuple(column for column in REQUIRED_COLUMNS if column != "iso3")
    report.extend(validate_columns(columns, required=required, warn_missing_optional=require_iso3))
    if not report.is_valid:
        return [], report
    if not rows:
        report.add("empty", "error", "The incidence file contains no rows.")
        return [], report

    records: list[IncidenceRecord] = []
    for index, raw in enumerate(rows, start=2):
        iso3 = str(raw.get("iso3") or "").strip().upper()
        year_value = _parse_year(raw.get("year"))
        if require_iso3 or iso3:
            if not ISO3_PATTERN.match(iso3):
                report.add("invalid_iso3", "error", f"Row {index}: invalid ISO3 code {raw.get('iso3')!r}.", year=year_value)
                continue
        if year_value is None:
            report.add("invalid_year", "error", f"Row {index}: invalid year {raw.get('year')!r}.", iso3=iso3 or None)
            continue
        numbers: dict[str, float | None] = {}
        row_ok = True
        for column in NUMERIC_COLUMNS:
            parsed, ok = _parse_number(raw.get(column))
            if not ok:
                report.add(
                    "non_numeric",
                    "error",
                    f"Row {index}: non-numeric {column} {raw.get(column)!r}.",
                    iso3=iso3 or None,
                    year=year_value,
                )
                row_ok = False
            numbers[column] = parsed
        if not row_ok:
            continue
        for column in REQUIRED_COLUMNS[2:]:
            if numbers[column] is None:
                report.add("missing_value", "error", f"Row {index}: {column} is missing.", iso3=iso3 or None, year=year_value)
                row_ok = False
        if not row_ok:
            continue
        records.append(
            IncidenceRecord(
                iso3=iso3,
                country=str(raw.get("country") or "").strip(),
                year=year_value,
                incidence_per_100k=numbers["incidence_per_100k"],
                incidence_per_100k_lo=numbers["incidence_per_100k_lo"],
                incidence_per_100k_hi=numbers["incidence_per_100k_hi"],
                incident_cases=numbers["incident_cases"],
                incident_cases_lo=numbers["incident_cases_lo"],
                incident_cases_hi=numbers["incident_cases_hi"],
                population=numbers["population"],
                estimation_method=str(raw.get("estimation_method") or "").strip(),
            )
        )
    report.extend(validate_records(records))
    return records, report


def validate_records(records: list[IncidenceRecord]) -> ValidationReport:
    """Validate parsed records for duplicates, bounds, gaps and discontinuities."""
    report = ValidationReport()
    seen: set[tuple[str, int]] = set()
    by_country: dict[str, list[IncidenceRecord]] = {}
    for record in records:
        key = (record.iso3, record.year)
        if key in seen:
            report.add("duplicate_row", "error", f"Duplicate row for {record.iso3} {record.year}.", iso3=record.iso3, year=record.year)
        seen.add(key)
        by_country.setdefault(record.iso3, []).append(record)
        if not MIN_YEAR <= record.year <= MAX_YEAR:
            report.add("invalid_year", "error", f"Year {record.year} is out of range.", iso3=record.iso3, year=record.year)
        for column in NUMERIC_COLUMNS:
            value = getattr(record, column)
            if value is not None and value < 0:
                report.add("negative_value", "error", f"Negative {column}.", iso3=record.iso3, year=record.year)
        for estimate_col, lo_col, hi_col in BOUND_TRIPLES:
            estimate, lower, upper = (getattr(record, col) for col in (estimate_col, lo_col, hi_col))
            if estimate is None:
                continue
            if lower is not None and lower > estimate:
                report.add("lower_exceeds_estimate", "error", f"{lo_col} exceeds {estimate_col}.", iso3=record.iso3, year=record.year)
            if upper is not None and upper < estimate:
                report.add("upper_below_estimate", "error", f"{hi_col} is below {estimate_col}.", iso3=record.iso3, year=record.year)

    for iso3, country_records in by_country.items():
        years = sorted({record.year for record in country_records})
        expected = set(range(years[0], years[-1] + 1))
        missing = sorted(expected - set(years))
        if missing:
            report.add("missing_years", "warning", f"{iso3}: missing year(s) {_format_years(missing)}.", iso3=iso3)
        ordered = sorted(country_records, key=lambda record: record.year)
        for previous, current in zip(ordered, ordered[1:]):
            if current.year != previous.year + 1:
                continue
            a, b = previous.incidence_per_100k, current.incidence_per_100k
            if max(a, b) < DISCONTINUITY_MIN_RATE:
                continue
            ratio = (max(a, b) / min(a, b)) if min(a, b) > 0 else math.inf
            if ratio > DISCONTINUITY_RATIO:
                report.add(
                    "implausible_discontinuity",
                    "warning",
                    f"{iso3}: incidence changes by a factor of {ratio:.1f} between {previous.year} and {current.year}.",
                    iso3=iso3,
                    year=current.year,
                )
    return report


def validate_manifest(manifest: dict[str, Any]) -> ValidationReport:
    report = ValidationReport()
    missing = [key for key in MANIFEST_REQUIRED_FIELDS if key not in manifest]
    if missing:
        report.add("manifest_incomplete", "error", f"Manifest is missing: {', '.join(missing)}.")
        return report
    if manifest["schemaVersion"] != SNAPSHOT_SCHEMA_VERSION:
        report.add(
            "schema_changed",
            "error",
            f"Unsupported snapshot schema {manifest['schemaVersion']!r}; expected {SNAPSHOT_SCHEMA_VERSION!r}.",
        )
    if tuple(manifest["columns"]) != RECORD_COLUMNS:
        report.add("schema_changed", "error", "Manifest column list does not match the application schema.")
    data_file = manifest.get("dataFile") or {}
    for key in ("filename", "sha256", "rows"):
        if key not in data_file:
            report.add("manifest_incomplete", "error", f"Manifest dataFile is missing {key}.")
    return report


def check_vintage_compatibility(manifests: list[dict[str, Any]]) -> ValidationReport:
    """Series from different report vintages or schemas must not be combined silently.

    WHO revises its entire retrospective series each report round, so values
    for the same year differ between vintages.
    """
    report = ValidationReport()
    vintages = {(m.get("sourceDataset"), m.get("sourceReportYear")) for m in manifests}
    schemas = {m.get("schemaVersion") for m in manifests}
    transforms = {m.get("transformationVersion") for m in manifests}
    if len(vintages) > 1:
        described = ", ".join(sorted(f"{dataset} ({year})" for dataset, year in vintages))
        report.add("incompatible_vintage", "error", f"Snapshots come from different data vintages: {described}.")
    if len(schemas) > 1:
        report.add("schema_changed", "error", "Snapshots use different schema versions.")
    if len(transforms) > 1:
        report.add("incompatible_vintage", "warning", "Snapshots were produced by different transformation versions.")
    return report


def _parse_year(value: Any) -> int | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(str(value).strip())
    except ValueError:
        return None
    if not math.isfinite(number) or number != int(number):
        return None
    return int(number)


def _parse_number(value: Any) -> tuple[float | None, bool]:
    if value is None:
        return None, True
    if isinstance(value, bool):
        return None, False
    if isinstance(value, (int, float)):
        number = float(value)
        if math.isnan(number):
            return None, True
        return (number, True) if math.isfinite(number) else (None, False)
    text = str(value).strip()
    if text in {"", "NA", "NaN", "nan"}:
        return None, True
    try:
        number = float(text)
    except ValueError:
        return None, False
    return (number, True) if math.isfinite(number) else (None, False)


def _format_years(years: list[int]) -> str:
    if len(years) <= 6:
        return ", ".join(str(year) for year in years)
    return f"{years[0]}-{years[-1]} ({len(years)} years)"
