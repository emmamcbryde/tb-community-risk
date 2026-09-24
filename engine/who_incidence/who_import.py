"""Maintainer-only import of the public WHO TB burden-estimates CSV.

Never called by ordinary application sessions. The importer reads local files
supplied by a maintainer, validates them, optionally cross-checks against a
locally exported copy of the Global TB Report estimates, and writes a
deterministic, checksummed snapshot plus manifest.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
import hashlib
import io
import json
import math
from pathlib import Path
from typing import Any, Iterable

from engine.who_incidence.contract import (
    COLUMNS,
    COUNT_TRIPLE,
    HIV_TRIPLE,
    IMPORTER_VERSION,
    ISO2_PATTERN,
    ISO3_PATTERN,
    NUMERIC_COLUMNS,
    RATE_TRIPLE,
    REQUIRED_WHO_COLUMNS,
    SNAPSHOT_CONTRACT_VERSION,
    WHO_COLUMN_MAP,
    WHO_REGIONS,
    SnapshotRow,
    field_definitions,
    rounding_half_unit,
)


WHO_ESTIMATES_URL = "https://extranet.who.int/tme/generateCSV.asp?ds=estimates"
WHO_DICTIONARY_URL = "https://extranet.who.int/tme/generateCSV.asp?ds=dictionary"
WHO_TB_DATA_PAGE = "https://www.who.int/teams/global-programme-on-tuberculosis-and-lung-health/data"
WHO_TERMS_URL = "https://www.who.int/about/policies/publishing/data-policy/terms-and-conditions"
GTB_REPORT_UPSTREAM = "https://github.com/GTB-TME/gtbreport2025"
FIRST_EXPECTED_YEAR = 2000
# Rate x population / 100 000 versus published count: tolerance in multiples of the
# combined published rounding half-units.
CONSISTENCY_WARNING_FACTOR = 1.0
CONSISTENCY_ERROR_FACTOR = 3.0
# Cross-check: |published - unrounded| must not exceed this multiple of the half-unit.
CROSSCHECK_TOLERANCE_FACTOR = 1.0 + 1e-9
CROSSCHECK_COLUMNS = {
    "incidence_per_100k": "inc",
    "incidence_per_100k_lo": "inc.lo",
    "incidence_per_100k_hi": "inc.hi",
    "incident_cases": "inc.num",
    "incident_cases_lo": "inc.lo.num",
    "incident_cases_hi": "inc.hi.num",
    "population": "pop",
}
# The WHO dictionary does not define the calendar-year key itself.
DICTIONARY_EXEMPT = {"year"}
CATEGORIES = ("fatal", "warning", "expected_missingness", "incomplete_series", "info")


class ImportFailed(ValueError):
    """Raised when fatal validation findings prevent writing a snapshot."""

    def __init__(self, report: "ImportReport") -> None:
        self.report = report
        first = "; ".join(item.message for item in report.fatal[:5])
        super().__init__(f"WHO import failed with {len(report.fatal)} fatal finding(s): {first}")


@dataclass(frozen=True)
class Finding:
    category: str
    code: str
    message: str
    iso3: str | None = None
    year: int | None = None

    def to_dict(self) -> dict[str, Any]:
        return {"category": self.category, "code": self.code, "message": self.message, "iso3": self.iso3, "year": self.year}


@dataclass
class ImportReport:
    findings: list[Finding] = field(default_factory=list)

    def add(self, category: str, code: str, message: str, *, iso3: str | None = None, year: int | None = None) -> None:
        if category not in CATEGORIES:
            raise ValueError(f"Unknown finding category {category!r}")
        self.findings.append(Finding(category, code, message, iso3, year))

    def of(self, category: str) -> list[Finding]:
        return [item for item in self.findings if item.category == category]

    @property
    def fatal(self) -> list[Finding]:
        return self.of("fatal")

    @property
    def is_valid(self) -> bool:
        return not self.fatal

    def codes(self, category: str | None = None) -> set[str]:
        return {item.code for item in self.findings if category is None or item.category == category}

    def counts(self) -> dict[str, int]:
        return {category: len(self.of(category)) for category in CATEGORIES}

    def to_dict(self) -> dict[str, Any]:
        return {"isValid": self.is_valid, "counts": self.counts(), "findings": [item.to_dict() for item in self.findings]}

    def summary_markdown(self, *, title: str = "WHO incidence import validation") -> str:
        lines = [f"# {title}", "", f"Result: {'PASSED' if self.is_valid else 'FAILED'}", ""]
        for category in CATEGORIES:
            items = self.of(category)
            lines.append(f"## {category.replace('_', ' ').capitalize()} ({len(items)})")
            if not items:
                lines.append("None.")
            grouped: dict[str, list[Finding]] = {}
            for item in items:
                grouped.setdefault(item.code, []).append(item)
            for code, group in sorted(grouped.items()):
                lines.append(f"- `{code}` x{len(group)}: {group[0].message}")
            lines.append("")
        return "\n".join(lines)


def decode_csv(payload: bytes, report: ImportReport, *, label: str) -> list[dict[str, str]] | None:
    if b"\x00" in payload:
        report.add("fatal", "corrupt_encoding", f"{label} contains NUL bytes; it is not a text CSV.")
        return None
    try:
        text = payload.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        report.add("fatal", "corrupt_encoding", f"{label} is not valid UTF-8 ({exc.reason} at byte {exc.start}).")
        return None
    try:
        reader = csv.DictReader(io.StringIO(text, newline=""))
        rows = list(reader)
    except csv.Error as exc:
        report.add("fatal", "corrupt_csv", f"{label} could not be parsed as CSV: {exc}.")
        return None
    if reader.fieldnames is None:
        report.add("fatal", "empty_file", f"{label} has no header row.")
        return None
    if any(None in row for row in rows):
        report.add("fatal", "corrupt_csv", f"{label} has rows with more fields than the header.")
        return None
    return rows


def parse_who_estimates(payload: bytes, *, report_year: int) -> tuple[list[SnapshotRow], ImportReport]:
    """Validate the WHO burden-estimates CSV and map it to the snapshot contract."""
    report = ImportReport()
    raw_rows = decode_csv(payload, report, label="WHO estimates CSV")
    if raw_rows is None:
        return [], report
    header = list(raw_rows[0].keys()) if raw_rows else []
    missing = [column for column in REQUIRED_WHO_COLUMNS if column not in header]
    if not raw_rows:
        report.add("fatal", "empty_file", "WHO estimates CSV has no data rows.")
        return [], report
    if missing:
        report.add("fatal", "schema_changed", f"Upstream schema changed; missing required column(s): {', '.join(missing)}.")
        return [], report
    extra = sorted(set(header) - set(REQUIRED_WHO_COLUMNS))
    report.add("info", "unused_columns", f"{len(extra)} published column(s) are not part of this snapshot contract.")

    rows: list[SnapshotRow] = []
    last_year = report_year - 1
    for line_number, raw in enumerate(raw_rows, start=2):
        values = {column: (raw.get(who) or "").strip() for who, column in WHO_COLUMN_MAP.items()}
        iso3 = values["iso3"]
        if not ISO3_PATTERN.match(iso3):
            report.add("fatal", "invalid_iso3", f"Line {line_number}: invalid or absent ISO3 code {iso3!r}.")
            continue
        if values["iso2"] and not ISO2_PATTERN.match(values["iso2"]):
            report.add("fatal", "invalid_iso2", f"Line {line_number}: invalid ISO2 code {values['iso2']!r}.", iso3=iso3)
        if values["who_region"] not in WHO_REGIONS:
            report.add("fatal", "invalid_region", f"Line {line_number}: unknown WHO region {values['who_region']!r}.", iso3=iso3)
        year = _parse_int(values["year"])
        if year is None:
            report.add("fatal", "invalid_year", f"Line {line_number}: invalid year {values['year']!r}.", iso3=iso3)
            continue
        if year > last_year:
            report.add(
                "fatal",
                "mixed_vintage",
                f"{iso3} {year}: year is later than {last_year}, the last estimate year of the {report_year} report round.",
                iso3=iso3,
                year=year,
            )
        elif year < FIRST_EXPECTED_YEAR:
            report.add("warning", "unexpected_year", f"{iso3} {year}: year precedes {FIRST_EXPECTED_YEAR}.", iso3=iso3, year=year)
        numbers: dict[str, float | None] = {}
        row_ok = True
        for column in NUMERIC_COLUMNS:
            parsed, ok = _parse_number(values[column])
            if not ok:
                report.add("fatal", "non_numeric", f"{iso3} {year}: non-numeric {column} {values[column]!r}.", iso3=iso3, year=year)
                row_ok = False
            elif parsed is not None and parsed < 0:
                report.add("fatal", "negative_value", f"{iso3} {year}: negative {column}.", iso3=iso3, year=year)
                row_ok = False
            numbers[column] = parsed
        if not row_ok:
            continue
        status = _triple_status(numbers, RATE_TRIPLE, iso3, year, report, "incidence rate")
        _triple_status(numbers, COUNT_TRIPLE, iso3, year, report, "incident cases")
        _triple_status(numbers, HIV_TRIPLE, iso3, year, report, "HIV-positive incidence", allow_partial=True)
        if numbers["population"] is None:
            report.add("fatal", "missing_population", f"{iso3} {year}: population is missing.", iso3=iso3, year=year)
        elif numbers["population"] == 0 and (numbers["incident_cases"] or 0) > 0:
            report.add("fatal", "zero_population", f"{iso3} {year}: zero population with non-zero incident cases.", iso3=iso3, year=year)
        _check_rate_count_consistency(values, numbers, iso3, year, report)
        values["year"] = str(year)
        values["incidence_status"] = "estimated" if status == "complete" else "not_estimated"
        rows.append(SnapshotRow(raw=tuple((column, values[column]) for column in COLUMNS)))

    _check_keys_and_mappings(rows, report)
    _check_coverage(rows, report, report_year=report_year)
    return sorted(rows, key=lambda row: (row.iso3, row.year)), report


def check_dictionary(payload: bytes, report: ImportReport) -> dict[str, str]:
    """Confirm every imported WHO variable is defined in the WHO data dictionary."""
    rows = decode_csv(payload, report, label="WHO data dictionary")
    if rows is None:
        return {}
    if not rows or not {"variable_name", "definition"} <= set(rows[0]):
        report.add("fatal", "dictionary_schema", "WHO data dictionary lacks variable_name/definition columns.")
        return {}
    definitions = {row["variable_name"].strip(): row["definition"].strip() for row in rows}
    for who in REQUIRED_WHO_COLUMNS:
        if who in DICTIONARY_EXEMPT:
            continue
        if who not in definitions:
            report.add("fatal", "dictionary_missing_variable", f"Variable {who} is not defined in the WHO data dictionary.")
    return {who: definitions[who] for who in REQUIRED_WHO_COLUMNS if who in definitions}


def load_crosscheck_export(payload: bytes, report: ImportReport) -> dict[tuple[str, int], dict[str, float | None]]:
    """Parse a CSV exported from est.rda (columns iso3, year, inc, inc.lo, ...)."""
    rows = decode_csv(payload, report, label="Cross-check export")
    if rows is None:
        return {}
    needed = {"iso3", "year", *CROSSCHECK_COLUMNS.values()}
    header = set(rows[0].keys()) if rows else set()
    if not needed <= header:
        report.add("fatal", "crosscheck_schema", f"Cross-check export is missing: {', '.join(sorted(needed - header))}.")
        return {}
    out: dict[tuple[str, int], dict[str, float | None]] = {}
    for row in rows:
        key = (row["iso3"].strip(), int(float(row["year"])))
        out[key] = {column: _parse_number(row[source])[0] for column, source in CROSSCHECK_COLUMNS.items()}
    return out


def cross_check(rows: list[SnapshotRow], reference: dict[tuple[str, int], dict[str, float | None]], report: ImportReport) -> dict[str, Any]:
    """Compare published values against unrounded Global TB Report estimates."""
    keys = {(row.iso3, row.year) for row in rows}
    missing_in_reference = sorted(keys - set(reference))
    extra_in_reference = sorted(set(reference) - keys)
    if missing_in_reference or extra_in_reference:
        report.add(
            "fatal",
            "crosscheck_coverage",
            f"Country-year keys differ: {len(missing_in_reference)} only in WHO CSV, {len(extra_in_reference)} only in cross-check data.",
        )
    compared = disagreements = withheld = 0
    worst: list[dict[str, Any]] = []
    for row in rows:
        ref = reference.get((row.iso3, row.year))
        if ref is None:
            continue
        for column in CROSSCHECK_COLUMNS:
            text = row.text(column)
            ref_value = ref[column]
            if text == "" and (ref_value is None or (isinstance(ref_value, float) and math.isnan(ref_value))):
                continue
            compared += 1
            if text == "":
                withheld += 1
                report.add(
                    "warning",
                    "crosscheck_unpublished",
                    f"{row.iso3} {row.year} {column}: present in the report-round analysis output but not published in the public WHO dataset; the public (blank) value is kept.",
                    iso3=row.iso3,
                    year=row.year,
                )
                continue
            if ref_value is None:
                disagreements += 1
                report.add("fatal", "crosscheck_disagreement", f"{row.iso3} {row.year} {column}: published by WHO but absent from the cross-check data.", iso3=row.iso3, year=row.year)
                continue
            published = float(text)
            allowed = rounding_half_unit(text) * CROSSCHECK_TOLERANCE_FACTOR + 1e-9 * abs(ref_value)
            difference = abs(published - ref_value)
            if difference > allowed:
                disagreements += 1
                report.add(
                    "fatal",
                    "crosscheck_disagreement",
                    f"{row.iso3} {row.year} {column}: published {text} vs unrounded {ref_value:.6g} (tolerance {allowed:.3g}).",
                    iso3=row.iso3,
                    year=row.year,
                )
                worst.append({"iso3": row.iso3, "year": row.year, "column": column, "published": text, "reference": ref_value})
    if compared and not disagreements:
        report.add(
            "info",
            "crosscheck_passed",
            f"All {compared - withheld} values published in both sources agree within published rounding precision.",
        )
    return {
        "comparedValues": compared,
        "disagreements": disagreements,
        "withheldFromPublicDataset": withheld,
        "tolerance": "absolute difference <= half the unit of the last published digit",
        "examples": worst[:20],
    }


def snapshot_bytes(rows: list[SnapshotRow]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(COLUMNS)
    for row in sorted(rows, key=lambda item: (item.iso3, item.year)):
        writer.writerow([value for _, value in row.raw])
    return buffer.getvalue().encode("utf-8")


def build_manifest(
    rows: list[SnapshotRow],
    report: ImportReport,
    *,
    data_bytes: bytes,
    report_year: int,
    source_files: list[dict[str, Any]],
    crosscheck: dict[str, Any] | None,
    access_date: str,
    importer_commit: str | None,
    kind: str,
    subset_iso3: Iterable[str] | None = None,
) -> dict[str, Any]:
    data_sha = hashlib.sha256(data_bytes).hexdigest()
    countries = sorted({row.iso3 for row in rows})
    estimated = sorted({row.iso3 for row in rows if row.is_estimated})
    years = sorted({row.year for row in rows})
    by_country: dict[str, list[SnapshotRow]] = {}
    for row in rows:
        by_country.setdefault(row.iso3, []).append(row)
    exclusions = [
        {"iso3": iso3, "reason": "No incidence estimates published in this report round."}
        for iso3, group in sorted(by_country.items())
        if not any(item.is_estimated for item in group)
    ]
    snapshot_id = f"who-gtb{report_year}-incidence-{'fixture-' if kind == 'test_fixture' else ''}{data_sha[:12]}"
    manifest = {
        "contractVersion": SNAPSHOT_CONTRACT_VERSION,
        "snapshotId": snapshot_id,
        "snapshotKind": kind,
        "isCompleteDataset": kind == "production" and subset_iso3 is None,
        "subsetIso3": sorted(subset_iso3) if subset_iso3 is not None else None,
        "reportYear": report_year,
        "seriesDescription": (
            f"Retrospective country series from the {report_year} WHO Global Tuberculosis Report round "
            f"({years[0]}-{years[-1]}). Rounds are never combined."
        ),
        "source": {
            "publisher": "World Health Organization",
            "dataset": "WHO TB burden estimates (country level)",
            "url": WHO_ESTIMATES_URL,
            "dictionaryUrl": WHO_DICTIONARY_URL,
            "dataPage": WHO_TB_DATA_PAGE,
        },
        "sourceFiles": source_files,
        "crossCheck": crosscheck,
        "importerVersion": IMPORTER_VERSION,
        "importerCommit": importer_commit,
        "transformation": (
            "Selected published columns renamed to the snapshot contract; values copied as published text without "
            "re-rounding; one derived column (incidence_status); rows sorted by ISO3 and year; no values imputed, "
            "smoothed or combined across report rounds."
        ),
        "dataFile": {"filename": f"{snapshot_id}.csv", "sha256": data_sha, "rows": len(rows), "bytes": len(data_bytes)},
        "coverage": {
            "countriesAndAreas": len(countries),
            "withEstimates": len(estimated),
            "yearRange": [years[0], years[-1]] if years else None,
            "incompleteSeries": sorted(
                {item.iso3 for item in report.of("incomplete_series") if item.iso3}
            ),
        },
        "fields": field_definitions(),
        "citation": (
            f"World Health Organization. WHO TB burden estimates (country level), Global Tuberculosis Report {report_year} "
            f"data. Geneva: WHO; {report_year}. Accessed {access_date}. Estimates are based on data reported by "
            "countries and areas to WHO."
        ),
        "terms": {
            "status": "reuse_permitted_with_conditions",
            "url": WHO_TERMS_URL,
            "summary": (
                "WHO dataset terms grant a royalty-free, non-exclusive right to use, reproduce, extract, copy, "
                "distribute and include the datasets for public health purposes, with attribution; no commercial "
                "promotion; no implied WHO endorsement; no use of the WHO name or emblem beyond attribution; "
                "only minimal alteration without prior written WHO authorization."
            ),
            "noEndorsement": "WHO has not reviewed or endorsed this application or its use of the data.",
            "openQuestion": (
                "Whether renaming columns and storing a selected subset of published values counts as 'minimal "
                "alteration' under the WHO dataset terms; values themselves are unaltered."
            ),
        },
        "validationSummary": report.counts(),
        "knownExclusions": exclusions,
        "volatile": {"accessDate": access_date},
    }
    manifest["analyticalHash"] = analytical_hash(manifest)
    return manifest


def analytical_hash(manifest: dict[str, Any]) -> str:
    """Hash of the manifest excluding volatile metadata (access date)."""
    payload = {key: value for key, value in manifest.items() if key not in {"volatile", "analyticalHash", "citation"}}
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def write_snapshot_files(out_dir: Path, data_bytes: bytes, manifest: dict[str, Any], report: ImportReport) -> dict[str, Path]:
    if not report.is_valid:
        raise ImportFailed(report)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = manifest["snapshotId"]
    paths = {
        "data": out_dir / f"{stem}.csv",
        "manifest": out_dir / f"{stem}_manifest.json",
        "validation": out_dir / f"{stem}_validation.json",
        "summary": out_dir / f"{stem}_validation.md",
    }
    paths["data"].write_bytes(data_bytes)
    paths["manifest"].write_text(json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    paths["validation"].write_text(json.dumps(report.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    paths["summary"].write_text(report.summary_markdown().rstrip() + "\n", encoding="utf-8", newline="\n")
    return paths


def _triple_status(
    numbers: dict[str, float | None],
    triple: tuple[str, str, str],
    iso3: str,
    year: int,
    report: ImportReport,
    label: str,
    *,
    allow_partial: bool = False,
) -> str:
    estimate, lower, upper = (numbers[column] for column in triple)
    present = [value is not None for value in (estimate, lower, upper)]
    if not any(present):
        return "absent"
    if not all(present):
        if allow_partial and estimate is not None:
            report.add("warning", "incomplete_bounds", f"{iso3} {year}: {label} published without complete bounds.", iso3=iso3, year=year)
            return "partial"
        report.add("fatal", "incomplete_triple", f"{iso3} {year}: {label} estimate and bounds are only partly published.", iso3=iso3, year=year)
        return "partial"
    if lower > estimate:
        report.add("fatal", "lower_exceeds_estimate", f"{iso3} {year}: {label} lower bound exceeds the estimate.", iso3=iso3, year=year)
    if upper < estimate:
        report.add("fatal", "upper_below_estimate", f"{iso3} {year}: {label} upper bound is below the estimate.", iso3=iso3, year=year)
    return "complete"


def _check_rate_count_consistency(values: dict[str, str], numbers: dict[str, float | None], iso3: str, year: int, report: ImportReport) -> None:
    rate, cases, population = numbers["incidence_per_100k"], numbers["incident_cases"], numbers["population"]
    if rate is None or cases is None or not population:
        return
    implied = rate * population / 100_000
    slack = rounding_half_unit(values["incidence_per_100k"]) * population / 100_000 + rounding_half_unit(values["incident_cases"])
    difference = abs(implied - cases)
    if difference > slack * CONSISTENCY_ERROR_FACTOR + 1e-9:
        report.add(
            "fatal",
            "rate_count_inconsistent",
            f"{iso3} {year}: rate x population implies {implied:.4g} cases but {values['incident_cases']} are published.",
            iso3=iso3,
            year=year,
        )
    elif difference > slack * CONSISTENCY_WARNING_FACTOR + 1e-9:
        report.add(
            "warning",
            "rate_count_rounding",
            f"{iso3} {year}: rate x population ({implied:.4g}) and published cases ({values['incident_cases']}) differ by more than published rounding.",
            iso3=iso3,
            year=year,
        )


def _check_keys_and_mappings(rows: list[SnapshotRow], report: ImportReport) -> None:
    seen: set[tuple[str, int]] = set()
    identity: dict[str, tuple[str, str, str, str]] = {}
    names: dict[str, str] = {}
    iso2s: dict[str, str] = {}
    for row in rows:
        key = (row.iso3, row.year)
        if key in seen:
            report.add("fatal", "duplicate_row", f"Duplicate row for {row.iso3} {row.year}.", iso3=row.iso3, year=row.year)
        seen.add(key)
        ident = (row.text("country"), row.text("iso2"), row.text("iso_numeric"), row.text("who_region"))
        if identity.setdefault(row.iso3, ident) != ident:
            report.add("fatal", "inconsistent_identity", f"{row.iso3}: name, codes or region differ between years.", iso3=row.iso3)
        if names.setdefault(row.text("country"), row.iso3) != row.iso3:
            report.add("fatal", "duplicate_iso3_mapping", f"Country name {row.text('country')!r} maps to more than one ISO3 code.", iso3=row.iso3)
        if row.text("iso2") and iso2s.setdefault(row.text("iso2"), row.iso3) != row.iso3:
            report.add("fatal", "duplicate_iso3_mapping", f"ISO2 {row.text('iso2')!r} maps to more than one ISO3 code.", iso3=row.iso3)


def _check_coverage(rows: list[SnapshotRow], report: ImportReport, *, report_year: int) -> None:
    by_country: dict[str, list[SnapshotRow]] = {}
    for row in rows:
        by_country.setdefault(row.iso3, []).append(row)
    if not by_country:
        return
    first = min(row.year for row in rows)
    last = max(row.year for row in rows)
    if last != report_year - 1:
        report.add("fatal", "unexpected_year_range", f"Latest year is {last}; the {report_year} round should end in {report_year - 1}.")
    for iso3, group in sorted(by_country.items()):
        estimated_years = sorted(row.year for row in group if row.is_estimated)
        all_years = sorted(row.year for row in group)
        if not estimated_years:
            report.add("expected_missingness", "no_estimates", f"{iso3}: population published but no incidence estimates.", iso3=iso3)
            continue
        gaps = sorted(set(range(estimated_years[0], estimated_years[-1] + 1)) - set(estimated_years))
        if gaps:
            report.add("warning", "internal_gap", f"{iso3}: estimates missing inside the series for {gaps}.", iso3=iso3)
        if all_years[0] > first or all_years[-1] < last:
            report.add(
                "incomplete_series",
                "short_series",
                f"{iso3}: series covers {all_years[0]}-{all_years[-1]} rather than {first}-{last}.",
                iso3=iso3,
            )


def _parse_int(text: str) -> int | None:
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value) or value != int(value):
        return None
    return int(value)


def _parse_number(text: Any) -> tuple[float | None, bool]:
    if text is None:
        return None, True
    text = str(text).strip()
    if text in {"", "NA"}:
        return None, True
    try:
        value = float(text)
    except ValueError:
        return None, False
    if math.isnan(value):
        return None, True
    return (value, True) if math.isfinite(value) else (None, False)
