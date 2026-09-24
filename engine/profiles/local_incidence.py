"""Validated local or subnational incidence uploads.

Uploaded data stay in the user's session; nothing is sent to external services.
An upload is always recorded as a local file, never as WHO data, and the
``measure`` column must declare estimated incidence so that notification counts
are not silently treated as incidence.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field, replace
import hashlib
import io
import json
import math
import re
from typing import Any

from engine.profiles.country import INCIDENCE_LINK_NOTE, ConflictRequiresChoice, KEEP_CURRENT
from engine.profiles.population_profile import (
    IncidenceData,
    IncidencePoint,
    Location,
    LocationKind,
    PopulationProfile,
    Provenance,
)


TEMPLATE_COLUMNS = ("location", "iso3", "year", "measure", "incidence_per_100k", "lower", "upper", "population", "source", "notes")
REQUIRED_COLUMNS = ("location", "year", "measure", "incidence_per_100k", "source")
ACCEPTED_MEASURES = {"estimated_incidence"}
REJECTED_MEASURES = {"notification", "notifications", "notified", "notification_rate", "case_notification"}
ISO3_PATTERN = re.compile(r"^[A-Z]{3}$")
MAX_BYTES = 2_000_000
TEMPLATE_CSV = (
    ",".join(TEMPLATE_COLUMNS)
    + "\n"
    + "Example district,,2022,estimated_incidence,48,39,58,250000,Example local estimate (replace),Illustrative row - replace\n"
    + "Example district,,2023,estimated_incidence,45,36,55,252000,Example local estimate (replace),Illustrative row - replace\n"
)
FIELD_GUIDE = {
    "location": "Name of the country, district or area (required, one location per file).",
    "iso3": "Optional ISO3 country code; leave blank for subnational areas.",
    "year": "Calendar year (required, one row per year).",
    "measure": "Must be 'estimated_incidence'. Notification rates are not accepted as incidence.",
    "incidence_per_100k": "Estimated TB disease incidence per 100,000 population per year (required).",
    "lower": "Optional lower uncertainty bound (per 100,000). Give both bounds or neither.",
    "upper": "Optional upper uncertainty bound (per 100,000).",
    "population": "Optional population of the area for that year.",
    "source": "Citation or description of where the estimates come from (required).",
    "notes": "Optional notes.",
}


@dataclass
class UploadResult:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    location: str = ""
    iso3: str | None = None
    points: tuple[IncidencePoint, ...] = ()
    sources: tuple[str, ...] = ()
    populations: tuple[tuple[int, float], ...] = ()
    file_sha256: str = ""
    filename: str = ""

    @property
    def is_valid(self) -> bool:
        return not self.errors

    @property
    def bounds_available(self) -> bool:
        return bool(self.points) and all(p.lower is not None and p.upper is not None for p in self.points)

    @property
    def content_hash(self) -> str:
        payload = {
            "location": self.location,
            "iso3": self.iso3,
            "points": [p.to_dict() for p in self.points],
            "sources": list(self.sources),
            "populations": [list(item) for item in self.populations],
        }
        return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()

    def preview_rows(self) -> list[dict[str, Any]]:
        populations = dict(self.populations)
        return [
            {
                "Year": p.year,
                "Estimate": p.estimate,
                "Lower": p.lower,
                "Upper": p.upper,
                "Population": populations.get(p.year),
            }
            for p in self.points
        ]


def parse_local_incidence(payload: bytes, *, filename: str = "upload.csv") -> UploadResult:
    result = UploadResult(filename=filename, file_sha256=hashlib.sha256(payload).hexdigest())
    if len(payload) > MAX_BYTES:
        result.errors.append("The file is larger than 2 MB; a single-location incidence file should be much smaller.")
        return result
    try:
        text = payload.decode("utf-8-sig")
    except UnicodeDecodeError:
        result.errors.append("The file is not UTF-8 encoded text. Save it as CSV (UTF-8).")
        return result
    try:
        reader = csv.DictReader(io.StringIO(text, newline=""))
        rows = list(reader)
    except csv.Error as exc:
        result.errors.append(f"The file could not be read as CSV: {exc}.")
        return result
    header = [column.strip() for column in (reader.fieldnames or [])]
    missing = [column for column in REQUIRED_COLUMNS if column not in header]
    if missing:
        result.errors.append(f"Missing required column(s): {', '.join(missing)}. Download the template for the expected layout.")
        return result
    unknown = [column for column in header if column not in TEMPLATE_COLUMNS]
    if unknown:
        result.warnings.append(f"Ignored unrecognised column(s): {', '.join(unknown)}.")
    if not rows:
        result.errors.append("The file contains no data rows.")
        return result

    locations = {(row.get("location") or "").strip() for row in rows}
    if "" in locations:
        result.errors.append("Every row needs a location.")
    if len(locations - {""}) > 1:
        result.errors.append("The file contains more than one location; upload one location per file.")
    result.location = next(iter(sorted(locations - {""})), "")
    codes = {(row.get("iso3") or "").strip().upper() for row in rows}
    if len(codes) > 1:
        result.errors.append("ISO3 codes differ between rows.")
    code = next(iter(codes), "")
    if code and not ISO3_PATTERN.match(code):
        result.errors.append(f"Invalid ISO3 code {code!r}.")
    result.iso3 = code or None

    points: list[IncidencePoint] = []
    populations: list[tuple[int, float]] = []
    sources: set[str] = set()
    seen: set[int] = set()
    for line, row in enumerate(rows, start=2):
        measure = (row.get("measure") or "").strip().lower()
        if measure in REJECTED_MEASURES:
            result.errors.append(f"Line {line}: notification data cannot be used as estimated incidence.")
            continue
        if measure not in ACCEPTED_MEASURES:
            result.errors.append(f"Line {line}: measure must be 'estimated_incidence', not {measure!r}.")
            continue
        year = _int(row.get("year"))
        if year is None or not 1900 <= year <= 2100:
            result.errors.append(f"Line {line}: invalid year {row.get('year')!r}.")
            continue
        if year in seen:
            result.errors.append(f"Line {line}: duplicate year {year}.")
            continue
        seen.add(year)
        values = {}
        for column in ("incidence_per_100k", "lower", "upper", "population"):
            value, ok = _number(row.get(column))
            if not ok:
                result.errors.append(f"Line {line}: {column} is not a number ({row.get(column)!r}).")
            elif value is not None and value < 0:
                result.errors.append(f"Line {line}: {column} is negative.")
            values[column] = value
        estimate, lower, upper = values["incidence_per_100k"], values["lower"], values["upper"]
        if estimate is None:
            result.errors.append(f"Line {line}: incidence_per_100k is required.")
            continue
        if (lower is None) != (upper is None):
            result.errors.append(f"Line {line}: give both lower and upper bounds, or neither.")
            continue
        if lower is not None and lower > estimate:
            result.errors.append(f"Line {line}: lower bound exceeds the estimate.")
        if upper is not None and upper < estimate:
            result.errors.append(f"Line {line}: upper bound is below the estimate.")
        source = (row.get("source") or "").strip()
        if not source:
            result.errors.append(f"Line {line}: source is required.")
        sources.add(source)
        points.append(IncidencePoint(year=year, estimate=estimate, lower=lower, upper=upper))
        if values["population"] is not None:
            populations.append((year, values["population"]))
    result.points = tuple(sorted(points, key=lambda p: p.year))
    result.populations = tuple(sorted(populations))
    result.sources = tuple(sorted(sources - {""}))
    if points and not result.bounds_available:
        result.warnings.append("Uncertainty bounds are missing for some years; propagated trend uncertainty will be unavailable.")
    if any("who" in source.lower() for source in result.sources):
        result.warnings.append(
            "The file cites WHO. It is still treated as a local file, not as the checksummed WHO snapshot."
        )
    return result


def apply_local_incidence(
    profile: PopulationProfile,
    upload: UploadResult,
    *,
    resolutions: dict[str, str] | None = None,
) -> PopulationProfile:
    if not upload.is_valid:
        raise ValueError("Cannot apply an invalid upload.")
    resolutions = resolutions or {}
    replacing_who = profile.incidence.provenance is Provenance.WHO_SNAPSHOT
    replacing_local = profile.incidence.provenance is Provenance.LOCAL_UPLOAD
    if (replacing_who or replacing_local) and "incidence" not in resolutions:
        raise ConflictRequiresChoice(["Estimated TB incidence series"])
    if resolutions.get("incidence") == KEEP_CURRENT:
        return profile
    kind = LocationKind.SUBNATIONAL if not upload.iso3 else LocationKind.COUNTRY
    latest_population = upload.populations[-1] if upload.populations else None
    return replace(
        profile,
        location=Location(
            name=upload.location,
            kind=kind,
            iso3=upload.iso3,
            national_population=None if latest_population is None else latest_population[1],
            national_population_year=None if latest_population is None else latest_population[0],
            national_population_source="" if latest_population is None else f"Local file {upload.filename}",
        ),
        incidence=IncidenceData(
            source=f"Local file: {upload.filename}",
            snapshot_id=None,
            provenance=Provenance.LOCAL_UPLOAD,
            series=upload.points,
            notes=INCIDENCE_LINK_NOTE,
            source_detail=tuple(
                sorted(
                    {
                        "filename": upload.filename,
                        "fileSha256": upload.file_sha256,
                        "contentHash": upload.content_hash,
                        "citedSources": tuple(upload.sources),
                        "boundsAvailable": upload.bounds_available,
                    }.items()
                )
            ),
        ),
        data_vintage=f"Local incidence file {upload.filename} (sha256 {upload.file_sha256[:12]})",
    )


def _int(value: Any) -> int | None:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return int(number) if math.isfinite(number) and number == int(number) else None


def _number(value: Any) -> tuple[float | None, bool]:
    text = "" if value is None else str(value).strip()
    if text == "":
        return None, True
    try:
        number = float(text)
    except ValueError:
        return None, False
    return (number, True) if math.isfinite(number) else (None, False)
