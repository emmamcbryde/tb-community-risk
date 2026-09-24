"""Read and write versioned, checksummed incidence snapshots (offline only)."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import io
import json
from pathlib import Path
from typing import Any

from engine.profiles.population_profile import IncidencePoint
from engine.who_incidence.schema import (
    RECORD_COLUMNS,
    SNAPSHOT_SCHEMA_VERSION,
    TRANSFORMATION_VERSION,
    IncidenceRecord,
    ValidationReport,
    parse_rows,
    validate_manifest,
    validate_records,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
BUNDLED_SNAPSHOT_DIR = REPO_ROOT / "data" / "who_incidence"
BUNDLED_SNAPSHOT_ID = "who-gtb2025-fixture-v1"
BUNDLED_MANIFEST = BUNDLED_SNAPSHOT_DIR / "fixture" / f"{BUNDLED_SNAPSHOT_ID}_manifest.json"


CRLF = bytes([13, 10])
LF = bytes([10])


class SnapshotError(ValueError):
    """Raised when a snapshot fails integrity or schema validation."""


@dataclass(frozen=True)
class IncidenceSnapshot:
    manifest: dict[str, Any]
    records: tuple[IncidenceRecord, ...]
    validation: dict[str, Any]

    @property
    def snapshot_id(self) -> str:
        return str(self.manifest["snapshotId"])

    @property
    def is_complete_dataset(self) -> bool:
        return bool(self.manifest.get("isCompleteDataset"))

    def countries(self) -> list[dict[str, str]]:
        names: dict[str, str] = {}
        for record in self.records:
            names.setdefault(record.iso3, record.country or record.iso3)
        return [{"iso3": iso3, "name": names[iso3]} for iso3 in sorted(names, key=lambda code: names[code])]

    def records_for(self, iso3: str) -> list[IncidenceRecord]:
        return sorted((record for record in self.records if record.iso3 == iso3), key=lambda record: record.year)

    def series_for(self, iso3: str) -> tuple[IncidencePoint, ...]:
        return tuple(
            IncidencePoint(
                year=record.year,
                estimate=record.incidence_per_100k,
                lower=record.incidence_per_100k_lo,
                upper=record.incidence_per_100k_hi,
            )
            for record in self.records_for(iso3)
        )

    def provenance_summary(self) -> dict[str, Any]:
        upstream = self.manifest.get("upstream") or {}
        return {
            "snapshotId": self.snapshot_id,
            "snapshotKind": self.manifest.get("snapshotKind"),
            "isCompleteDataset": self.is_complete_dataset,
            "sourceDataset": self.manifest.get("sourceDataset"),
            "sourceReportYear": self.manifest.get("sourceReportYear"),
            "sourceUrl": self.manifest.get("sourceUrl"),
            "upstreamRepository": upstream.get("repository"),
            "upstreamCommit": upstream.get("commit"),
            "extractionDate": self.manifest.get("extractionDate"),
            "transformationVersion": self.manifest.get("transformationVersion"),
            "dataSha256": (self.manifest.get("dataFile") or {}).get("sha256"),
            "validationStatus": self.manifest.get("validationStatus"),
            "citation": self.manifest.get("citation"),
            "licenceStatus": (self.manifest.get("licence") or {}).get("status"),
        }


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(Path(path).read_bytes())


def load_snapshot(manifest_path: Path | str) -> IncidenceSnapshot:
    """Load a snapshot from a local manifest, verifying checksum and schema."""
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    report = validate_manifest(manifest)
    if not report.is_valid:
        raise SnapshotError(_describe(report))
    data_path = manifest_path.parent / manifest["dataFile"]["filename"]
    payload = data_path.read_bytes()
    if CRLF in payload:  # tolerate a CRLF checkout of the LF-committed text file
        payload = payload.replace(CRLF, LF)
    actual = sha256_bytes(payload)
    if actual != manifest["dataFile"]["sha256"]:
        raise SnapshotError(f"Checksum mismatch for {data_path.name}: expected {manifest['dataFile']['sha256']}, found {actual}.")
    reader = csv.DictReader(io.StringIO(payload.decode("utf-8")))
    records, parse_report = parse_rows(reader)
    if tuple(reader.fieldnames or ()) != RECORD_COLUMNS:
        parse_report.add("schema_changed", "error", "Data file columns do not match the application schema.")
    if not parse_report.is_valid:
        raise SnapshotError(_describe(parse_report))
    if len(records) != int(manifest["dataFile"]["rows"]):
        raise SnapshotError("Row count does not match the manifest.")
    return IncidenceSnapshot(manifest=manifest, records=tuple(records), validation=parse_report.to_dict())


@lru_cache(maxsize=4)
def _load_cached(manifest_path: str, mtime: float) -> IncidenceSnapshot:
    return load_snapshot(manifest_path)


def load_bundled_snapshot(manifest_path: Path | str | None = None) -> IncidenceSnapshot:
    """Load the snapshot bundled with the application (no network access)."""
    path = Path(manifest_path or BUNDLED_MANIFEST)
    return _load_cached(str(path), path.stat().st_mtime)


def records_to_csv_bytes(records: list[IncidenceRecord]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(RECORD_COLUMNS), lineterminator="\n")
    writer.writeheader()
    for record in sorted(records, key=lambda item: (item.iso3, item.year)):
        writer.writerow({key: _csv_value(value) for key, value in record.to_row().items()})
    return buffer.getvalue().encode("utf-8")


def write_snapshot(
    records: list[IncidenceRecord],
    *,
    out_dir: Path,
    data_filename: str,
    manifest_filename: str,
    manifest_fields: dict[str, Any],
) -> dict[str, Any]:
    """Validate records and write a CSV plus manifest; refuses invalid data."""
    report = validate_records(records)
    if not report.is_valid:
        raise SnapshotError(_describe(report))
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = records_to_csv_bytes(records)
    (out_dir / data_filename).write_bytes(payload)
    years = sorted({record.year for record in records})
    manifest = {
        "schemaVersion": SNAPSHOT_SCHEMA_VERSION,
        "transformationVersion": TRANSFORMATION_VERSION,
        **manifest_fields,
        "columns": list(RECORD_COLUMNS),
        "dataFile": {"filename": data_filename, "sha256": sha256_bytes(payload), "rows": len(records)},
        "countries": sorted({record.iso3 for record in records}),
        "yearRange": [years[0], years[-1]] if years else None,
        "validationStatus": "passed_with_warnings" if report.warnings else "passed",
        "validationWarnings": [issue.message for issue in report.warnings],
    }
    text = json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    (out_dir / manifest_filename).write_text(text, encoding="utf-8", newline="\n")
    return manifest


def _csv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def _describe(report: ValidationReport) -> str:
    return "; ".join(issue.message for issue in report.errors[:10])
