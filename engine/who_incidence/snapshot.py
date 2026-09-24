"""Offline loading of versioned, checksummed WHO incidence snapshots (contract v2).

Ordinary application sessions only read local files. A production snapshot in
``data/who_incidence/snapshots/`` is preferred; the small test fixture in
``data/who_incidence/fixture/`` is used when no production snapshot is present.
"""

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
from engine.who_incidence.contract import (
    COLUMNS,
    MANIFEST_REQUIRED_FIELDS,
    SNAPSHOT_CONTRACT_VERSION,
    WHO_COLUMN_MAP,
    SnapshotRow,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT_ROOT = REPO_ROOT / "data" / "who_incidence"
PRODUCTION_DIR = SNAPSHOT_ROOT / "snapshots"
FIXTURE_DIR = SNAPSHOT_ROOT / "fixture"
UPDATE_GUIDE = "docs/who_data_update_guide.md"
CRLF = bytes([13, 10])
LF = bytes([10])


class SnapshotError(ValueError):
    """Raised when a snapshot fails integrity or contract validation."""


class SnapshotUnavailable(SnapshotError):
    """Raised when no snapshot is installed; the message is maintainer-facing."""


@dataclass(frozen=True)
class IncidenceSnapshot:
    manifest: dict[str, Any]
    rows: tuple[SnapshotRow, ...]
    manifest_path: str

    @property
    def snapshot_id(self) -> str:
        return str(self.manifest["snapshotId"])

    @property
    def report_year(self) -> int:
        return int(self.manifest["reportYear"])

    @property
    def is_complete_dataset(self) -> bool:
        return bool(self.manifest.get("isCompleteDataset"))

    @property
    def data_sha256(self) -> str:
        return str(self.manifest["dataFile"]["sha256"])

    def countries(self) -> list[dict[str, Any]]:
        info: dict[str, dict[str, Any]] = {}
        for row in self.rows:
            entry = info.setdefault(
                row.iso3,
                {"iso3": row.iso3, "name": row.text("country"), "region": row.text("who_region"), "estimatedYears": 0},
            )
            entry["estimatedYears"] += 1 if row.is_estimated else 0
        return sorted(info.values(), key=lambda item: item["name"])

    def rows_for(self, iso3: str) -> list[SnapshotRow]:
        return sorted((row for row in self.rows if row.iso3 == iso3), key=lambda row: row.year)

    def series_for(self, iso3: str) -> tuple[IncidencePoint, ...]:
        """Estimated years only; unpublished years are omitted, never imputed."""
        return tuple(
            IncidencePoint(
                year=row.year,
                estimate=row.number("incidence_per_100k"),
                lower=row.number("incidence_per_100k_lo"),
                upper=row.number("incidence_per_100k_hi"),
            )
            for row in self.rows_for(iso3)
            if row.is_estimated
        )

    def country_summary(self, iso3: str) -> dict[str, Any]:
        rows = self.rows_for(iso3)
        if not rows:
            raise KeyError(f"{iso3} is not in snapshot {self.snapshot_id}.")
        estimated = [row for row in rows if row.is_estimated]
        latest = estimated[-1] if estimated else None
        return {
            "iso3": iso3,
            "name": rows[0].text("country"),
            "region": rows[0].text("who_region"),
            "years": [rows[0].year, rows[-1].year],
            "estimatedYears": [row.year for row in estimated],
            "latest": None
            if latest is None
            else {
                "year": latest.year,
                "estimate": latest.number("incidence_per_100k"),
                "lower": latest.number("incidence_per_100k_lo"),
                "upper": latest.number("incidence_per_100k_hi"),
                "cases": latest.number("incident_cases"),
                "population": latest.number("population"),
            },
        }

    def provenance_summary(self) -> dict[str, Any]:
        source = self.manifest.get("source") or {}
        cross = self.manifest.get("crossCheck") or {}
        return {
            "snapshotId": self.snapshot_id,
            "snapshotKind": self.manifest.get("snapshotKind"),
            "isCompleteDataset": self.is_complete_dataset,
            "reportYear": self.report_year,
            "dataset": source.get("dataset"),
            "sourceUrl": source.get("url"),
            "accessDate": (self.manifest.get("volatile") or {}).get("accessDate"),
            "crossCheckRepository": cross.get("repository"),
            "crossCheckCommit": cross.get("commit"),
            "importerVersion": self.manifest.get("importerVersion"),
            "dataSha256": self.data_sha256,
            "analyticalHash": self.manifest.get("analyticalHash"),
            "coverage": self.manifest.get("coverage"),
            "citation": self.manifest.get("citation"),
            "terms": self.manifest.get("terms"),
        }


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(Path(path).read_bytes())


def validate_manifest(manifest: Any) -> list[str]:
    if not isinstance(manifest, dict):
        return ["Manifest must be a JSON object."]
    problems = [f"Manifest is missing '{key}'." for key in MANIFEST_REQUIRED_FIELDS if key not in manifest]
    if problems:
        return problems
    if manifest["contractVersion"] != SNAPSHOT_CONTRACT_VERSION:
        problems.append(f"Unsupported snapshot contract {manifest['contractVersion']!r}; expected {SNAPSHOT_CONTRACT_VERSION!r}.")
    data_file = manifest.get("dataFile")
    if not isinstance(data_file, dict) or not {"filename", "sha256", "rows"} <= set(data_file):
        problems.append("Manifest dataFile must record filename, sha256 and rows.")
    if [item.get("name") for item in manifest.get("fields") or []] != list(COLUMNS):
        problems.append("Manifest field list does not match the snapshot contract.")
    return problems


def load_snapshot(manifest_path: Path | str) -> IncidenceSnapshot:
    """Load a snapshot, verifying manifest, checksum, columns and every row."""
    from engine.who_incidence.who_import import parse_who_estimates

    manifest_path = Path(manifest_path)
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SnapshotError(f"Cannot read snapshot manifest {manifest_path.name}: {exc}") from exc
    problems = validate_manifest(manifest)
    if problems:
        raise SnapshotError(f"Invalid snapshot manifest {manifest_path.name}: " + " ".join(problems))
    data_path = manifest_path.parent / manifest["dataFile"]["filename"]
    try:
        payload = data_path.read_bytes()
    except OSError as exc:
        raise SnapshotError(f"Snapshot data file {data_path.name} is missing: {exc}") from exc
    if CRLF in payload:  # tolerate a CRLF checkout of the LF-committed text file
        payload = payload.replace(CRLF, LF)
    actual = sha256_bytes(payload)
    if actual != manifest["dataFile"]["sha256"]:
        raise SnapshotError(
            f"Checksum mismatch for {data_path.name}: manifest records {manifest['dataFile']['sha256']}, file has {actual}. "
            "The snapshot may be corrupted or edited; re-run the importer."
        )
    rows = _read_rows(payload)
    if len(rows) != int(manifest["dataFile"]["rows"]):
        raise SnapshotError("Snapshot row count does not match the manifest.")
    # Re-validate by mapping back to the published column names and applying the importer's rules.
    reparsed, report = parse_who_estimates(_as_who_csv(rows), report_year=int(manifest["reportYear"]))
    if not report.is_valid:
        raise SnapshotError("Snapshot failed validation: " + "; ".join(item.message for item in report.fatal[:5]))
    if tuple(reparsed) != tuple(rows):
        raise SnapshotError("Snapshot rows are not in canonical form.")
    return IncidenceSnapshot(manifest=manifest, rows=tuple(rows), manifest_path=str(manifest_path))


def find_manifests(directory: Path) -> list[Path]:
    if not directory.is_dir():
        return []
    return sorted(directory.glob("*_manifest.json"))


def default_manifest_path() -> Path:
    """Production snapshot if installed, otherwise the bundled test fixture."""
    for directory in (PRODUCTION_DIR, FIXTURE_DIR):
        manifests = find_manifests(directory)
        if len(manifests) > 1:
            raise SnapshotError(
                f"More than one snapshot is installed in {directory.relative_to(REPO_ROOT).as_posix()}; "
                "keep exactly one report round per directory."
            )
        if manifests:
            return manifests[0]
    raise SnapshotUnavailable(
        "No WHO incidence snapshot is installed under data/who_incidence/. "
        f"Build one with scripts/import_who_incidence.py (see {UPDATE_GUIDE})."
    )


@lru_cache(maxsize=4)
def _load_cached(manifest_path: str, mtime: float, size: int) -> IncidenceSnapshot:
    return load_snapshot(manifest_path)


def load_bundled_snapshot(manifest_path: Path | str | None = None) -> IncidenceSnapshot:
    """Load the installed snapshot (no network access)."""
    path = Path(manifest_path) if manifest_path else default_manifest_path()
    if not path.is_file():
        raise SnapshotUnavailable(f"Snapshot manifest {path.name} was not found. See {UPDATE_GUIDE}.")
    stat = path.stat()
    return _load_cached(str(path), stat.st_mtime, stat.st_size)


def _read_rows(payload: bytes) -> list[SnapshotRow]:
    reader = csv.reader(io.StringIO(payload.decode("utf-8"), newline=""))
    header = next(reader, None)
    if tuple(header or ()) != COLUMNS:
        raise SnapshotError("Snapshot data columns do not match the snapshot contract.")
    rows = []
    for values in reader:
        if len(values) != len(COLUMNS):
            raise SnapshotError("Snapshot data file has a malformed row.")
        rows.append(SnapshotRow(raw=tuple(zip(COLUMNS, values))))
    return rows


def _as_who_csv(rows: list[SnapshotRow]) -> bytes:
    inverse = {column: who for who, column in WHO_COLUMN_MAP.items()}
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    who_columns = [inverse[column] for column in COLUMNS if column in inverse]
    writer.writerow(who_columns)
    for row in rows:
        values = dict(row.raw)
        writer.writerow([values[WHO_COLUMN_MAP[who]] for who in who_columns])
    return buffer.getvalue().encode("utf-8")
