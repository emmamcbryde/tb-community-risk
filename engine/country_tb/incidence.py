"""Strict frozen WHO loader and versioned national/local incidence profiles."""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal
import json
from pathlib import Path
import re

from engine.country_tb.data_common import (
    DataValidationError, LOADER_VERSION, PROFILE_VERSION, REPO_ROOT, canonical,
    contained_file, csv_rows, decimal, digest, file_bytes, numeric_text, require_hash, seal,
    publish_immutable,
)
from engine.who_incidence.contract import COLUMNS

WHO_2025 = "who-gtb2025-incidence-6997e3b087d2"
WHO_ROOT = REPO_ROOT / "data/who_incidence/snapshots"
RATE_UNIT = "cases_per_100000_person_year"
MODES = {"national_who_snapshot", "subnational_scaled", "subnational_user_series", "custom_population_series"}


def _manifest_hash(manifest):
    return digest(canonical({k: v for k, v in manifest.items() if k not in {"volatile", "analyticalHash", "citation"}}))


def _triple(row, columns):
    values = [row[column] for column in columns]
    if all(value == "" for value in values):
        return None
    if any(value == "" for value in values):
        raise DataValidationError(f"Partial uncertainty triple: {columns}")
    central, lower, upper = map(decimal, values)
    if not lower <= central <= upper:
        raise DataValidationError("Uncertainty ordering must be lower <= central <= upper.")
    return [numeric_text(value) for value in (central, lower, upper)]


def load_who_snapshot(snapshot_id, directory=WHO_ROOT):
    """Load only an explicitly identified production snapshot, never a fixture."""
    if not isinstance(snapshot_id, str) or not re.fullmatch(r"who-gtb\d{4}-incidence-[0-9a-f]{12}", snapshot_id):
        raise DataValidationError("An explicit production WHO snapshot identifier is required.")
    directory = Path(directory)
    if (directory / snapshot_id).is_dir():
        directory = directory / snapshot_id
    try:
        path = contained_file(directory, snapshot_id + "_manifest.json")
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError) as exc:
        raise DataValidationError(f"Missing/invalid manifest for {snapshot_id}.") from exc
    filename = manifest.get("dataFile", {}).get("filename", "")
    payload = file_bytes(contained_file(directory, filename), normalize_lf=True)
    return validate_who_bundle(snapshot_id, manifest, payload)


def validate_who_bundle(snapshot_id, manifest, payload):
    """Pure semantic validation of already acquired bytes, also used at import."""
    required = {"contractVersion", "snapshotId", "snapshotKind", "isCompleteDataset", "reportYear",
                "source", "sourceFiles", "dataFile", "fields", "coverage", "importerVersion",
                "importerCommit", "transformation", "analyticalHash", "terms", "citation", "volatile"}
    if not required <= manifest.keys():
        raise DataValidationError("Incomplete WHO manifest.")
    if manifest["snapshotId"] != snapshot_id or manifest["contractVersion"] != "who_incidence_snapshot_v2":
        raise DataValidationError("Snapshot identity/schema mismatch; explicit migration required.")
    if manifest["snapshotKind"] != "production" or manifest["isCompleteDataset"] is not True:
        raise DataValidationError("Test fixtures/incomplete imports cannot supply national profiles.")
    if manifest["importerVersion"] != "who_incidence_importer_v2" or not manifest["transformation"]:
        raise DataValidationError("Unsupported transformation version.")
    if manifest["analyticalHash"] != _manifest_hash(manifest):
        raise DataValidationError("Manifest analytical checksum mismatch.")
    if [field["name"] for field in manifest["fields"]] != list(COLUMNS):
        raise DataValidationError("Manifest column schema mismatch.")
    if type(manifest["reportYear"]) is not int:
        raise DataValidationError("Invalid report year.")
    if not re.fullmatch(r"[0-9a-f]{40}", manifest["importerCommit"]):
        raise DataValidationError("Missing importer commit.")
    for source in manifest["sourceFiles"]:
        require_hash(source["sha256"])
    data_info = manifest["dataFile"]
    payload = payload.replace(b"\r\n", b"\n")
    if digest(payload) != require_hash(data_info["sha256"]):
        raise DataValidationError("WHO data checksum mismatch.")
    if not snapshot_id.endswith(digest(payload)[:12]):
        raise DataValidationError("WHO content identifier mismatch.")
    columns, rows = csv_rows(payload)
    if columns != COLUMNS:
        raise DataValidationError("Missing/unexpected WHO columns.")
    if len(rows) != data_info["rows"]:
        raise DataValidationError("WHO row-count mismatch.")
    seen, identities = set(), {}
    report_year = manifest["reportYear"]
    for row in rows:
        if not re.fullmatch(r"[A-Z]{3}", row["iso3"]) or not re.fullmatch(r"[A-Z]{2}", row["iso2"]):
            raise DataValidationError("Invalid stable country identifier.")
        if not re.fullmatch(r"\d{1,3}", row["iso_numeric"]):
            raise DataValidationError("Invalid numeric country identifier.")
        if not row["country"] or row["who_region"] not in {"AFR", "AMR", "EMR", "EUR", "SEA", "WPR"}:
            raise DataValidationError("Missing country name or invalid WHO region.")
        if not re.fullmatch(r"\d{4}", row["year"]):
            raise DataValidationError("Invalid observation year.")
        year = int(row["year"])
        if not 2000 <= year < report_year:
            raise DataValidationError("Observation year outside report coverage.")
        key = row["iso3"], year
        if key in seen:
            raise DataValidationError(f"Duplicate country/year: {key}")
        seen.add(key)
        identity = row["iso2"], row["iso_numeric"], row["who_region"]
        if row["iso3"] in identities and identities[row["iso3"]] != identity:
            raise DataValidationError("Inconsistent country identifiers.")
        identities[row["iso3"]] = identity
        if decimal(row["population"]) <= 0:
            raise DataValidationError("Population denominator must be positive.")
        rate = _triple(row, ("incidence_per_100k", "incidence_per_100k_lo", "incidence_per_100k_hi"))
        count = _triple(row, ("incident_cases", "incident_cases_lo", "incident_cases_hi"))
        _triple(row, ("hiv_positive_incidence_per_100k", "hiv_positive_incidence_per_100k_lo", "hiv_positive_incidence_per_100k_hi"))
        if row["incidence_status"] != ("estimated" if rate is not None else "not_estimated"):
            raise DataValidationError("Missingness status disagrees with published values.")
        if (rate is None) != (count is None):
            raise DataValidationError("Incidence rate/count missingness mismatch.")
    years = sorted({int(row["year"]) for row in rows})
    coverage = manifest["coverage"]
    if not years or coverage["yearRange"] != [years[0], years[-1]] or years[-1] != report_year - 1:
        raise DataValidationError("Report vintage/observation endpoint mismatch.")
    if len(identities) != coverage["countriesAndAreas"]:
        raise DataValidationError("Country coverage mismatch.")
    estimated_countries = {r["iso3"] for r in rows if r["incidence_status"] == "estimated"}
    if len(estimated_countries) != coverage["withEstimates"]:
        raise DataValidationError("Estimated-country coverage mismatch.")
    if snapshot_id == WHO_2025 and (len(rows), len(identities), years) != (5347, 217, list(range(2000, 2025))):
        raise DataValidationError("Pinned WHO 2025 coverage changed.")
    return {"snapshot_id": snapshot_id, "report_year": report_year, "final_observation_year": years[-1],
            "source_data_hash": digest(payload), "manifest": deepcopy(manifest),
            "loader_version": LOADER_VERSION, "rows": sorted(rows, key=lambda r: (r["iso3"], int(r["year"])))}


def freeze_who_bundle(manifest, payload, directory):
    """Offline annual-update boundary; no replacement or automatic migration."""
    snapshot_id = manifest["snapshotId"]
    validate_who_bundle(snapshot_id, manifest, payload)
    return publish_immutable(directory, snapshot_id, {
        snapshot_id + "_manifest.json": canonical(manifest) + b"\n",
        manifest["dataFile"]["filename"]: payload.replace(b"\r\n", b"\n"),
    })


def country_index(snapshot):
    result = []
    for iso3 in sorted({row["iso3"] for row in snapshot["rows"]}):
        rows = sorted((row for row in snapshot["rows"] if row["iso3"] == iso3), key=lambda row: int(row["year"]))
        available = [int(row["year"]) for row in rows if row["incidence_status"] == "estimated"]
        expected = set(range(2000, snapshot["final_observation_year"] + 1))
        warnings = []
        if not available:
            warnings.append("no_published_incidence_estimates")
        elif set(available) != expected:
            warnings.append("incomplete_incidence_history")
        result.append({"country_id": iso3, "display": {"name": rows[-1]["country"]},
                       "who_region": rows[0]["who_region"], "available_incidence_years": available,
                       "record_years": [int(row["year"]) for row in rows],
                       "completeness": "complete" if set(available) == expected else "unavailable" if not available else "partial",
                       "warnings": warnings, "snapshot_id": snapshot["snapshot_id"]})
    return result


def _provenance(snapshot):
    m = snapshot["manifest"]
    return {"snapshot_id": snapshot["snapshot_id"], "report_year": snapshot["report_year"],
            "final_observation_year": snapshot["final_observation_year"],
            "input_source_hash": snapshot["source_data_hash"], "manifest_analytical_hash": m["analyticalHash"],
            "source_files": deepcopy(m["sourceFiles"]), "source": deepcopy(m["source"]),
            "source_commit": None, "importer_commit": m["importerCommit"],
            "crosscheck_repository": (m.get("crossCheck") or {}).get("repository"),
            "crosscheck_commit": (m.get("crossCheck") or {}).get("commit"),
            "import_date": m["volatile"]["accessDate"], "transformation_version": m["importerVersion"],
            "transformation": m["transformation"], "terms": deepcopy(m["terms"]), "citation": m["citation"]}


def national_profile(snapshot, country_id):
    index = next((entry for entry in country_index(snapshot) if entry["country_id"] == country_id), None)
    if index is None:
        raise DataValidationError(f"Unknown country/area: {country_id}")
    rows = []
    for source in snapshot["rows"]:
        if source["iso3"] != country_id:
            continue
        rows.append({"year": int(source["year"]), "central": source["incidence_per_100k"] or None,
                     "lower": source["incidence_per_100k_lo"] or None, "upper": source["incidence_per_100k_hi"] or None,
                     "numerator": source["incident_cases"] or None, "denominator": source["population"],
                     "original_rate": source["incidence_per_100k"] or None, "unit": RATE_UNIT,
                     "quantity": "estimated_tb_incidence", "status": source["incidence_status"]})
    return seal({"schema_id": "country_tb.incidence_history.v1", "contract_version": PROFILE_VERSION,
                 "loader_version": LOADER_VERSION, "mode": "national_who_snapshot",
                 "geography_id": "who:" + country_id, "country_id": country_id,
                 "baseline_year": 2024, "report_year": snapshot["report_year"],
                 "final_observation_year": snapshot["final_observation_year"],
                 "rows": rows, "provenance": _provenance(snapshot), "display": index["display"],
                 "warnings": index["warnings"], "completeness": index["completeness"]})


def scaled_profile(national, geography_id, scale, *, rationale, provenance, review_status):
    if national["mode"] != "national_who_snapshot" or not geography_id.startswith(("iso3166-2:", "custom:")):
        raise DataValidationError("Scaling needs a national parent and explicit local geography.")
    if not rationale.strip() or review_status not in {"reviewed", "not_reviewed"} or not provenance:
        raise DataValidationError("Scaling requires rationale, provenance and review status.")
    years = [row["year"] for row in national["rows"]]
    factors = {year: numeric_text(decimal(scale)) for year in years} if not isinstance(scale, dict) else scale
    if set(factors) != set(years):
        raise DataValidationError("Scale factors must cover every parent year exactly.")
    rows = deepcopy(national["rows"])
    for row in rows:
        factor = decimal(factors[row["year"]])
        for key in ("central", "lower", "upper"):
            if row[key] is not None:
                row[key] = numeric_text(decimal(row[key]) * factor)
        # National count/population are retained only in the unmodified parent.
        row["numerator"], row["denominator"], row["original_rate"] = None, None, None
    result = {k: deepcopy(v) for k, v in national.items() if k != "profile_hash"}
    result.update(mode="subnational_scaled", geography_id=geography_id, rows=rows,
                  national_parent=deepcopy(national),
                  transformation={"rate_ratios": [{"year": y, "value": numeric_text(factors[y])} for y in years],
                                  "rationale": rationale, "provenance": provenance, "review_status": review_status})
    return seal(result)


def user_series_profile(mode, geography_id, rows, *, provenance, review_status, conversion=None):
    if mode not in {"subnational_user_series", "custom_population_series"} or not provenance:
        raise DataValidationError("Explicit user series mode and provenance are required.")
    output, years, warnings = [], set(), []
    for item in rows:
        row = deepcopy(item)
        year = row["year"]
        if type(year) is not int or not 1 <= year <= 2024 or year in years:
            raise DataValidationError("Invalid/duplicated user incidence year.")
        years.add(year)
        if row.get("unit") != RATE_UNIT or row.get("quantity") not in {"notifications", "estimated_tb_incidence"}:
            raise DataValidationError("Explicit quantity and incidence-rate unit required.")
        for field in ("original_rate", "numerator", "denominator", "lower", "upper"):
            row.setdefault(field, None)
            if row[field] is not None:
                row[field] = numeric_text(decimal(row[field]))
        if row["denominator"] is not None and decimal(row["denominator"]) == 0:
            raise DataValidationError("Denominator must be positive.")
        if row["original_rate"] is None:
            if row["numerator"] is None or row["denominator"] is None:
                raise DataValidationError("Need a rate or a numerator/denominator pair.")
            rate = decimal(row["numerator"]) / decimal(row["denominator"]) * 100000
        else:
            rate = decimal(row["original_rate"])
        if (row["lower"] is None) != (row["upper"] is None):
            raise DataValidationError("Both uncertainty bounds or neither required.")
        if row["lower"] is not None and not decimal(row["lower"]) <= rate <= decimal(row["upper"]):
            raise DataValidationError("Invalid uncertainty ordering.")
        if row["numerator"] is not None and row["denominator"] is not None and row["original_rate"] is not None:
            derived = decimal(row["numerator"]) / decimal(row["denominator"]) * 100000
            if abs(derived - rate) > max(Decimal("0.01"), rate * Decimal("0.01")):
                raise DataValidationError("Numerator/denominator/rate conflict beyond declared 1% or 0.01 tolerance.")
        row["central"] = numeric_text(rate)
        if row["quantity"] == "notifications":
            if conversion is None:
                warnings.append("notifications_not_estimated_incidence")
            else:
                if conversion.get("method") != "ascertainment_fraction" or conversion.get("review_status") != "reviewed" or not conversion.get("source"):
                    raise DataValidationError("A reviewed explicit notification conversion is required.")
                fraction = decimal(conversion["fraction"])
                if not 0 < fraction <= 1:
                    raise DataValidationError("Ascertainment fraction must be in (0,1].")
                row["central"] = numeric_text(rate / fraction)
                for field in ("lower", "upper"):
                    if row[field] is not None:
                        row[field] = numeric_text(decimal(row[field]) / fraction)
                row["original_quantity"] = "notifications"
                row["quantity"] = "estimated_tb_incidence"
        output.append(row)
    if not years or review_status not in {"reviewed", "not_reviewed"}:
        raise DataValidationError("Empty series or invalid review status.")
    gaps = sorted(set(range(min(years), max(years) + 1)) - years)
    if gaps:
        warnings.append("missing_incidence_years")
    return seal({"schema_id": "country_tb.incidence_history.v1", "contract_version": PROFILE_VERSION,
                 "loader_version": LOADER_VERSION, "mode": mode, "geography_id": geography_id,
                 "baseline_year": 2024, "report_year": None, "final_observation_year": max(years),
                 "rows": sorted(output, key=lambda row: row["year"]), "provenance": deepcopy(provenance),
                 "review_status": review_status, "conversion": conversion, "missing_years": gaps,
                 "rate_reconciliation_tolerance": {"relative": "0.01", "absolute_per_100k": "0.01"},
                 "warnings": sorted(set(warnings)), "completeness": "partial" if gaps else "complete"})


def compare_snapshots(old, new):
    def lookup(snapshot):
        return {(row["iso3"], int(row["year"])): row for row in snapshot["rows"]}
    before, after = lookup(old), lookup(new)
    changed = []
    for key in sorted(before.keys() & after.keys()):
        a, b = before[key], after[key]
        fields = {}
        for field in ("population", "incidence_per_100k", "incidence_per_100k_lo", "incidence_per_100k_hi",
                      "incident_cases", "incident_cases_lo", "incident_cases_hi",
                      "hiv_positive_incidence_per_100k", "hiv_positive_incidence_per_100k_lo", "hiv_positive_incidence_per_100k_hi"):
            if a[field] != b[field]:
                fields[field] = {"old": a[field] or None, "new": b[field] or None,
                                 "difference": numeric_text(decimal(b[field]) - decimal(a[field])) if a[field] and b[field] else None}
        if fields:
            changed.append({"country_id": key[0], "year": key[1], "changes": fields})
    countries_old, countries_new = {k[0] for k in before}, {k[0] for k in after}
    return {"old_snapshot": old["snapshot_id"], "new_snapshot": new["snapshot_id"],
            "countries_added": sorted(countries_new - countries_old), "countries_removed": sorted(countries_old - countries_new),
            "records_added": [list(k) for k in sorted(after.keys() - before.keys())],
            "records_removed": [list(k) for k in sorted(before.keys() - after.keys())], "historical_revisions": changed}
