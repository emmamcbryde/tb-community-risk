"""Read-only baseline data composition and developer preview exports."""
from __future__ import annotations

from decimal import Decimal
import json
from pathlib import Path

from engine.country_tb.data_common import (
    DataValidationError, LOADER_VERSION, PROFILE_VERSION, REPO_ROOT, canonical,
    csv_bytes, decimal, digest, numeric_text, publish_immutable, seal, profile_hash,
)
from engine.country_tb.incidence import WHO_2025, load_who_snapshot, national_profile
from engine.country_tb.demography import load_demography_snapshot, demographic_profile


def load_candidates(catalog_path=REPO_ROOT / "data/demography/catalog.json"):
    catalog_path = Path(catalog_path)
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    if catalog.get("schema_id") != "country_tb.demography_catalog.v1":
        raise DataValidationError("Unsupported demographic catalogue.")
    identifiers = catalog["snapshot_ids"]
    if len(identifiers) != len(set(identifiers)):
        raise DataValidationError("Duplicate demographic catalogue identifiers.")
    return [load_demography_snapshot(identifier, catalog_path.parent / "snapshots") for identifier in identifiers]


def baseline_data_profile(country_id, snapshot, candidates):
    if snapshot["snapshot_id"] != WHO_2025 or snapshot["report_year"] != 2025 or snapshot["final_observation_year"] != 2024:
        raise DataValidationError("This baseline explicitly binds WHO 2025, final observation 2024.")
    incidence = national_profile(snapshot, country_id)
    demographic = demographic_profile("who:" + country_id, candidates, year=2024)
    return seal({"schema_id": "country_tb.baseline_data_profile.v1", "contract_version": PROFILE_VERSION,
                 "loader_version": LOADER_VERSION, "country_id": country_id, "baseline_year": 2024,
                 "who_report_year": 2025, "incidence_final_observation_year": 2024,
                 "scope": "frozen_data_only", "incidence": incidence, "demography": demographic,
                 "warnings": sorted(set(incidence["warnings"] + demographic["warnings"])),
                 "display": incidence["display"]})


def preview_files(profile):
    country = profile["country_id"]
    incidence, demography = profile["incidence"], profile["demography"]
    files = {"profile.json": canonical(profile) + b"\n"}
    files["incidence.csv"] = csv_bytes(incidence["rows"], tuple(incidence["rows"][0]))
    people = [{k: row[k] for k in ("sex", "age_lower", "age_upper_exclusive", "population", "year", "estimate_status", "projection_scenario")}
              for row in demography["rows"]]
    files["population_by_sex_age.csv"] = csv_bytes(people, tuple(people[0]))
    for field in ("combined_rows", "display_bands", "published_both_checks"):
        rows = demography[field]
        files[field + ".csv"] = csv_bytes(rows, tuple(rows[0])) if rows else b"unavailable\n"
    latest = next((row for row in reversed(incidence["rows"]) if row["central"] is not None), None)
    lines = [
        "# " + country + " frozen baseline data preview", "",
        "Baseline demographic year: 2024. WHO report vintage: 2025. Final incidence observation: 2024.", "",
        "Scope: data review only. No infection pressure, TBI, individual population generation or projection.", "",
        "Selected demographic snapshot: " + demography["provenance"]["snapshot_id"],
        "Source vintage: " + demography["provenance"]["source_vintage"],
        "Combined population from canonical records: " + str(demography["combined_population"]),
        "Sex totals: " + json.dumps(demography["sex_totals"], sort_keys=True),
        "Sex-specific modelling available: " + str(demography["sex_specific_available"]),
        "Baseline resolution: " + json.dumps(demography["source_resolution"]["baseline_resolution"], sort_keys=True),
        "Incidence coverage: " + incidence["completeness"],
        "Latest published incidence: " + json.dumps(latest, sort_keys=True),
        "Warnings: " + ", ".join(profile["warnings"]), "",
        "Scientific profile SHA-256: " + profile["profile_hash"],
        "WHO input SHA-256: " + incidence["provenance"]["input_source_hash"],
        "Demography input SHA-256: " + demography["provenance"]["input_source_hash"],
        "Demography transformed SHA-256: " + demography["provenance"]["transformed_file_hash"], "",
        "See profile.json for source candidates, rejected alternatives, transformations, terms and all checks.",
    ]
    files["README.md"] = ("\n".join(lines) + "\n").encode("utf-8")
    return files


def write_preview(profile, output=REPO_ROOT / "abm/output/country_tb_m2a_preview"):
    output = Path(output).resolve()
    if not output.is_relative_to((REPO_ROOT / "abm/output").resolve()):
        raise DataValidationError("Developer preview output must stay under abm/output.")
    # Receipt identity includes presentation changes; scientific identity excludes them.
    identifier = profile["country_id"] + "-" + digest(canonical(profile))[:12]
    return publish_immutable(output, identifier, preview_files(profile))


def restore_baseline_data(payload, *, who_directory=None, demography_directory=None):
    """Restore an analysis against its original IDs, without reranking a catalogue."""
    result = json.loads(payload) if isinstance(payload, (str, bytes)) else json.loads(canonical(payload))
    if result.get("schema_id") != "country_tb.baseline_data_profile.v1" or result.get("contract_version") != PROFILE_VERSION:
        raise DataValidationError("Explicit baseline data-profile migration required.")
    if result.get("profile_hash") != profile_hash(result):
        raise DataValidationError("Saved baseline data-profile checksum mismatch.")
    if (result.get("baseline_year"), result.get("who_report_year"), result.get("incidence_final_observation_year"),
        result.get("scope")) != (2024, 2025, 2024, "frozen_data_only"):
        raise DataValidationError("Saved baseline data-profile scope/year mismatch.")
    for component in ("incidence", "demography"):
        if result[component].get("profile_hash") != profile_hash(result[component]):
            raise DataValidationError("Saved component profile checksum mismatch.")
    incidence = result["incidence"]
    args = {} if who_directory is None else {"directory": who_directory}
    snapshot = load_who_snapshot(incidence["provenance"]["snapshot_id"], **args)
    expected = national_profile(snapshot, result["country_id"])
    if expected["profile_hash"] != incidence["profile_hash"]:
        raise DataValidationError("Saved incidence does not resolve to its original frozen source.")
    d = result["demography"]
    args = {} if demography_directory is None else {"directory": demography_directory}
    source = load_demography_snapshot(d["provenance"]["snapshot_id"], **args)
    rebuilt = demographic_profile(d["geography_id"], [source], year=d["baseline_year"],
        scenario=d["source_resolution"]["scenario"],
        year_rule=d["source_resolution"]["baseline_resolution"].get("reviewed_rule"))
    for key in ("rows", "sex_totals", "combined_population", "provenance"):
        if rebuilt[key] != d[key]:
            raise DataValidationError(f"Saved demography no longer matches original source: {key}")
    return result
