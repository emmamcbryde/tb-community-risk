"""Frozen demographic loading, complete-key validation and source resolution."""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from decimal import Decimal
import json
from pathlib import Path
import re

from engine.country_tb.data_common import (
    DataValidationError, LOADER_VERSION, PROFILE_VERSION, REPO_ROOT, canonical,
    contained_file, csv_rows, decimal, digest, numeric_text, require_hash, seal,
)

DEMOGRAPHY_ROOT = REPO_ROOT / "data/demography/snapshots"
DEMOGRAPHIC_KEY = (
    "source_id", "source_vintage", "geography_id", "geography_level", "year",
    "estimate_status", "projection_scenario", "sex", "age_lower",
    "age_upper_exclusive", "unit", "population_variant",
)
CSV_COLUMNS = (
    "geography_id", "geography_level", "year", "estimate_status", "projection_scenario",
    "sex", "age_lower", "age_upper_exclusive", "population", "unit", "population_variant",
)
DISPLAY_BANDS = ((0, 5), (5, 15), (15, 45), (45, 65), (65, None))


def complete_key(row):
    try:
        return tuple(row[name] for name in DEMOGRAPHIC_KEY)
    except KeyError as exc:
        raise DataValidationError(f"Missing demographic key dimension: {exc}") from exc


def normalize_records(records):
    """Deduplicate identical complete records only; any same-key conflict fails."""
    unique, duplicates = {}, 0
    for source in records:
        row = deepcopy(source)
        if row.get("sex") not in {"female", "male", "both", "unknown"}:
            raise DataValidationError("Sex must be female, male, both or explicitly unknown.")
        if row.get("geography_level") not in {"national", "subnational", "custom"}:
            raise DataValidationError("Invalid geography level.")
        if row.get("estimate_status") not in {"estimate", "projection", "unknown"}:
            raise DataValidationError("Invalid estimate/projection status.")
        if row.get("unit") != "persons":
            raise DataValidationError("Canonical population counts must use persons.")
        if not row.get("projection_scenario") or not row.get("source_id") or not row.get("source_vintage"):
            raise DataValidationError("Missing source/vintage/scenario dimension.")
        if type(row.get("year")) is not int or row["year"] < 1:
            raise DataValidationError("Reference year must be an integer.")
        lo, hi = row.get("age_lower"), row.get("age_upper_exclusive")
        if type(lo) is not int or lo < 0 or (hi is not None and (type(hi) is not int or hi <= lo)):
            raise DataValidationError("Invalid half-open age interval.")
        if row["estimate_status"] == "estimate" and row["projection_scenario"] != "historical":
            raise DataValidationError("Historical estimates cannot be labelled a projection scenario.")
        if row["estimate_status"] == "projection" and row["projection_scenario"] in {"historical", "unknown"}:
            raise DataValidationError("Projection scenario must be explicit.")
        row["population"] = numeric_text(decimal(row["population"]))
        key = complete_key(row)
        if key in unique:
            # Display labels have no bearing on scientific record equality.
            a = {k: v for k, v in row.items() if k != "display"}
            b = {k: v for k, v in unique[key].items() if k != "display"}
            if a != b:
                raise DataValidationError(f"Conflicting demographic values/provenance under complete key: {key}")
            duplicates += 1
        else:
            unique[key] = row
    result = sorted(unique.values(), key=lambda row: canonical({name: row[name] for name in DEMOGRAPHIC_KEY}))
    # Overlap is checked within each complete series, preserving variants/vintages.
    groups = defaultdict(list)
    for row in result:
        key = tuple(row[k] for k in DEMOGRAPHIC_KEY if k not in {"age_lower", "age_upper_exclusive"})
        groups[key].append(row)
    for rows in groups.values():
        ordered = sorted(rows, key=lambda row: row["age_lower"])
        for previous, current in zip(ordered, ordered[1:]):
            if previous["age_upper_exclusive"] is None or current["age_lower"] < previous["age_upper_exclusive"]:
                raise DataValidationError("Overlapping incompatible demographic age intervals.")
    return result, duplicates


def age_coverage(rows):
    warnings, gaps = [], []
    for sex in sorted({row["sex"] for row in rows}):
        ordered = sorted((r for r in rows if r["sex"] == sex), key=lambda r: r["age_lower"])
        expected = 0
        for row in ordered:
            if row["age_lower"] != expected:
                gaps.append({"sex": sex, "age_lower": expected, "age_upper_exclusive": row["age_lower"]})
            expected = row["age_upper_exclusive"]
        if expected is not None:
            gaps.append({"sex": sex, "age_lower": expected, "age_upper_exclusive": None})
    if gaps:
        warnings.append("incomplete_age_coverage")
    return {"warnings": warnings, "gaps": gaps}


def canonical_sexes(records, *, absolute_tolerance="1", relative_tolerance="0.000000001"):
    """Never includes published both-sex checks in male/female population sums."""
    rows, duplicate_count = normalize_records(records)
    if len({tuple(row[k] for k in DEMOGRAPHIC_KEY if k not in {"sex", "age_lower", "age_upper_exclusive"}) for row in rows}) != 1:
        raise DataValidationError("Select one geography/year/vintage/status/scenario/measure before canonical sex selection.")
    sexes = {row["sex"] for row in rows}
    if not rows:
        raise DataValidationError("No demographic records.")
    warnings, checks = [], []
    if "unknown" in sexes and len(sexes) > 1:
        raise DataValidationError("Unspecified-sex records cannot be combined with labelled sex totals.")
    if {"male", "female"} <= sexes:
        active = [row for row in rows if row["sex"] in {"male", "female"}]
        by_age = defaultdict(dict)
        for row in rows:
            by_age[(row["age_lower"], row["age_upper_exclusive"])][row["sex"]] = row
        for bounds, values in sorted(by_age.items()):
            if not {"male", "female"} <= values.keys():
                warnings.append("incomplete_sex_age_pairs")
                continue
            derived = decimal(values["male"]["population"]) + decimal(values["female"]["population"])
            if "both" in values:
                published = decimal(values["both"]["population"])
                difference = derived - published
                tolerance = decimal(absolute_tolerance) + decimal(relative_tolerance) * abs(published)
                if abs(difference) > tolerance:
                    raise DataValidationError(f"Male+female disagree with published both-sex total at {bounds}: {difference}")
                checks.append({"age_lower": bounds[0], "age_upper_exclusive": bounds[1],
                               "derived": numeric_text(derived), "published": numeric_text(published),
                               "difference": numeric_text(difference),
                               "relative_difference": numeric_text(difference / published) if published else None})
        sex_specific = True
    elif "both" in sexes:
        active = [row for row in rows if row["sex"] == "both"]
        warnings.append("sex_specific_modelling_unavailable")
        if sexes != {"both"}:
            warnings.append("incomplete_sex_records_preserved_but_inactive")
        sex_specific = False
    else:
        active = rows
        warnings.append("sex_unknown" if sexes == {"unknown"} else "incomplete_sex_coverage")
        sex_specific = False
    coverage = age_coverage(active)
    warnings += coverage["warnings"]
    return {"rows": sorted(active, key=lambda row: (row["sex"], row["age_lower"])),
            "raw_records": rows, "exact_duplicates_removed": duplicate_count,
            "sex_specific_available": sex_specific, "warnings": sorted(set(warnings)),
            "missing_age_intervals": coverage["gaps"], "published_both_checks": checks,
            "reconciliation_tolerance": {"absolute_persons": numeric_text(absolute_tolerance),
                                         "relative": numeric_text(relative_tolerance)}}


def aggregate_display(rows):
    """Exact decimal count aggregation; no disaggregation of coarse/open bands."""
    totals = defaultdict(Decimal)
    for row in rows:
        matches = [(lo, hi) for lo, hi in DISPLAY_BANDS if row["age_lower"] >= lo and
                   (hi is None or (row["age_upper_exclusive"] is not None and row["age_upper_exclusive"] <= hi))]
        if len(matches) != 1:
            raise DataValidationError("An age interval crosses a display boundary; explicit transformation required.")
        lo, hi = matches[0]
        totals[(row["sex"], lo, hi)] += decimal(row["population"])
    result = []
    for sex in sorted({row["sex"] for row in rows}):
        for lo, hi in DISPLAY_BANDS:
            result.append({"sex": sex, "age_lower": lo, "age_upper_exclusive": hi,
                           "population": numeric_text(totals[(sex, lo, hi)])})
    if sum((decimal(r["population"]) for r in result), Decimal(0)) != sum((decimal(r["population"]) for r in rows), Decimal(0)):
        raise DataValidationError("Display aggregation does not reconcile.")
    return result


def combined_rows(rows):
    active_sexes = {row["sex"] for row in rows}
    if "both" in active_sexes and len(active_sexes) != 1:
        raise DataValidationError("Both-sex records must not be counted with sex-specific rows.")
    if active_sexes not in ({"female", "male"}, {"both"}, {"unknown"}):
        raise DataValidationError("Cannot derive a combined population from incomplete sex coverage.")
    if active_sexes == {"female", "male"}:
        male = {(r["age_lower"], r["age_upper_exclusive"]) for r in rows if r["sex"] == "male"}
        female = {(r["age_lower"], r["age_upper_exclusive"]) for r in rows if r["sex"] == "female"}
        if male != female:
            raise DataValidationError("Sex-specific age supports must match before deriving combined rows.")
    totals = defaultdict(Decimal)
    for row in rows:
        totals[(row["age_lower"], row["age_upper_exclusive"])] += decimal(row["population"])
    return [{"sex": "both" if active_sexes != {"unknown"} else "unknown",
             "age_lower": bounds[0], "age_upper_exclusive": bounds[1], "population": numeric_text(total)}
            for bounds, total in sorted(totals.items())]


def load_demography_snapshot(snapshot_id, directory=DEMOGRAPHY_ROOT):
    if not isinstance(snapshot_id, str) or not re.fullmatch(r"[a-zA-Z0-9_.-]+", snapshot_id):
        raise DataValidationError("Explicit safe demographic snapshot ID required.")
    folder = Path(directory) / snapshot_id
    try:
        manifest = json.loads(contained_file(folder, "manifest.json").read_text(encoding="utf-8"))
    except (ValueError, OSError) as exc:
        raise DataValidationError("Missing/invalid demographic manifest.") from exc
    required = {"snapshot_id", "schema_version", "transformation_version", "source_id", "source_vintage",
                "source_category", "source_files", "data_file", "review_status", "terms", "citation",
                "source_url", "import_date", "projection_boundary", "transformations", "coverage", "manifest_hash"}
    if not required <= manifest.keys() or manifest["snapshot_id"] != snapshot_id:
        raise DataValidationError("Incomplete demographic manifest.")
    if manifest["schema_version"] != "country_tb.demography_snapshot.v1" or manifest["transformation_version"] != "country_tb_demography_import_1.0.0":
        raise DataValidationError("Unsupported demographic schema/transformation; explicit migration required.")
    if manifest["manifest_hash"] != digest(canonical({k: v for k, v in manifest.items() if k != "manifest_hash"})):
        raise DataValidationError("Demographic manifest checksum mismatch.")
    if not manifest["source_files"] or manifest["source_category"] not in {"user", "abs", "un_wpp", "fallback"}:
        raise DataValidationError("Missing source files or unsupported source hierarchy category.")
    if manifest["review_status"] not in {"reviewed", "not_reviewed", "rejected"}:
        raise DataValidationError("Invalid demographic review status.")
    for source in manifest["source_files"]:
        require_hash(source["sha256"])
    for evidence in manifest.get("evidence_files", []):
        if digest(contained_file(folder, evidence["filename"]).read_bytes()) != require_hash(evidence["sha256"]):
            raise DataValidationError("Demographic published-source evidence checksum mismatch.")
    info = manifest["data_file"]
    payload = contained_file(folder, info["filename"]).read_bytes()
    if digest(payload) != require_hash(info["sha256"]):
        raise DataValidationError("Demographic data checksum mismatch.")
    columns, data = csv_rows(payload, compressed=info["filename"].endswith(".gz"))
    if columns != CSV_COLUMNS or len(data) != info["rows"]:
        raise DataValidationError("Demographic CSV schema/row-count mismatch.")
    rows = []
    for record in data:
        row = dict(record)
        try:
            row["year"], row["age_lower"] = int(row["year"]), int(row["age_lower"])
            row["age_upper_exclusive"] = int(row["age_upper_exclusive"]) if row["age_upper_exclusive"] else None
        except ValueError as exc:
            raise DataValidationError("Invalid demographic year/age.") from exc
        row.update(source_id=manifest["source_id"], source_vintage=manifest["source_vintage"],
                   source_checksum=manifest["source_files"][0]["sha256"],
                   review_status=manifest["review_status"])
        boundary = manifest["projection_boundary"]
        if boundary is not None and (row["estimate_status"] == "projection") != (row["year"] >= boundary):
            raise DataValidationError("Historical/projection boundary mismatch.")
        rows.append(row)
    rows, duplicate_count = normalize_records(rows)
    if duplicate_count:
        raise DataValidationError("Frozen transformed snapshot contains duplicates.")
    if len({r["geography_id"] for r in rows}) != manifest["coverage"]["geographies"]:
        raise DataValidationError("Demographic country coverage mismatch.")
    if sorted({r["year"] for r in rows}) != manifest["coverage"]["years"]:
        raise DataValidationError("Demographic year coverage mismatch.")
    return {"snapshot_id": snapshot_id, "manifest": manifest, "rows": rows,
            "input_source_hash": manifest["source_files"][0]["sha256"], "transformed_file_hash": info["sha256"]}


def resolve_source(geography_id, candidates, *, year=2024, scenario="medium", preferred_snapshot_id=None, year_rule=None):
    """Pure deterministic ranking; equal eligible priorities require explicit choice."""
    country = geography_id.removeprefix("who:") if geography_id.startswith("who:") else "AUS" if geography_id.startswith("iso3166-2:AU-") else None
    ranks = {"user": 0, "abs": 1, "un_wpp": 2, "fallback": 3} if country == "AUS" else {"user": 0, "un_wpp": 1, "fallback": 2}
    considered, eligible = [], []
    for candidate in candidates:
        m = candidate["manifest"]
        rows = [r for r in candidate["rows"] if r["geography_id"] == geography_id and r["year"] == year
                and (r["estimate_status"] == "estimate" or r["projection_scenario"] == scenario)]
        source_year = year
        if not rows and year_rule is not None:
            if year_rule.get("method") != "nearest" or year_rule.get("review_status") != "reviewed" or not year_rule.get("rationale"):
                raise DataValidationError("Only an explicitly reviewed nearest-year rule is implemented; interpolation/extrapolation is not automatic.")
            possible = [r for r in candidate["rows"] if r["geography_id"] == geography_id and
                        (r["estimate_status"] == "estimate" or r["projection_scenario"] == scenario)]
            if possible:
                available_years = {r["year"] for r in possible}
                distances = {y: abs(year-y) for y in available_years}
                closest = [y for y in available_years if distances[y] == min(distances.values())]
                if len(closest) != 1:
                    raise DataValidationError("Equally close reference years require an explicit source year.")
                source_year = closest[0]
                if abs(source_year-year) <= year_rule.get("maximum_year_displacement", 0):
                    rows = [r for r in possible if r["year"] == source_year]
        reason = None
        if m["source_category"] not in ranks:
            reason = "not_in_approved_source_hierarchy"
        elif m["review_status"] != "reviewed":
            reason = "source_not_reviewed"
        elif not rows:
            reason = "requested_geography_year_scenario_unavailable"
        item = {"snapshot_id": candidate["snapshot_id"], "source_id": m["source_id"],
                "source_category": m["source_category"], "eligible": reason is None,
                "reason": reason, "rank": ranks.get(m["source_category"]), "source_year": source_year,
                "age_completeness": "partial" if rows and age_coverage(rows)["gaps"] else "complete" if rows else "unavailable"}
        considered.append(item)
        if reason is None:
            eligible.append((ranks[m["source_category"]], candidate["snapshot_id"], candidate, rows))
    if preferred_snapshot_id is not None:
        eligible = [item for item in eligible if item[1] == preferred_snapshot_id]
    if not eligible:
        raise DataValidationError("No reviewed source covers the exact geography/year/scenario; no cached fallback or relabelling.")
    eligible.sort(key=lambda item: (item[0], item[1]))
    if preferred_snapshot_id is None and len(eligible) > 1 and eligible[0][0] == eligible[1][0]:
        raise DataValidationError("Equally ranked frozen vintages require an explicit snapshot selection.")
    _, selected_id, selected, selected_rows = eligible[0]
    selected_year = selected_rows[0]["year"]
    warnings = []
    if selected["manifest"]["source_category"] == "fallback":
        warnings.append("documented_fallback_source_selected")
    if year == 2024 and any(row["estimate_status"] == "projection" for row in selected_rows):
        warnings.append("baseline_2024_is_source_projection")
    if selected_year != year:
        warnings.append("baseline_uses_explicit_nearest_year_proxy")
    for item in considered:
        if item["snapshot_id"] == selected_id:
            item["reason"] = "explicit_snapshot_restoration" if preferred_snapshot_id else "highest_priority_reviewed_exact_year_source"
        elif item["eligible"]:
            item["reason"] = "lower_priority_or_unselected_frozen_source"
    return selected, selected_rows, {"policy_version": "country_tb_source_hierarchy_1.0.0",
        "requested_year": year, "scenario": scenario, "selected_snapshot_id": selected_id,
        "candidates": sorted(considered, key=lambda x: x["snapshot_id"]), "warnings": warnings,
        "baseline_resolution": {"method": "exact_source_year" if selected_year == year else "reviewed_nearest_year_proxy",
                                "source_years": [selected_year], "reviewed_rule": year_rule,
                                "adjustment_persons": "0", "relative_count_adjustment": "0",
                                "year_displacement": selected_year-year}}


def demographic_profile(geography_id, candidates, *, year=2024, scenario="medium", preferred_snapshot_id=None, year_rule=None):
    selected, rows, resolution = resolve_source(geography_id, candidates, year=year, scenario=scenario,
                                               preferred_snapshot_id=preferred_snapshot_id, year_rule=year_rule)
    m = selected["manifest"]
    result = canonical_sexes(rows, **m.get("reconciliation", {}))
    warnings = sorted(set(result["warnings"] + resolution["warnings"]))
    try:
        combined = combined_rows(result["rows"])
    except DataValidationError:
        combined = []
        warnings.append("combined_population_unavailable")
    try:
        displayed = aggregate_display(result["rows"])
    except DataValidationError:
        displayed = []
        warnings.append("display_requires_explicit_age_disaggregation")
    totals = {sex: numeric_text(sum((decimal(r["population"]) for r in result["rows"] if r["sex"] == sex), Decimal(0)))
              for sex in sorted({r["sex"] for r in result["rows"]})}
    if sum((decimal(value) for value in totals.values()), Decimal(0)) <= 0:
        raise DataValidationError("Population must have a positive total.")
    inactive = [{"snapshot_id": c["snapshot_id"], "source_id": c["manifest"]["source_id"],
                 "years": sorted({r["year"] for r in c["rows"] if r["geography_id"] == geography_id}),
                 "scenarios": sorted({r["projection_scenario"] for r in c["rows"] if r["geography_id"] == geography_id})}
                for c in candidates]
    return seal({"schema_id": "country_tb.demography.v1", "contract_version": PROFILE_VERSION,
        "loader_version": LOADER_VERSION, "geography_id": geography_id, "geography_level": rows[0]["geography_level"],
        "baseline_year": year, "rows": result["rows"], "sex_totals": totals,
        "combined_population": numeric_text(sum((decimal(r["population"]) for r in result["rows"]), Decimal(0))) if combined else None,
        "combined_rows": combined, "display_bands": displayed,
        "published_both_checks": result["published_both_checks"],
        "sex_specific_available": result["sex_specific_available"],
        "missing_age_intervals": result["missing_age_intervals"], "warnings": warnings,
        "completeness": "complete" if not result["missing_age_intervals"] and "incomplete_sex_coverage" not in warnings else "partial",
        "source_resolution": resolution, "inactive_alternatives": sorted(inactive, key=lambda item: item["snapshot_id"]),
        "provenance": {"snapshot_id": selected["snapshot_id"], "source_id": m["source_id"],
            "source_vintage": m["source_vintage"], "input_source_hash": selected["input_source_hash"],
            "source_files": deepcopy(m["source_files"]) if "source_files" in m else [],
            "source_commit": m.get("source_commit"),
            "transformed_file_hash": selected["transformed_file_hash"], "manifest_hash": m["manifest_hash"],
            "schema_version": m["schema_version"], "transformation_version": m["transformation_version"],
            "transformations": m["transformations"], "review_status": m["review_status"],
            "source_url": m["source_url"], "import_date": m["import_date"], "terms": m["terms"], "citation": m["citation"],
            "projection_boundary": m["projection_boundary"], "reconciliation": result["reconciliation_tolerance"]}})
