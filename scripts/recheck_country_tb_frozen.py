"""Independent, offline Milestone 2B check; deliberately imports no model loaders.

Reconstructs demographic records from published source fields, not M2A reports.
Writes nothing unless --output is explicitly supplied under abm/output.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WHO = "who-gtb2025-incidence-6997e3b087d2"
ABS = "abs-erp-jun2024-07ca11bcb092"
WPP = "un-wpp2024-medium-8028a7a3e163"


def sha(payload):
    return hashlib.sha256(payload).hexdigest()


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                      allow_nan=False).encode("utf-8")


def independent_scientific_hash(payload):
    """Reproduce documented hash policy without importing data_common."""
    set_fields = {"rows", "source_files", "warnings", "candidates",
                  "inactive_alternatives", "published_both_checks"}
    def walk(value, key=""):
        if isinstance(value, dict):
            return {k: walk(v, k) for k, v in value.items()
                    if k not in {"display", "profile_hash", "receipt"}}
        if isinstance(value, list):
            values = [walk(v) for v in value]
            return sorted(values, key=encoded) if key in set_fields else values
        return value
    return sha(encoded(walk(payload)))


def text(value):
    result = format(Decimal(str(value)), "f")
    return result.rstrip("0").rstrip(".") if "." in result else result


def read_csv(payload):
    return list(csv.DictReader(io.StringIO(payload.decode("utf-8-sig"))))


def legacy_check(payload):
    """All five stored fields are tested; absent upstream dimensions not invented."""
    rows = read_csv(payload)
    groups = defaultdict(list)
    for row in rows:
        groups[(row["iso_code"], row["year"], row["age"])].append(row)
    unique = []
    repeated = Counter()
    countries = Counter()
    for key, values in sorted(groups.items()):
        signatures = {tuple(sorted(row.items())) for row in values}
        if len(signatures) != 1:
            raise ValueError(f"Conflicting complete stored key: {key}")
        unique.append(dict(next(iter(signatures))))
        if len(values) > 1:
            repeated[len(values)] += 1
            countries[key[0]] += len(values) - 1
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=("location", "iso_code", "year", "age", "population"),
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(sorted(unique, key=lambda r: (r["iso_code"], int(r["year"]), int(r["age"]))))
    canonical_bytes = stream.getvalue().encode()
    return {
        "source_sha256": sha(payload), "source_rows": len(rows), "unique_keys": len(unique),
        "duplicate_groups": sum(repeated.values()), "groups_by_multiplicity": dict(sorted(repeated.items())),
        "excess_rows": sum((n - 1) * count for n, count in repeated.items()),
        "per_country_excess": dict(sorted(countries.items())),
        "all_repeated_stored_fields_identical": True, "maximum_difference": "0",
        "canonical_sha256": sha(canonical_bytes),
    }, canonical_bytes


def who_check(root):
    directory = root / "data/who_incidence/snapshots"
    path = directory / (WHO + "_manifest.json")
    manifest_bytes = path.read_bytes()
    manifest = json.loads(manifest_bytes)
    analytical = sha(encoded({k: v for k, v in manifest.items()
                              if k not in {"volatile", "analyticalHash", "citation"}}))
    assert analytical == manifest["analyticalHash"]
    payload = (directory / manifest["dataFile"]["filename"]).read_bytes()
    assert sha(payload) == manifest["dataFile"]["sha256"]
    rows = read_csv(payload)
    assert len({(r["iso3"], r["year"]) for r in rows}) == len(rows)
    triples = (("incidence_per_100k", "incidence_per_100k_lo", "incidence_per_100k_hi"),
               ("incident_cases", "incident_cases_lo", "incident_cases_hi"),
               ("hiv_positive_incidence_per_100k", "hiv_positive_incidence_per_100k_lo",
                "hiv_positive_incidence_per_100k_hi"))
    missing = Counter()
    for row in rows:
        assert Decimal(row["population"]).is_finite() and Decimal(row["population"]) > 0
        for central, lower, upper in triples:
            values = [row[k] for k in (lower, central, upper)]
            if all(v == "" for v in values):
                missing[central] += 1
                continue
            assert all(v != "" for v in values)
            lo, mid, hi = map(Decimal, values)
            assert all(v.is_finite() and v >= 0 for v in (lo, mid, hi))
            assert lo <= mid <= hi
        assert row["incidence_status"] in {"estimated", "not_estimated"}
        assert (row["incidence_status"] == "not_estimated") == (row["incidence_per_100k"] == "")
    return {
        "snapshot_id": WHO, "report_year": manifest["reportYear"],
        "years": sorted({int(r["year"]) for r in rows}), "rows": len(rows),
        "countries_areas": len({r["iso3"] for r in rows}),
        "with_estimates": len({r["iso3"] for r in rows if r["incidence_per_100k"] != ""}),
        "transformed_file_sha256": sha(payload), "manifest_file_sha256": sha(manifest_bytes),
        "manifest_analytical_hash": analytical, "transformation_version": manifest["importerVersion"],
        "source_files": manifest["sourceFiles"], "missing_triples": dict(missing),
        "PRK_all_missing": all(r["incidence_per_100k"] == "" for r in rows if r["iso3"] == "PRK"),
        "entity_examples": {code: {"name": next(r["country"] for r in rows if r["iso3"] == code),
                                 "years": sorted(int(r["year"]) for r in rows if r["iso3"] == code)}
                            for code in ("ANT", "SCG", "CUW", "SXM", "VNM", "PRK")},
    }


def fingerprint_records(rows):
    """Order-independent multiset fingerprint plus explicit complete-key uniqueness."""
    total = 0
    keys = set()
    examples = defaultdict(dict)
    count = 0
    for r in rows:
        key = (r["geography_id"], int(r["year"]), r["sex"], int(r["age_lower"]))
        assert key not in keys
        keys.add(key)
        row = (key, r["age_upper_exclusive"], text(r["population"]), r["estimate_status"],
               r["projection_scenario"], r["unit"], r["population_variant"], r["geography_level"])
        total = (total + int(sha(encoded(row)), 16)) % (1 << 256)
        count += 1
        if int(r["year"]) == 2024 and r["geography_id"] in {"who:AUS", "who:VNM", "who:PRK"}:
            examples[r["geography_id"]][(r["sex"], int(r["age_lower"]))] = Decimal(r["population"])
    return (count, f"{total:064x}"), examples


def native_wpp(path, *, full_archive=False):
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as stream:
        for r in csv.DictReader(stream):
            if full_archive and (not r["ISO3_code"] or r["LocTypeName"] != "Country/Area"
                                 or not 2024 <= int(r["Time"]) <= 2030):
                continue
            assert r["Variant"] == "Medium" and r["LocTypeName"] == "Country/Area"
            assert r["ISO3_code"] and 2024 <= int(r["Time"]) <= 2030
            for sex, field in (("female", "PopFemale"), ("male", "PopMale"), ("both", "PopTotal")):
                yield {"geography_id": "who:" + r["ISO3_code"], "geography_level": "national",
                       "year": r["Time"], "sex": sex, "age_lower": r["AgeGrpStart"],
                       "age_upper_exclusive": "" if r["AgeGrp"] == "100+" and r["AgeGrpSpan"] == "-1"
                       else str(int(r["AgeGrpStart"]) + int(r["AgeGrpSpan"])),
                       "population": text(Decimal(r[field]) * 1000), "unit": "persons",
                       "population_variant": "resident_population_1_july",
                       "estimate_status": "projection", "projection_scenario": "medium"}


def native_abs(path):
    import openpyxl  # Read-only source evidence; not used by ordinary loaders.
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        for sheet in ("Data1", "Data2"):
            values = list(workbook[sheet].values)
            for row in values[10:]:
                if not hasattr(row[0], "year") or row[0].year != 2024:
                    continue
                for index, heading in enumerate(values[0][1:], 1):
                    if not heading:
                        continue
                    _, sex, age, *_ = [part.strip() for part in heading.split(";")]
                    lo = int(age.split()[0])
                    yield {"geography_id": "who:AUS", "geography_level": "national", "year": "2024",
                           "sex": {"Female": "female", "Male": "male", "Persons": "both"}[sex],
                           "age_lower": str(lo), "age_upper_exclusive": "" if "over" in age else str(lo + 1),
                           "population": text(row[index]), "unit": "persons",
                           "population_variant": "estimated_resident_population_30_june",
                           "estimate_status": "estimate", "projection_scenario": "historical"}
    finally:
        workbook.close()


def demographic_check(root, identifier):
    directory = root / "data/demography/snapshots" / identifier
    manifest_bytes = (directory / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    actual_hash = sha(encoded({k: v for k, v in manifest.items() if k != "manifest_hash"}))
    assert actual_hash == manifest["manifest_hash"]
    for evidence in manifest["evidence_files"]:
        assert sha((directory / evidence["filename"]).read_bytes()) == evidence["sha256"]
    transformed_path = directory / "population.csv.gz"
    assert sha(transformed_path.read_bytes()) == manifest["data_file"]["sha256"]
    with gzip.open(transformed_path, "rt", encoding="utf-8-sig", newline="") as stream:
        observed, examples = fingerprint_records(csv.DictReader(stream))
    evidence_path = directory / manifest["evidence_files"][0]["filename"]
    native = native_abs(evidence_path) if identifier == ABS else native_wpp(evidence_path)
    reconstructed, native_examples = fingerprint_records(native)
    assert observed == reconstructed and examples == native_examples
    result = {
        "snapshot_id": identifier, "source_id": manifest["source_id"], "source_vintage": manifest["source_vintage"],
        "source_files": manifest["source_files"], "transformed_file_hash": manifest["data_file"]["sha256"],
        "manifest_hash": actual_hash, "manifest_file_sha256": sha(manifest_bytes),
        "transformation_version": manifest["transformation_version"], "coverage": manifest["coverage"],
        "independently_reconstructed_rows": observed[0], "native_transformed_multiset_hash": observed[1],
        "examples": {},
    }
    bands = ((0, 5), (5, 15), (15, 45), (45, 65), (65, 101))
    for geo, cells in examples.items():
        totals = {s: sum(v for (sex, _), v in cells.items() if sex == s) for s in ("female", "male", "both")}
        combined = totals["female"] + totals["male"]
        discrepancy = [cells["female", a] + cells["male", a] - cells["both", a] for a in range(101)]
        per_sex = {s: [text(sum(cells[s, a] for a in range(lo, hi))) for lo, hi in bands]
                   for s in ("female", "male")}
        assert all(sum(map(Decimal, values)) == totals[s] for s, values in per_sex.items())
        result["examples"][geo] = {
            "sex_totals": {k: text(v) for k, v in totals.items()}, "derived_combined": text(combined),
            "derived_minus_published_both": text(combined - totals["both"]),
            "max_cell_absolute_difference": text(max(map(abs, discrepancy))),
            "display_bands_by_sex": per_sex,
            "estimate_status": "estimate" if identifier == ABS else "projection",
            "scenario": "historical" if identifier == ABS else "medium",
            "reference_date": "2024-06-30" if identifier == ABS else "2024-07-01",
        }
    # Full archive is intentionally not necessary for ordinary loading/rechecking.
    archive = root / "abm/output/country_tb_m2a_audit/sources" / manifest["source_files"][0]["filename"]
    result["full_original_archive_present_locally"] = archive.is_file()
    if archive.is_file():
        result["full_original_archive_sha256"] = sha(archive.read_bytes())
        assert result["full_original_archive_sha256"] == manifest["source_files"][0]["sha256"]
        if identifier == WPP:
            original_records, original_examples = fingerprint_records(native_wpp(archive, full_archive=True))
            assert original_records == observed and original_examples == examples
            result["full_original_selected_records_match"] = True
    return result


def recheck(root=ROOT):
    source = (root / "data/population_age_latest.csv").read_bytes()
    legacy, unique = legacy_check(source)
    second, again = legacy_check(unique)
    assert again == unique and second["excess_rows"] == 0
    assert legacy_check(b"\n".join([source.splitlines()[0]] + list(reversed(source.splitlines()[1:]))) + b"\n")[1] == unique
    assert unique == (root / "data/demography/legacy_audit/legacy_unique_quarantine.csv").read_bytes()
    return {"schema_id": "country_tb.independent_frozen_recheck.v1", "method": "source-artifact reconstruction; no production loader imports",
            "who": who_check(root), "legacy": legacy,
            "demography": [demographic_check(root, identifier) for identifier in (ABS, WPP)]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    result = recheck()
    if args.output:
        path = args.output.resolve()
        if not path.is_relative_to((ROOT / "abm/output").resolve()):
            raise ValueError("Recheck output must stay under abm/output.")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(json.dumps(result, indent=2, ensure_ascii=False).encode() + b"\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
