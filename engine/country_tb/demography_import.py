"""Maintainer-only OFFLINE demographic import boundary.

Accepts already downloaded files. Never downloads or appends data.
"""
from __future__ import annotations

import csv
from decimal import Decimal
import gzip
import io
from pathlib import Path

from engine.country_tb.data_common import canonical, csv_bytes, csv_rows, decimal, digest, numeric_text, publish_immutable
from engine.country_tb.demography import CSV_COLUMNS, normalize_records

TRANSFORMATION_VERSION = "country_tb_demography_import_1.0.0"


def _manifest(source_id, vintage, category, source_name, source_bytes, url, records, *,
              projection_boundary, transformations, terms, citation, import_date):
    columns = CSV_COLUMNS
    records = sorted(records, key=lambda r: (r["geography_id"], int(r["year"]), r["projection_scenario"], r["sex"], int(r["age_lower"])))
    # Prove uniqueness and overlap before any publication.
    check = []
    for record in records:
        row = dict(record, source_id=source_id, source_vintage=vintage)
        row["age_upper_exclusive"] = None if row["age_upper_exclusive"] == "" else int(row["age_upper_exclusive"])
        check.append(row)
    validated, duplicates = normalize_records(check)
    if duplicates:
        raise ValueError("Importer produced duplicates.")
    plain = csv_bytes(records, columns)
    compressed = gzip.compress(plain, mtime=0)
    snapshot_id = source_id.replace("_", "-") + "-" + digest(compressed)[:12]
    manifest = {
        "snapshot_id": snapshot_id, "schema_version": "country_tb.demography_snapshot.v1",
        "transformation_version": TRANSFORMATION_VERSION, "source_id": source_id,
        "source_vintage": vintage, "source_category": category,
        "source_url": url, "source_commit": None, "import_date": import_date,
        "source_files": [{"filename": source_name, "sha256": digest(source_bytes), "url": url}],
        "data_file": {"filename": "population.csv.gz", "sha256": digest(compressed), "rows": len(records)},
        "review_status": "reviewed", "review_scope": "technical source, unit, age and sex mapping; no epidemiological calibration",
        "terms": terms, "citation": citation, "projection_boundary": projection_boundary,
        "transformations": transformations,
        "reconciliation": {"absolute_tolerance": "1", "relative_tolerance": "0.000000001"},
        "coverage": {"geographies": len({r["geography_id"] for r in records}),
                     "years": sorted({r["year"] for r in records}),
                     "sexes": sorted({r["sex"] for r in records}),
                     "source_estimate_projection_status": sorted({r["estimate_status"] for r in records})},
    }
    return manifest, compressed


def freeze_wpp(source_path, directory, *, import_date, years=(2024, 2025, 2026, 2027, 2028, 2029, 2030)):
    source_path = Path(source_path)
    source = source_path.read_bytes()
    stream = io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(source)), encoding="utf-8-sig", newline="")
    reader = csv.DictReader(stream)
    header = reader.fieldnames
    selected = [r for r in reader if r["ISO3_code"] and r["LocTypeName"] == "Country/Area" and int(r["Time"]) in years]
    if not selected:
        raise ValueError("No selected WPP country/area rows.")
    records = []
    for row in selected:
        if row["Variant"] != "Medium":
            raise ValueError("Expected an explicit Medium WPP source.")
        lo = int(row["AgeGrpStart"])
        upper = "" if row["AgeGrp"] == "100+" and row["AgeGrpSpan"] == "-1" else lo + int(row["AgeGrpSpan"])
        for sex, field in (("female", "PopFemale"), ("male", "PopMale"), ("both", "PopTotal")):
            records.append({"geography_id": "who:" + row["ISO3_code"], "geography_level": "national",
                "year": int(row["Time"]), "estimate_status": "projection", "projection_scenario": "medium",
                "sex": sex, "age_lower": lo, "age_upper_exclusive": upper,
                "population": numeric_text(decimal(row[field]) * 1000), "unit": "persons",
                "population_variant": "resident_population_1_july"})
    url = "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/CSV_FILES/" + source_path.name
    m, data = _manifest("un_wpp2024_medium", "WPP2024", "un_wpp", source_path.name, source, url, records,
        projection_boundary=2024, import_date=import_date,
        transformations=["Select countries/areas with ISO3 and years " + ",".join(map(str, years)) + " from the original WPP2024 Medium file.",
                         "Map published female/male/total columns to sex; retain both for checking only.",
                         "Multiply published thousands by 1000 exactly; preserve 100+ as an open age interval.",
                         "Retain WPP projection status from 2024; no interpolation, extrapolation or relabelling."],
        terms={"status": "CC_BY_3.0_IGO", "url": "https://creativecommons.org/licenses/by/3.0/igo/",
               "source": "UN WPP official downloads catalogue copyright/footer"},
        citation="United Nations, DESA, Population Division (2024). World Population Prospects 2024, Online Edition.")
    selected.sort(key=lambda r: (r["ISO3_code"], int(r["Time"]), int(r["AgeGrpStart"])))
    evidence = gzip.compress(csv_bytes(selected, header), mtime=0)
    m["evidence_files"] = [{"filename": "published_country_rows.csv.gz", "sha256": digest(evidence),
                            "description": "Unaltered source fields/values for every selected country/year/age row; filtered rows only."}]
    m["manifest_hash"] = digest(canonical(m))
    return publish_immutable(directory, m["snapshot_id"], {
        "manifest.json": canonical(m) + b"\n", "population.csv.gz": data, "published_country_rows.csv.gz": evidence})


def freeze_abs(source_path, directory, *, import_date, years=(2024,)):
    # Optional importer dependency only. Offline loaders do not import openpyxl.
    import openpyxl
    source_path = Path(source_path)
    source = source_path.read_bytes()
    workbook = openpyxl.load_workbook(io.BytesIO(source), read_only=True, data_only=True)
    records = []
    for sheet_name in ("Data1", "Data2"):
        rows = list(workbook[sheet_name].iter_rows(values_only=True))
        headers = rows[0]
        for values in rows[10:]:
            when = values[0]
            if not hasattr(when, "year") or when.year not in years:
                continue
            for column, label in enumerate(headers[1:], 1):
                if not label:
                    continue
                _, sex, age, *_ = [v.strip() for v in label.split(";")]
                sex = {"Female": "female", "Male": "male", "Persons": "both"}[sex]
                lower = int(age.split()[0])
                upper = "" if age == "100 and over" else lower + 1
                records.append({"geography_id": "who:AUS", "geography_level": "national",
                    "year": when.year, "estimate_status": "estimate", "projection_scenario": "historical",
                    "sex": sex, "age_lower": lower, "age_upper_exclusive": upper,
                    "population": numeric_text(decimal(values[column])), "unit": "persons",
                    "population_variant": "estimated_resident_population_30_june"})
    url = "https://www.abs.gov.au/statistics/people/population/national-state-and-territory-population/jun-2024/3101059.xlsx"
    m, data = _manifest("abs_erp_jun2024", "ABS_2024-12-12_release", "abs", source_path.name, source, url, records,
        projection_boundary=None, import_date=import_date,
        transformations=["Read Table 59 Data1 and Data2; exact 2024 source row, published persons by sex/single age.",
                         "Workbook annual timestamps say June 1; publication defines the stock at 30 June.",
                         "Preserve 100 and over as an open interval; retain Persons as audit-only published both-sex counts."],
        terms={"status": "CC_BY_4.0_with_ABS_exceptions", "url": "https://www.abs.gov.au/website-privacy-copyright-and-disclaimer",
               "no_endorsement": True},
        citation="Australian Bureau of Statistics (2024). National, state and territory population, June 2024, Table 59; released 12 December 2024.")
    m["reconciliation"] = {"absolute_tolerance": "0", "relative_tolerance": "0"}
    m["evidence_files"] = [{"filename": "published_source.xlsx", "sha256": digest(source)}]
    m["manifest_hash"] = digest(canonical(m))
    return publish_immutable(directory, m["snapshot_id"], {
        "manifest.json": canonical(m) + b"\n", "population.csv.gz": data, "published_source.xlsx": source})


def freeze_user_demography(source_path, metadata, directory, *, import_date):
    """Import authoritative entered national/subnational/custom counts, offline.

    Input CSV follows CSV_COLUMNS in persons; no inferred sex/age/year conversion.
    Review is supplied explicitly and is never promoted by the importer.
    """
    source_path = Path(source_path)
    source = source_path.read_bytes()
    header, rows = csv_rows(source)
    if header != CSV_COLUMNS:
        raise ValueError("User demographic CSV must follow the documented canonical columns.")
    for row in rows:
        row["year"], row["age_lower"] = int(row["year"]), int(row["age_lower"])
        row["age_upper_exclusive"] = int(row["age_upper_exclusive"]) if row["age_upper_exclusive"] else ""
    if metadata.get("review_status") not in {"reviewed", "not_reviewed"}:
        raise ValueError("User source needs an explicit review status.")
    m, data = _manifest(metadata["source_id"], metadata["source_vintage"], "user", source_path.name,
        source, metadata["source_url"], rows, projection_boundary=metadata.get("projection_boundary"), import_date=import_date,
        transformations=["User counts already in persons and half-open source age intervals; no age/sex split or year substitution."],
        terms=metadata["terms"], citation=metadata["citation"])
    m["review_status"] = metadata["review_status"]
    m["review_scope"] = metadata.get("review_scope", "User-declared source review")
    m["evidence_files"] = [{"filename": "published_user_source.csv", "sha256": digest(source)}]
    # Unlike the original built-in datasets, user sources can share identical
    # counts across different vintages; identity must include source provenance.
    identity = {k: v for k, v in m.items() if k not in {"snapshot_id", "import_date", "manifest_hash"}}
    m["snapshot_id"] = m["source_id"].replace("_", "-") + "-" + digest(canonical(identity))[:12]
    m["manifest_hash"] = digest(canonical(m))
    return publish_immutable(directory, m["snapshot_id"], {
        "manifest.json": canonical(m) + b"\n", "population.csv.gz": data,
        "published_user_source.csv": source})
