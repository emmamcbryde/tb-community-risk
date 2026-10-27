"""Forensic classification of the preserved legacy demographic CSV."""
from __future__ import annotations

from collections import Counter, defaultdict
from decimal import Decimal
import json

from engine.country_tb.data_common import DataValidationError, csv_rows, csv_bytes, decimal, digest, numeric_text

CATEGORIES = {
    1: "exact_duplicate_rows", 2: "sex_omitted_from_key", 3: "both_and_sex_specific",
    4: "different_vintages", 5: "historical_vs_projection", 6: "different_scenarios",
    7: "overlapping_age_bands", 8: "repeated_cache_ingestion",
    9: "conflicting_complete_key", 10: "other_explicit_category",
}
AUDIT_COLUMNS = ("iso_code", "year", "age", "category", "category_name", "rows",
                 "excess_rows", "max_population_difference", "relative_difference",
                 "source_line_numbers", "row_sha256", "resolution")


def audit_legacy(payload):
    columns, rows = csv_rows(payload)
    if columns != ("location", "iso_code", "year", "age", "population"):
        raise DataValidationError("Unexpected legacy source schema; re-audit required.")
    grouped = defaultdict(list)
    for line, row in enumerate(rows, 2):
        decimal(row["population"])
        grouped[(row["iso_code"], row["year"], row["age"])].append((line, row))
    audit, canonical_rows = [], []
    for key, values in sorted(grouped.items(), key=lambda item: (item[0][0], int(item[0][1]), int(item[0][2]))):
        signatures = {tuple(row[column] for column in columns) for _, row in values}
        if len(values) > 1:
            counts = [decimal(row["population"]) for _, row in values]
            discrepancy = max(counts) - min(counts)
            category = 1 if len(signatures) == 1 else 9
            signature_hash = digest(json.dumps(sorted(signatures), separators=(",", ":"), ensure_ascii=False).encode())
            audit.append({
                "iso_code": key[0], "year": key[1], "age": key[2], "category": category,
                "category_name": CATEGORIES[category], "rows": len(values), "excess_rows": len(values) - 1,
                "max_population_difference": numeric_text(discrepancy),
                "relative_difference": numeric_text(discrepancy / max(counts)) if max(counts) else "0",
                "source_line_numbers": ";".join(str(line) for line, _ in values), "row_sha256": signature_hash,
                "resolution": "retain_one_identical_complete_legacy_row" if category == 1 else "fail_unresolved_conflict",
            })
        if len(signatures) != 1:
            raise DataValidationError(f"Legacy key conflict requiring manual review: {key}")
        # Reconstruct from its unique VALUE signature, never first/last file row.
        canonical_rows.append(dict(zip(columns, next(iter(signatures)))))
    categories = []
    for category, name in CATEGORIES.items():
        entries = [entry for entry in audit if entry["category"] == category]
        categories.append({"category": category, "name": name, "groups": len(entries),
                           "rows": sum(e["rows"] for e in entries), "excess_rows": sum(e["excess_rows"] for e in entries),
                           "countries": sorted({e["iso_code"] for e in entries}),
                           "years": sorted({int(e["year"]) for e in entries}),
                           "ages": sorted({int(e["age"]) for e in entries}),
                           "maximum_difference": max((decimal(e["max_population_difference"]) for e in entries), default=Decimal(0)).to_eng_string(),
                           "typical_difference_median": "0" if category == 1 else None,
                           "examples": entries[:2]})
    summary = {"schema_id": "country_tb.legacy_duplicate_audit.v1", "source_sha256": digest(payload),
               "source_rows": len(rows), "unique_legacy_keys": len(canonical_rows),
               "duplicate_groups": len(audit), "involved_rows": sum(e["rows"] for e in audit),
               "excess_rows": sum(e["excess_rows"] for e in audit), "categories": categories,
               "per_country_excess": dict(sorted(Counter({code: sum(e["excess_rows"] for e in audit if e["iso_code"] == code)
                                                         for code in {e["iso_code"] for e in audit}}).items())),
               "unresolved_numerical_conflicts": 0,
               "provenance_warning": "Original generator and exact external vintage unrecorded; quarantine only; sex unspecified."}
    return audit, summary, csv_bytes(canonical_rows, columns)
