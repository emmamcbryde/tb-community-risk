"""Maintainer command: build a versioned WHO incidence snapshot from local source files.

This command never downloads data. Obtain the source explicitly first, then run
one of:

    python scripts/import_who_incidence.py who-csv --source-file TB_burden_countries.csv \
        --report-year 2025 --snapshot-id who-gtb2025-yyyy-mm-dd --kind production

    python scripts/import_who_incidence.py gtbreport-est --gtbreport-root ../gtbreport2025 \
        --report-year 2025 --snapshot-id ... --countries AUS BRA ... --kind test_fixture

The gtbreport route needs the optional ``pyreadr`` package. The application and
its tests read only the resulting local snapshot.
"""

from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.who_incidence.adapters import (  # noqa: E402
    GtbReportEstimatesAdapter,
    WhoBurdenCsvAdapter,
    transform,
)
from engine.who_incidence.snapshot import BUNDLED_SNAPSHOT_DIR, write_snapshot  # noqa: E402


WHO_REPORT_URL = "https://www.who.int/teams/global-programme-on-tuberculosis-and-lung-health/tb-reports"
LICENCE = {
    "status": "to_be_confirmed",
    "statement": (
        "WHO states that data are provided in accordance with WHO's data policy and their use is "
        "subject to WHO's terms and conditions. The Global TB Report repository has no licence file; "
        "only WHO estimate values are extracted here, no repository code."
    ),
    "policyUrl": "https://www.who.int/about/policies/publishing/data-policy",
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", choices=["who-csv", "gtbreport-est"])
    parser.add_argument("--source-file", type=Path, help="Local WHO burden-estimates CSV (who-csv).")
    parser.add_argument("--gtbreport-root", type=Path, help="Local checkout of the Global TB Report repository.")
    parser.add_argument("--report-year", type=int, required=True)
    parser.add_argument("--snapshot-id", required=True)
    parser.add_argument("--kind", choices=["test_fixture", "production"], required=True)
    parser.add_argument("--countries", nargs="*", help="ISO3 codes to keep (default: all).")
    parser.add_argument("--rate-decimals", type=int, default=2)
    parser.add_argument("--out-dir", type=Path, default=BUNDLED_SNAPSHOT_DIR / "snapshots")
    parser.add_argument("--extraction-date", default=date.today().isoformat())
    args = parser.parse_args(argv)

    if args.source == "who-csv":
        if not args.source_file:
            parser.error("who-csv requires --source-file")
        adapter = WhoBurdenCsvAdapter(args.source_file, report_year=args.report_year)
        dataset = "WHO TB burden estimates (country level)"
    else:
        if not args.gtbreport_root:
            parser.error("gtbreport-est requires --gtbreport-root")
        adapter = GtbReportEstimatesAdapter(args.gtbreport_root)
        dataset = f"WHO Global Tuberculosis Report {args.report_year} country incidence estimates"

    records, report = transform(adapter, iso3_filter=args.countries, rate_decimals=args.rate_decimals)
    for issue in report.issues:
        print(f"[{issue.severity}] {issue.code}: {issue.message}")
    if not report.is_valid:
        print("Validation failed; no snapshot written.")
        return 1

    upstream = adapter.upstream()
    upstream["files"] = [{"path": item.path, "sha256": item.sha256, "role": item.role} for item in adapter.source_files()]
    stem = args.snapshot_id
    manifest = write_snapshot(
        records,
        out_dir=args.out_dir,
        data_filename=f"{stem}.csv",
        manifest_filename=f"{stem}_manifest.json",
        manifest_fields={
            "snapshotId": args.snapshot_id,
            "snapshotKind": args.kind,
            "isCompleteDataset": args.kind == "production" and not args.countries,
            "sourceDataset": dataset,
            "sourceReportYear": args.report_year,
            "sourceUrl": f"{WHO_REPORT_URL}/global-tuberculosis-report-{args.report_year}",
            "publicDownloadUrl": "https://extranet.who.int/tme/generateCSV.asp?ds=estimates",
            "seriesType": "retrospective series from a single report round",
            "adapter": adapter.source_id,
            "upstream": upstream,
            "extractionDate": args.extraction_date,
            "licence": LICENCE,
            "citation": f"World Health Organization. Global tuberculosis report {args.report_year}. Geneva: WHO; {args.report_year}.",
        },
    )
    print(f"Wrote snapshot {manifest['snapshotId']}: {manifest['dataFile']['rows']} rows, sha256 {manifest['dataFile']['sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
