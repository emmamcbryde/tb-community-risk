"""Maintainer command: build a WHO TB incidence snapshot from local source files.

This command never downloads data and is never run by the application. Obtain the
source files first (see docs/who_data_update_guide.md), then run for example:

    python scripts/import_who_incidence.py \
        --estimates TB_burden_countries_2026-09-24.csv \
        --dictionary TB_data_dictionary_2026-09-24.csv \
        --report-year 2025 --access-date 2026-09-24 --kind production \
        --crosscheck-export gtbreport2025_est_export.csv --crosscheck-commit 666088c...

The optional cross-check export is a CSV of iso3, year, inc, inc.lo, inc.hi,
inc.num, inc.lo.num, inc.hi.num and pop from inc_mort/analysis/est.rda of the
Global TB Report repository (scripts/export_gtbreport_crosscheck.py).
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import subprocess
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from engine.who_incidence.who_import import (  # noqa: E402
    GTB_REPORT_UPSTREAM,
    WHO_DICTIONARY_URL,
    WHO_ESTIMATES_URL,
    ImportFailed,
    build_manifest,
    check_dictionary,
    cross_check,
    load_crosscheck_export,
    parse_who_estimates,
    snapshot_bytes,
    write_snapshot_files,
)
from engine.who_incidence.snapshot import FIXTURE_DIR, PRODUCTION_DIR  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--estimates", type=Path, required=True, help="Local WHO TB burden-estimates CSV.")
    parser.add_argument("--dictionary", type=Path, required=True, help="Local WHO TB data dictionary CSV.")
    parser.add_argument("--report-year", type=int, required=True, help="Global TB Report round, e.g. 2025.")
    parser.add_argument("--access-date", required=True, help="Date the WHO files were downloaded (YYYY-MM-DD).")
    parser.add_argument("--kind", choices=["production", "test_fixture"], required=True)
    parser.add_argument("--subset", nargs="*", help="ISO3 codes to keep (test fixtures only).")
    parser.add_argument("--crosscheck-export", type=Path, help="CSV exported from est.rda for cross-checking.")
    parser.add_argument("--crosscheck-commit", help="Commit of the Global TB Report repository used for the export.")
    parser.add_argument("--expected-sha256", help="Fail unless the estimates file has this SHA-256.")
    parser.add_argument("--out-dir", type=Path, help="Output directory (default by kind).")
    args = parser.parse_args(argv)
    if args.kind == "production" and args.subset:
        parser.error("A production snapshot must not be subset.")

    estimates = args.estimates.read_bytes()
    estimates_sha = hashlib.sha256(estimates).hexdigest()
    if args.expected_sha256 and args.expected_sha256.lower() != estimates_sha:
        print(f"Checksum mismatch: expected {args.expected_sha256}, found {estimates_sha}. No snapshot written.")
        return 1
    rows, report = parse_who_estimates(estimates, report_year=args.report_year)
    definitions = check_dictionary(args.dictionary.read_bytes(), report)
    crosscheck = None
    if args.crosscheck_export:
        reference = load_crosscheck_export(args.crosscheck_export.read_bytes(), report)
        crosscheck = {
            "repository": GTB_REPORT_UPSTREAM,
            "commit": args.crosscheck_commit,
            "sourceObjects": ["inc_mort/analysis/est.rda", "data/gtb/other/cty.rda"],
            "exportFile": args.crosscheck_export.name,
            "exportSha256": hashlib.sha256(args.crosscheck_export.read_bytes()).hexdigest(),
            **cross_check(rows, reference, report),
        }
    if args.subset:
        wanted = {code.upper() for code in args.subset}
        rows = [row for row in rows if row.iso3 in wanted]
    source_files = [
        {"role": "estimates", "filename": args.estimates.name, "url": WHO_ESTIMATES_URL, "sha256": estimates_sha},
        {
            "role": "dictionary",
            "filename": args.dictionary.name,
            "url": WHO_DICTIONARY_URL,
            "sha256": hashlib.sha256(args.dictionary.read_bytes()).hexdigest(),
            "definitions": definitions,
        },
    ]
    data = snapshot_bytes(rows)
    manifest = build_manifest(
        rows,
        report,
        data_bytes=data,
        report_year=args.report_year,
        source_files=source_files,
        crosscheck=crosscheck,
        access_date=args.access_date,
        importer_commit=_git_head(),
        kind=args.kind,
        subset_iso3=args.subset,
    )
    out_dir = args.out_dir or (PRODUCTION_DIR if args.kind == "production" else FIXTURE_DIR)
    print(report.summary_markdown())
    try:
        paths = write_snapshot_files(out_dir, data, manifest, report)
    except ImportFailed as exc:
        print(str(exc))
        return 1
    coverage = manifest["coverage"]
    print(
        f"Wrote {manifest['snapshotId']}: {manifest['dataFile']['rows']} rows, {coverage['countriesAndAreas']} countries/areas "
        f"({coverage['withEstimates']} with estimates), years {coverage['yearRange']}, sha256 {manifest['dataFile']['sha256']}"
    )
    for label, path in paths.items():
        print(f"  {label}: {path}")
    return 0


def _git_head() -> str | None:
    try:
        return subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"], check=True, capture_output=True, text=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


if __name__ == "__main__":
    raise SystemExit(main())
