# WHO incidence data: maintainer update guide

Status: **implemented and validated**. The application never downloads data;
maintainers run these steps by hand.

## 1. Obtain the source files

Download both files from WHO and keep them outside the repository:

* Estimates: `https://extranet.who.int/tme/generateCSV.asp?ds=estimates`
  (saved as `TB_burden_countries_YYYY-MM-DD.csv`)
* Dictionary: `https://extranet.who.int/tme/generateCSV.asp?ds=dictionary`
  (saved as `TB_data_dictionary_YYYY-MM-DD.csv`)

Record the download date. Check the WHO TB data page to confirm which Global
Tuberculosis Report round the files belong to. In a round published in year R,
the latest estimate year is R - 1.

## 2. Optional cross-check export

From a pinned local checkout of the Global TB Report repository for the same
round (read-only):

```bash
python scripts/export_gtbreport_crosscheck.py ../gtbreport2025 gtbreport_est_export.csv
git -C ../gtbreport2025 rev-parse HEAD      # record this commit
```

This needs the optional `pyreadr` package. Do not commit the export.

## 3. Build the snapshot

```bash
python scripts/import_who_incidence.py \
  --estimates TB_burden_countries_2026-09-24.csv \
  --dictionary TB_data_dictionary_2026-09-24.csv \
  --report-year 2025 --access-date 2026-09-24 --kind production \
  --crosscheck-export gtbreport_est_export.csv \
  --crosscheck-commit 666088cac1e20dd1e9e52016ea857a9c48e0ba8d
```

Optionally add `--expected-sha256 <hash>` to confirm that the estimates file is the
one expected. The importer writes to `data/who_incidence/snapshots/`:

* `<snapshotId>.csv`: the data, byte-identical for identical inputs;
* `<snapshotId>_manifest.json`: the contract, provenance, terms, coverage and
  checksums;
* `<snapshotId>_validation.json` and `_validation.md`: the validation findings.

The snapshot ID, `who-gtb<round>-incidence-<first 12 hex of data SHA-256>`, is
derived from the content, so the same data always get the same ID.

## 4. Review the validation summary

* **Fatal** findings stop the import; nothing is written.
* **Warnings** need review, for example values present in the analysis output but
  withheld from the public dataset.
* **Expected missingness** means an area has population data but no incidence
  estimates.
* **Incomplete series** means an area covers only part of the year range.

## 5. Replace the old snapshot

Keep exactly one snapshot per directory. The loader refuses to run if there are
several, so that report rounds are never mixed. Delete the previous round's four
files, add the new ones, then run:

```bash
python -m pytest -q tests/test_general_who_incidence.py tests/test_general_trend_estimation.py
```

Saved profiles keep their original snapshot ID, so analyses can still be traced
to the round they used.

## 6. Test fixture

The test fixture in `data/who_incidence/fixture/` is a small subset used by the
fast tests. Regenerate it only when the contract changes:

```bash
python scripts/import_who_incidence.py ... --kind test_fixture --subset AUS BRA GBR IDN PHL ZAF PRK ANT
```

## Troubleshooting

| Message | Action |
| --- | --- |
| "No WHO incidence snapshot is installed" | Run the importer (step 3). |
| "Checksum mismatch" | The data file was edited or corrupted. Re-run the importer. |
| "More than one snapshot is installed" | Remove the older round (step 5). |
| `schema_changed` | WHO renamed or removed a column. Update `engine/who_incidence/contract.py` and bump the contract version. |
| `mixed_vintage` or `unexpected_year_range` | The report year is wrong, or the file is from a different round. |
