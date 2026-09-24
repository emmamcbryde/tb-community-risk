# WHO incidence snapshot schema (contract `who_incidence_snapshot_v2`)

Status: **implemented and validated**. Definitions: `engine/who_incidence/contract.py`.

## Data file (CSV, UTF-8, LF, sorted by `iso3`, `year`)

| Column | WHO source | Definition |
| --- | --- | --- |
| `iso3` | `iso3` | ISO 3166-1 alpha-3 code |
| `iso2` | `iso2` | ISO 3166-1 alpha-2 code |
| `iso_numeric` | `iso_numeric` | ISO numeric code as text (leading zeros kept) |
| `country` | `country` | Country or area name as published |
| `who_region` | `g_whoregion` | AFR, AMR, EMR, EUR, SEA or WPR |
| `year` | `year` | Calendar year |
| `population` | `e_pop_num` | Estimated total population (UN Population Division) |
| `incidence_per_100k`, `_lo`, `_hi` | `e_inc_100k`, `_lo`, `_hi` | Estimated incidence of TB disease (all forms) per 100,000, with bounds |
| `incident_cases`, `_lo`, `_hi` | `e_inc_num`, `_lo`, `_hi` | Estimated incident cases, with bounds |
| `hiv_positive_incidence_per_100k`, `_lo`, `_hi` | `e_inc_tbhiv_100k`, `_lo`, `_hi` | TB incidence among people living with HIV, per 100,000 total population |
| `incidence_status` | derived | `estimated` when all three rate values are published; `not_estimated` when none are |

Values are copied as published text, with no re-rounding and no float
reformatting. WHO publishes rates and counts to about two significant figures.
World Bank income group and the estimation-method field are **not** in the public
CSV, so they are not included.

## Manifest fields

`contractVersion`, `snapshotId`, `snapshotKind` (`production` or
`test_fixture`), `isCompleteDataset`, `subsetIso3`, `reportYear`,
`seriesDescription`, `source` (publisher, dataset, URLs), `sourceFiles` (role,
filename, URL, SHA-256, dictionary definitions), `crossCheck` (repository,
commit, source objects, export SHA-256, values compared, disagreements, values
withheld from the public dataset, tolerance), `importerVersion`,
`importerCommit`, `transformation`, `dataFile` (filename, SHA-256, rows, bytes),
`coverage` (countries and areas, number with estimates, year range, incomplete
series), `fields` (definitions), `citation`, `terms` (status, URL, summary,
no-endorsement statement, open question), `validationSummary`,
`knownExclusions`, `volatile` (`accessDate`), `analyticalHash`.

`analyticalHash` is the SHA-256 of the canonical manifest without the volatile
access date and the citation, which contains that date. Re-importing identical
source files gives a byte-identical data file and the same `analyticalHash`.

## Loader guarantees (`engine/who_incidence/snapshot.py`)

The loader checks that:

* the manifest structure and contract version are valid;
* the data SHA-256 matches, after tolerating CRLF line endings from a checkout;
* the header matches the contract, and the row count matches;
* every row passes the importer's validation rules again, and rows are in
  canonical form.

It prefers a production snapshot, falls back to the fixture, refuses more than
one snapshot per directory, and gives a maintainer-facing error when none is
installed.

## Installed production snapshot

| Item | Value |
| --- | --- |
| Snapshot ID | `who-gtb2025-incidence-6997e3b087d2` |
| Rows | 5,347 |
| Countries and areas | 217 (216 with estimates; PRK has no published estimates) |
| Years | 2000-2024 |
| Incomplete series | 8 areas (for example ANT, 2000-2009) |
| Cross-check | 37,279 values published in both sources agree within published rounding with `est.rda` at commit `666088c`; 150 PRK values exist in `est.rda` but are withheld from the public dataset |
| Data SHA-256 | `6997e3b087d26788b21a47d09968eca87b5b03c3cc50644ae4dcdcd860acfd07` |
