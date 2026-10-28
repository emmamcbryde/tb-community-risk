# Independent frozen-data recheck, Milestone 2B

Rechecked 2026-10-28 from a clean branch at
a4fceee5af7c2871943061ad9635084fd18cec89. No frozen data or 2A loader was edited.

## Independent method

The [recheck script](../../scripts/recheck_country_tb_frozen.py) imports no
production loader. It reads bytes with stdlib CSV/Decimal and the ABS workbook
with openpyxl in read-only mode. It independently computes file and manifest
hashes, full-field duplicate signatures and deterministic unique rows.
It reconstructs all 303 ABS rows and all 502,677 WPP rows from native source
fields, checking complete row identities with an order-free multiset
fingerprint plus explicit key uniqueness.

The committed WPP native evidence contains 167,559 selected original rows.
This session additionally reselected them from the preserved full
66,954,143-byte original archive and verified identical reconstructed records.
The full archive is a local audit artifact, not a network dependency.
Other checkouts can verify every transformed record against committed native
evidence without that archive. The receipt reports its availability explicitly.

See [machine receipt](independent_recheck_2b.json). A separate test layer compares
the actual loaders with these independent results, recomputes scientific
hashes using an independent canonicalizer, blocks network calls and restores
saved profiles while forbidding current-catalogue lookup.

    python scripts/recheck_country_tb_frozen.py --output abm/output/country_tb_m2b_audit/independent_recheck.json

## WHO 2025

Confirmed from bytes: 5,347 rows, 217 stable country/area IDs, 216 entities
with any incidence estimate, years 2000–2024. Report vintage is 2025, not
observation year 2025. Country/year keys are unique. All available rate,
count and HIV-rate uncertainty triples are finite, nonnegative and ordered.
There are 25 missing total-rate triples and 25 missing case-count triples
(all PRK), plus 895 missing HIV-rate triples. Missing values retain the
not_estimated label; they are not replaced by zero.

ANT (2000–2009), SCG (2000–2004), CUW and SXM (2010–2024) remain separate.
No name matching, territory merger or successor-state inference occurs.
PRK has demography but no usable incidence series. The selector has 217 entries.

Snapshot: who-gtb2025-incidence-6997e3b087d2.
Transformation: who_incidence_importer_v2; schema who_incidence_snapshot_v2.
Source file/dictionary hashes and URLs are in the receipt. Original public WHO
source bytes are not committed as a second complete raw archive: this recheck
verifies the frozen snapshot/manifest, not WHO's estimation machinery or R pipeline.

| WHO hash type | SHA-256 |
| --- | --- |
| Transformed CSV | 6997e3b087d26788b21a47d09968eca87b5b03c3cc50644ae4dcdcd860acfd07 |
| Manifest bytes | 22d21158b75a47f237f5bee8ecbdb75516e579ccd17fad4a83a9d3ccd80d368a |
| Manifest analytical content | 32f228fc4338a9a61d56fa4770799c2274a3b9696f814d2304b3e2a3ebdee765 |

Ordinary loading requires no network. The evidence review used published
guidance separately; no current download replaced frozen inputs. WHO notes
that live country profiles can change, supporting
[frozen identities](https://www.who.int/teams/global-tuberculosis-programme/data/).

## Exact duplicate decomposition

| Observed class | Groups | Involved rows | Excess rows |
| --- | ---: | ---: | ---: |
| Two identical rows/age: BOL, BRN, CIV, KOR, LBY, MMR, RUS, TWN, VEN, VNM | 1,010 | 2,020 | 1,010 |
| Three identical rows/VCT age | 101 | 303 | 202 |
| Total | 1,111 | 2,323 | 1,212 |

All repeated stored fields match exactly: location, identifier, year, age and
population. Maximum and typical discrepancy are zero. 21,917 rows yield
20,705 unique stored keys. The reconstructed CSV equals the committed
quarantine byte-for-byte; deduplication is idempotent and reversing/shuffling
the input yields identical output. A changed stored value under the same key
fails. No average or first/last selection occurs.

These are complete **stored legacy row** duplicate groups. ISO/year/age is
not a complete modern demographic key: upstream sex, vintage and generator
remain unrecorded. No new evidence establishes a sex/scenario/ingestion
explanation. All 1,212 excess rows remain exact stored duplicates; the original
cache remains untouched and excluded from the active baseline.

The modern complete key is source ID, source vintage, geography ID,
geography level, year, estimate/projection status, projection scenario, sex,
age lower/upper-exclusive bounds, unit and population variant. Conflicting
values/provenance under that key fail; absent dimensions are not invented.

## Native-source 2024 examples and sex reconciliation

| Selected source | Female | Male | Derived combined | Published both | Difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| Australia: ABS June 2024 **estimate** | 13,698,889 | 13,505,920 | 27,204,809 | 27,204,809 | 0 |
| Viet Nam: WPP 2024 Medium **projection** | 51,526,284 | 49,461,398 | 100,987,682 | 100,987,686 | −4 |
| PRK: WPP 2024 Medium **projection** | 13,390,111 | 13,108,724 | 26,498,835 | 26,498,834 | +1 |

Difference means derived minus published. Australia's inactive WPP alternative
derives 26,713,198 versus published 26,713,206 (−8); ABS exceeds that derived
alternative by 491,611. Selection retains ABS rather than averaging.
The [ABS release](https://www.abs.gov.au/statistics/people/population/national-state-and-territory-population/jun-2024)
defines ERP at 30 June 2024; the workbook's June-1 timestamp is a storage convention.
WPP's 1 July 2024 population belongs to its 2024–2100 projection series,
consistent with the [UN summary](https://population.un.org/wpp/assets/Files/WPP2024_Summary-of-Results.pdf)
separating estimates through 2023 from projections starting 2024.

Each canonical example has 202 records: 101 female and 101 male age cells.
Published both rows remain audit-only. No invented sex split, double counting
or silent correction occurs. Detailed counts reconcile exactly to sex-specific
display bands and combined totals. The 100+ cell remains open; it does not
become invented single-year ages. WPP cell differences in these examples are
at most one person, within the original absolute 1-person plus relative 1e-9
tolerance. Aggregate discrepancies are still reported and never adjusted away.

There is no baseline-year substitution or interpolation. WPP 2025–2030 remains
Medium projection; another variant requires a distinct source identity.
Installed coverage ends 2030. Source ranking and original-snapshot restoration
remain unchanged.

## Provenance identity

| Snapshot | Source ID / vintage | Transformed SHA-256 | Manifest analytical SHA-256 |
| --- | --- | --- | --- |
| abs-erp-jun2024-07ca11bcb092 | abs_erp_jun2024 / ABS_2024-12-12_release | 07ca11bcb092f6bd62c1fc16d3229164680bc68033eac9c5c64404d89f31118a | e8682e4dddd49d27ab6a39cedd2a1a30bea6793c6ad66d6695db55de80d6f0ef |
| un-wpp2024-medium-8028a7a3e163 | un_wpp2024_medium / WPP2024 | 8028a7a3e163eaa6146554a9167a2b6e1a3a87c95d7ad8ed006f4421c36c7ddd | a988c3f6428ea4dcdfdd001e88082ad0c963a8bc2e96c2bbf5186d5bbde93046 |

Original ABS SHA-256: 9d436fca6a0261732f7c6d8f51c834d36707c279f3878d9a815dc4014c8d8086.
Original WPP SHA-256: 31804a296b663716236cd26c46415271cc9386fb3b2be56aff6467ad32283dc8.
Both retain country_tb_demography_import_1.0.0 and
country_tb.demography_snapshot.v1. Loaded profiles remain revision
1.1.0 / country_tb_data_loader_1.0.0.

| Baseline scientific-content hash | SHA-256 |
| --- | --- |
| AUS | 11c7e058cbe54db4af1d302c0a95044807bacdabf36a2ea5003449b7c8fc8a29 |
| VNM | f08621c9e3574eeddc01e83629578962bcd7fa4a60bec4eaf7b5d2127ba2f08a |
| PRK | 78fd3c8c8ce0c98faf2842fefcd0d6e2d4c5602312d197def6653d0f0d005954 |

Manifest byte hashes differ from analytical hashes by design: analytical
hashes exclude the hash itself and documented receipt fields. All are reported
separately. No discrepancy from 2A counts, hashes or classifications was found.
The small WPP differences above are reproduced, not suppressed.
