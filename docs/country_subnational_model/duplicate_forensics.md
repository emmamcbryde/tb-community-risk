# Demographic duplicate forensic report

## Original evidence

The untouched source is data/population_age_latest.csv, 21,917 records. SHA-256: e51101fa3f5054b9aeb11999d38f6b39a3b0ab4625daa00c5c7fd52e77e535d4. Its five columns are location, iso_code, year, age and population. It has no sex, source revision, projection scenario or estimate-status column.

The original creation commit dca898f93c177d0a9b3848419666820e31b21dea (13 November 2025) already contained all 1,212 excess records, among 21,816 rows. Its file SHA-256 is 6033c39194a530ef11c9bf53ee8954ac0f6fa705d4d19781cc47d41ba2f9678f. The only later file change, commit 97ebd2b321f07a39a13e2d0c7037d91b3935868e, added 101 South Australia Aboriginal records under SAA. It added no duplicates.

Neither commit contains a generator. Repository searches found consumers and an older live World Bank/OWID fetch helper, but no code that writes this cache. UI comments attribute it to OWID/UN population data; they do not identify the originating export or UN revision. The live helper uses a different age-group source shape and cannot establish this file's lineage. The exact upstream extraction and any pre-Git repeated ingestion remain unproven.

There are 205 distinct identifiers, not 205 verified national ISO3 geographies: six OWID regional/world aggregates and custom SAA occur alongside country codes. Therefore even a unique country/year/age table would not be a validated national demographic dataset.

## Exact decomposition

The earlier 1,212 figure counts excess records, not duplicate groups. There are 1,111 groups containing 2,323 records:

| Observed category | Groups | Involved rows | Excess rows | Countries | Years and ages | Maximum / median numerical discrepancy | Resolution |
| --- | ---: | ---: | ---: | --- | --- | --- | --- |
| 1 Exact identical stored rows | 1,111 | 2,323 | 1,212 | Eleven below | 2021; every age 0–100 | 0 / 0 persons; relative 0 | Retain one unique full-row signature in quarantine |
| 2 Sex omitted from key | 0 established | 0 | 0 | None established | None | Not applicable | No reconstructed sex split |
| 3 Female/male/both totals together | 0 established | 0 | 0 | None established | None | Not applicable | No such dimension in the cache |
| 4 Different source vintages | 0 established | 0 | 0 | None established | None | Not applicable | External vintage remains unknown |
| 5 Historical versus projection | 0 established | 0 | 0 | None established | None | Not applicable | Status remains unknown |
| 6 Different scenarios | 0 established | 0 | 0 | None established | None | Not applicable | Scenario remains unknown |
| 7 Overlapping age bands | 0 | 0 | 0 | None | None | Not applicable | Source has single-age labels |
| 8 Repeated cache ingestion | 0 established | 0 | 0 | None established | None | Not applicable | Already present on creation; no writer to repair retrospectively |
| 9 Conflicting observed complete records | 0 | 0 | 0 | None | None | 0 conflicts | Any future conflict fails |
| 10 Other described duplicates | 0 | 0 | 0 | None | None | Not applicable | No remainder |

Ten countries have each age duplicated: BOL, BRN, CIV, KOR, LBY, MMR, RUS, TWN, VEN and VNM. Each contributes 101 excess rows: 10 × 101 = 1,010. VCT has each age triplicated: 101 × 2 = 202 excess rows. Total 1,010 + 202 = 1,212.

Examples: Bolivia age 0 occurs twice with population 237,991; Saint Vincent and the Grenadines age 0 occurs three times with population 1,440. The machine audit records original CSV line numbers, full-row signature hash, numerical discrepancy and resolution for all 1,111 groups. Its category summary lists all ten categories, affected countries, years, ages and representative records.

These classifications describe the available stored records. They do not prove that an unrecorded earlier transformation never omitted a dimension. That provenance uncertainty is why the repaired view remains quarantined.

## Complete canonical key

The new complete key is:

    source_id, source_vintage, geography_id, geography_level, year,
    estimate_status, projection_scenario, sex, age_lower,
    age_upper_exclusive, unit, population_variant

Source checksum, review and transformation provenance accompany each record or its manifest. Different provenance under an otherwise identical complete key is a conflict rather than permission to choose a file-order winner.

## Applied and implemented rules

For the legacy cache, deduplication was applied only after equality of every stored column was verified. A canonical row is reconstructed from the unique full-row value signature; no first/last-row choice or average is used. Output has 20,705 unique legacy keys. The original file and checksum remain unchanged. No undocumented field, sex split or 2024 year was invented.

For new snapshots, exact records can be deduplicated before publication, while different sexes, vintages, historical/projection status and scenarios remain separate complete keys. Published both-sex rows remain raw audit evidence; active male/female records exclude those totals. Conflicts and overlapping incompatible age intervals fail validation.

New imports validate complete-key uniqueness, stage a complete directory and atomically publish a new immutable snapshot. Repeating identical content is a no-op. Changed content cannot append to or overwrite an installed snapshot. This establishes an idempotent writer for the new data layer without modifying legacy consumers or pretending an absent old writer was recovered.

## Machine evidence

- data/demography/legacy_audit/duplicate_groups.csv
- data/demography/legacy_audit/duplicate_summary.json
- data/demography/legacy_audit/legacy_unique_quarantine.csv

Tests reconstruct the decomposition from the original bytes, compare the quarantined output and verify input-order independence and conflict rejection. No unresolved numerical conflicts remain in the observed cache; its external source revision, generator and sex semantics remain unresolved.
