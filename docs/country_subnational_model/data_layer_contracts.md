# Executable frozen data contracts

## Versions and boundaries

The first executable data profiles use contract_version 1.1.0 and loader_version country_tb_data_loader_1.0.0. Stable family IDs are country_tb.incidence_history.v1 and country_tb.demography.v1. The separate composition ID country_tb.baseline_data_profile.v1 is a data review object, not the synthetic common population from the architecture.

Milestone 1's 1.0.0 schema remains an archived structural proposal. It is not silently accepted as an executable data profile. New approved mode/sex conventions are specified in schemas/country_tb/v1/data_layer.schema.json; saved-profile restore rejects unsupported revisions. Future incompatible analytical contracts require an explicit migration.

No data-layer code imports screening, TBI history, calibration, projection engines, economics, Starsim or Streamlit. The only shared WHO dependency is the existing pure column contract. Nothing imports the new layer into deployed apps.

## WHO snapshot loader

load_who_snapshot requires an explicit content-based production snapshot ID. It verifies manifest schema and analytical hash; transformation/importer version; importer commit; source checksums; data checksum after the documented LF normalization; required columns; stable ISO3/ISO2/numeric IDs; country/region consistency; years and country/year uniqueness; population denominators; finite nonnegative central/lower/upper triples; uncertainty ordering; missingness; report endpoint and row/country coverage. It never falls back to the test fixture.

The pinned WHO 2025 snapshot has 5,347 rows, 217 countries/areas, 216 with published estimates and global observation coverage 2000–2024. PRK missing estimates remain missing. Historical short territorial series remain short. Report year 2025 and final observation year 2024 are separate fields.

The country index uses stable IDs and supplies display name, WHO region, available estimated years, record years, completeness, warnings and snapshot ID. No undocumented World Bank grouping is invented. Territories and historical areas retain their separate source identifiers.

## Incidence profiles

Modes are national_who_snapshot, subnational_scaled, subnational_user_series and custom_population_series. National profiles retain year, central/lower/upper, cases numerator, population denominator, original rate, quantity and unit; full frozen provenance includes report/final year, snapshot, original/transformed checksums, source, source commit, import date, transformation version, citation and terms.

Subnational scaling retains an unmodified national_parent and its hash, year-complete nonnegative rate ratios, rationale, source and review status. National case counts/population are not presented as local counts: those remain only in the parent. Missing parent estimates remain missing when scaled.

User series preserve supplied numerator, denominator, original rate and quantity. Duplicate/invalid years fail; gaps are listed without interpolation. Rates are per 100,000 person-years. Notifications remain notifications unless an explicit reviewed ascertainment-fraction conversion, source and fraction in (0,1] are supplied. Original values/quantity remain recorded. Numerator/denominator/rate reconciliation uses a disclosed maximum of 0.01 per 100,000 or 1% rounding tolerance; contradictory values outside it fail.

## Demographic snapshots and profiles

Snapshot schema is country_tb.demography_snapshot.v1 and transformation version country_tb_demography_import_1.0.0. Each immutable folder has manifest.json, population.csv.gz and checksummed published-source evidence. The catalogue names immutable IDs explicitly. Load verifies manifest hash, source/evidence/transformed checksums, schema/version, columns, row count, geographic/year coverage, projection boundary and complete-key uniqueness. Invalid numbers, units, sex/status values, conflicts and overlapping intervals fail.

The source CSV columns are geography_id, geography_level, year, estimate_status, projection_scenario, sex, age_lower, age_upper_exclusive, population, unit and population_variant. Source ID/vintage/checksum/review and transformations come from the manifest and are carried into loaded rows/profiles.

Canonical profiles contain baseline year, detailed sex-specific rows, sex totals, combined population/age rows, display bands, source resolution, inactive source/scenario alternatives, published-both checks, gap diagnostics and full provenance. Counts are canonical decimal strings for precision and portable serialization; they are not rounded display strings. See the [source policy](demographic_source_policy.md) for sex/age/year rules.

## Hashing and restoration

Profiles record input/source hash, transformed-file hash, schema/loader/transformation versions and a scientific profile hash. Canonical JSON uses sorted string keys, UTF-8, compact separators, no nonfinite numbers and canonical plain decimal count/rate text.

The explicit display object, profile_hash itself and receipt metadata are excluded from scientific hashing. Set-valued rows/candidates/warnings/alternatives/checks are sorted canonically; transformation sequences retain order. Changing snapshot, source, year, sex, age, population, rate, bounds, scenario, transformation or review changes the scientific hash. Translating a display label or reordering input records does not. File checksums still identify exact frozen bytes.

restore_baseline_data verifies root and component hashes, contract and year/scope, reloads the original WHO/demography IDs, checks source integrity and compares saved active data with those frozen sources. It preserves saved candidate receipts and never reranks against a newer catalogue. Missing old sources fail explicitly. No automatic analysis migration is implemented.

## Developer previews

Run from the repository:

    python scripts/preview_country_tb_data.py --countries AUS VNM PRK

The CLI publishes immutable export folders under abm/output/country_tb_m2a_preview containing incidence uncertainty CSV, sex-specific detailed population, derived combined ages, display-band counts, published-both checks, complete profile JSON and a Markdown summary. It imports no protected workflow and is absent from Streamlit navigation.

## Data limitations

The full original WPP download is retained in the local ignored audit archive, while the checksummed selected native rows and transformed population are committed. Rebuilding from the original requires that archive or an exact-checksum download at the maintainer boundary. Analysis runs need only committed frozen inputs.

The old cache's exact upstream source/generator and sex semantics are unresolved and quarantined. No dataset in this milestone supplies infection pressure, TBI, medical-factor populations, mortality calibration or a projection-ready synthetic population.
