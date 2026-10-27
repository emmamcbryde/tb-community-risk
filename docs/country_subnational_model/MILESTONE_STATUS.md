# Country and subnational TB milestone status

Updated 2026-10-27. This file replaces the Milestone 1 status; historical versions remain in Git.

## Current milestone

Milestone 2A frozen incidence and demographic data loaders is implemented and validated, ready for a data-layer recheck. No infection-pressure/TBI model, individual synthetic population, projection, intervention, economics or production app is implemented.

## Trusted data layer

- Strict explicit-ID offline WHO loader: GTR 2025, observations 2000?2024, baseline 2024, 5,347 rows and 217 countries/areas. No fixture fallback.
- Country/territory index retains stable identifiers, display names, WHO region, years, completeness and warnings.
- Incidence modes: national_who_snapshot, subnational_scaled, subnational_user_series, custom_population_series. User numerators/denominators/rates and notification status are preserved.
- Demographic complete key retains source/vintage, geography/level, year, estimate/projection status, scenario, sex, age bounds, unit and population measure.
- Source ranking: reviewed user, ABS for Australia, UN WPP, documented fallback. Tied vintages require explicit choice.
- Australia uses the exact ABS 2024 national estimate. Viet Nam and other supported countries use UN WPP 2024 Medium; its 2024 value is explicitly a projection.
- Female/male remain separate. Published both-sex rows are audit-only when both sexes exist. Both-only/unknown sources remain visibly limited.
- Immutable imports, source/evidence/transformed checksums, deterministic scientific hashes, exact decimal aggregation and original-snapshot restoration are implemented.
- Medium is the default supplied projection scenario after 2024. Frozen WPP coverage ends in 2030; later years fail without an explicit new import.

## Duplicate forensic result

Original file unchanged: 21,917 rows, SHA-256 e51101fa3f5054b9aeb11999d38f6b39a3b0ab4625daa00c5c7fd52e77e535d4.

There are 1,111 duplicate groups involving 2,323 rows and 1,212 excess rows. Ten countries contribute 101 excess rows each; VCT contributes 202. Every observed duplicate is an identical full-row repeat, with zero numerical discrepancy. Other categories have no established groups.

The quarantined view has 20,705 unique legacy keys. The exact upstream generator/revision and sex semantics remain unresolved. No cache record is promoted to an active 2024 profile.

## Validation and outputs

139 focused tests passed in 116.65 seconds, including new data/duplicate/provenance/protection checks and relevant existing WHO/profile tests. Compilation, preview schemas/hashes, local links and Git whitespace checks passed. See milestone2a_validation.md for commands, earlier corrected failures and source/profile hashes.

Read-only preview folders under abm/output/country_tb_m2a_preview:

- AUS-f61fc25afb90
- VNM-cc13e19cc3c1
- PRK-de1392b6cfc8

Baseline totals are 27,204,809; 100,987,682; and 26,498,835 respectively. PRK incidence remains missing and visibly unavailable. Machine-readable duplicate and WPP reconciliation audits are committed under data/demography.

## Contracts and protection

Executable profiles use contract revision 1.1.0 and loader country_tb_data_loader_1.0.0. Family IDs are country_tb.incidence_history.v1, country_tb.demography.v1 and country_tb.baseline_data_profile.v1. The archived 1.0.0 proposal is not silently migrated.

Protected worktree fingerprints and refs matched the starting receipt before commit. Frozen SA Health artifacts, unfinished explicit-TBI changes, the generic model, Starsim, WHO sibling, MATLAB and deployed apps were not modified. The only new CLI is an isolated data export. The final response records commit, push and post-push protection/clean-state checks.

## Unresolved issues and Milestone 2B

Data limits: legacy cache lineage; year rules beyond the reviewed nearest-year proxy; coverage beyond 2030; additional variants/subnational ABS snapshots; and the existing WHO terms clarification. The full original WPP archive is preserved locally; selected native rows and transformed inputs are versioned in Git.

Recommended Milestone 2B: independently recheck frozen source/year semantics, then approve child/adult natural-history evidence, recency/reinfection/clearance, pre-WHO histories, assay/TBI anchors and the calibration identification design. Preserve sex strata. Begin no calibration or intervention engine until explicitly authorized.
