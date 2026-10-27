# Frozen demographic data

The explicit catalogue selects two immutable snapshots. No loader downloads data or uses the legacy cache as a fallback.

| Snapshot | Source | Frozen coverage | Published sex records |
| --- | --- | --- | --- |
| abs-erp-jun2024-07ca11bcb092 | ABS June 2024 release, Table 59, released 12 December 2024 | Australia, 2024, single ages 0–99 and open 100+ | Female, male and Persons for audit |
| un-wpp2024-medium-8028a7a3e163 | UN World Population Prospects 2024, original Medium CSV | 237 countries/areas, 2024–2030, single ages 0–99 and open 100+ | Female, male and total for audit |

ABS counts are persons at 30 June. WPP native counts are thousands at 1 July and are multiplied by 1000 exactly. WPP classifies 2024 onwards as projections; the 2024 baseline is not relabelled as a historical estimate.

Both snapshots preserve published both-sex evidence, but canonical profiles sum only female and male when both are available. The WPP 2024 reconciliation CSV records age-cell rounding discrepancies for every country/area; ABS reconciliation is exact.

Each snapshot has a versioned manifest, transformed CSV checksum, original-source checksum, transformation list, citation, source terms and import date. ABS includes the original workbook. WPP includes all original fields/values for the selected country/year rows in a compressed source extract. The original 66,954,143-byte WPP archive is preserved under ignored abm/output/country_tb_m2a_audit/sources; its SHA-256 and official URL are recorded in the manifest. It is not duplicated in Git.

legacy_audit preserves a deterministic deduplicated quarantine view and the full duplicate-group audit. The untouched original remains data/population_age_latest.csv. Quarantine records are not eligible for the source catalogue: the old external revision/generator and sex dimension are unrecorded, its year is 2021, and it includes regional aggregates and a custom South Australia population.

See [source policy](../../docs/country_subnational_model/demographic_source_policy.md), [forensic report](../../docs/country_subnational_model/duplicate_forensics.md) and [data contracts](../../docs/country_subnational_model/data_layer_contracts.md).
