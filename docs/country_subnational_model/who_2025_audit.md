# WHO 2025 data audit

Audit date: 2026-10-27. The sibling gtbreport2025 repository was read-only throughout. This audit selected existing frozen public data; no live WHO dataset was imported or modified.

Milestone 2A adds a separate strict offline loader under engine/country_tb/incidence.py; its [current contract](data_layer_contracts.md) implements the production guards recommended here. The pre-existing shared loader audited below remains unchanged.

## Repository identity and state

| Item | Observed |
| --- | --- |
| Branch | main |
| Commit | 666088cac1e20dd1e9e52016ea857a9c48e0ba8d |
| Origin fetch/push | https://github.com/emmamcbryde/gtbreport2025.git |
| Upstream remote fetch/push | https://github.com/GTB-TME/gtbreport2025.git |
| Tracking branch | origin/main |
| Working tree | Clean; no tracked or untracked changes |
| Report vintage | README identifies Global tuberculosis report 2025, published 12 November 2025 |

The remote named upstream is distinct from the branch's tracking upstream origin/main. No fetch, checkout, R workflow execution, commit or push was performed there. Local remote-tracking refs are observations, not a claim about current server heads.

## Directly inspected data objects

R serialization was read using rdata 1.1.0 in an isolated audit dependency directory under ignored abm/output. R scripts were not executed. Values below are object inspections, not guesses from filenames.

| Relative source in WHO repository | Inspected object and coverage | Suitability |
| --- | --- | --- |
| inc_mort/analysis/est.rda | est; 5,347 rows, 203 columns; 217 ISO3 entities; 2000–2024 | Report-round cross-check only. inc, inc.lo, inc.hi and inc.num/count bounds; HIV/non-HIV components, methods and population also present |
| data/gtb/other/cty.rda | cty; 217 rows, 10 columns; ISO2, ISO3, iso.numeric, multilingual names, WHO status and region | Country/area identifiers; historical entities need effective-date handling |
| disaggregation/output/db_estimates_country_all.Rdata | 240,615 rows, 10 columns; 217 entities; 2000–2024; measure inc, unit num; best/lo/hi and analysed | Age/sex-disaggregated incident counts, including overlapping aggregate strata; not independent additional observations |
| data/gtb/other/pop.rda | pop; 21,487 rows, 28 columns; 217 entities; years 1950–2050 | Total and male/female age-group population. Older estimates and future projections require year/vintage interpretation |
| inc_mort/analysis/rf.rda | rf; 1,393 rows, 11 columns; 199 entities; 2024; sex a; ages a, 15+, 18+ | Attributable incidence counts, not medical-factor prevalence by age |
| data/gtb/other/sdg.rda | sdg; 52,863 rows, 10 columns; 213 entities; years 2000–2025 plus 2030 | Indicator observations, sex and bounds; source-specific age/definition/vintage review needed |

The disaggregation age labels are 0-4, 5-9, 10-14, 5-14, 0-14, 15-19, 20-24, 15-24, 25-34, 35-44, 45-54, 55-64, 65plus, 15plus and a. Sex codes are m, f and a. Select a mutually exclusive set and use matching population denominators to derive rates; do not sum all available rows.

Population columns include male/female 0–4, 5–14, 0–14, 15–24, 25–34, 35–44, 45–54, 55–64, 65+, 15+ and total. These overlap and cannot all be summed. Their precision differs from the single-year local cache.

Risk-factor output includes HIV, alcohol, smoking, diabetes old/new and undernutrition old/new. The generating script 06c-attributable_cases_rf.R explicitly uses population attributable fractions and disease RRs; it notes updated 2022 diabetes prevalence and a 2024 diabetes review. These are not the requested eleven-factor prevalence marginals and must not be substituted for them. SDG/UNAIDS inputs may support selected defaults only after definition and age reconciliation.

Country-reported notification snapshots also exist under data/gtb/snapshot_2025-07-30, including annual_notifs.csv and sub_annual_notifs.csv. They are not incidence estimates or a complete validated subnational incidence loader.

## Frozen public snapshot selected

The installed system is beyond the earlier six-country demonstration. There are two separate kinds:

| Artifact | Finding |
| --- | --- |
| data/who_incidence/snapshots/who-gtb2025-incidence-6997e3b087d2.csv | Production, complete country-level dataset flag true; 5,347 rows, 217 countries/areas, 216 with published incidence, 2000–2024 |
| Corresponding manifest and validation JSON/Markdown | Contract who_incidence_snapshot_v2; source hashes, dictionary, report round, importer version/commit, transformations, terms and validation |
| data/who_incidence/fixture/who-gtb2025-incidence-fixture-1709dd393fec.csv | Test fixture, complete flag false; 185 rows, eight codes: AUS, BRA, GBR, IDN, PHL, ZAF plus edge cases ANT and PRK |

The six ordinary countries in the fixture do not define production coverage. Its manifest coverage.incompleteSeries lists all eight incomplete entities from the full import rather than only the subset, another reason not to use fixture metadata for country coverage.

Selected source vintage is GTR 2025, estimate years 2000–2024, public download/access date 2026-09-24. No new round is selected. Public download URLs are mutable, so retaining the frozen snapshot and original source archive is necessary.

| Integrity item | SHA-256 |
| --- | --- |
| Production CSV, LF canonical bytes | 6997e3b087d26788b21a47d09968eca87b5b03c3cc50644ae4dcdcd860acfd07 |
| Production analytical manifest hash | 32f228fc4338a9a61d56fa4770799c2274a3b9696f814d2304b3e2a3ebdee765 |
| Original public estimates file recorded by manifest | 6849d50d91579eaefa5833aa46d19d34e68dd6408a3768875991013f05586c30 |
| Original dictionary recorded by manifest | 2b0aecd442bffabf2f313d2eaa478ffaea48d09fe78670a136a1956c042011bf |
| Report-round est.rda inspected | 4639b9f878b6a1ec487a18e5fab13a2ca46d1aa5855fbd1d21d055e4929ccd10 |
| Age/sex incidence Rdata inspected | 1bd9eed47005546327bd65db06523b6689b8ec62962ce4311bcbcf9bb55fc1b1 |
| Population Rdata inspected | 1d4fe97d9d558ffa38d29a73f2c3a6a24224ee2ab6422d29429fb44bccaca028 |

PRK has population but no publicly published incidence estimates; retain missingness. The import validation reports 37,279 public values agreeing with the report-round output within publication rounding, and 150 PRK values present in the R output but withheld publicly. The manifest comparedValues field counts 37,429 comparisons including those withheld cells; it must not be described as 37,429 published agreements.

Incomplete territorial series: ANT, CUW, MNE, SCG, SRB, SSD, SXM and TLS. Do not fabricate continuity across changing borders. Bounds are published uncertainty limits, not a joint distribution over years or independent Gaussian standard errors.

The snapshot includes total rate/count estimates and lower/upper bounds; HIV-positive incidence is per 100,000 total population, not per 100,000 people living with HIV. It has no age/sex incidence tables, age demography, complete medical-factor prevalence, life tables or TBI observations.

## Loader and manifest reuse assessment

engine/who_incidence/snapshot.py validates schema fields, CSV header, normalized-LF data checksum, row count and canonical row validation; ordinary use is offline. The importer rejects a subset marked production. Deterministic CSV serialization and analytical manifest hashing are valuable reusable patterns.

However, the loader chooses a fixture when production is absent, validates only part of manifest semantics, does not verify analyticalHash in load_snapshot, and does not bind a caller's expected snapshot ID/hash to a saved analysis. A production country adapter must explicitly require snapshotKind production, complete coverage flag, expected report/snapshot/content hash, valid territorial/year coverage and published estimates. These are proposed changes for a later milestone, not edits to the shared loader now.

The snapshot is suitable as a frozen national disease-incidence input after these boundary checks. It is not a complete production country model dataset.

## Licences and reproducibility

No LICENSE/LICENCE, package lockfile or repository licence statement was found in the WHO checkout. Public GitHub visibility is not a code/data reuse licence. Audit it read-only and do not copy repository code or redistribute its R objects without clarified rights.

The [WHO TB data page](https://www.who.int/teams/global-programme-on-tuberculosis-and-lung-health/data) links public CSVs and terms, and warns that profiles can change after publication. The [WHO dataset terms](https://www.who.int/about/policies/publishing/data-policy/terms-and-conditions) allow public-health use with attribution and no endorsement, restrict modification and commercial promotion, and reserve third-party rights. The existing local terms note leaves written clarification of renamed/subset redistribution open. Preserve that open question; this audit grants no new licence.

Missing reproducibility components include import/load_gtb.R (referenced by incidence and disaggregation scripts), disaggregation/local intermediates, inc_mort/estimates2023/old.rda referenced by the risk script, a locked R/Stan environment and a documented end-to-end rebuild command with all private/external inputs. An older disaggregation input exists, but is not a replacement for missing intermediates. The repository is an analysis archive, not a self-contained production loader.

The main project deliberately keeps original public CSV/dictionary and the cross-check export outside Git. Their recorded hashes help traceability but do not make them retrievable. Milestone 2 should specify a durable, access-controlled or redistributable source archive, archive checksums and import environment, plus frozen public age/sex/demography sources under their own terms. Do not re-download mutable URLs during an analysis run or silently replace this snapshot.
