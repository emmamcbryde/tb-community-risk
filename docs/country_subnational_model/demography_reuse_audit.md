# Demographic loader reuse audit

Audit basis: general-community commit 24b632390ce75e24c15088309e45960a990e66ad. No existing component was edited.

This is the historical Milestone 1 audit. [Milestone 2A forensics](duplicate_forensics.md) now classifies all 1,212 excess records as identical stored-row repeats, with unresolved upstream lineage. The [current source policy](demographic_source_policy.md) uses new frozen ABS/WPP data rather than this cache.

| Component | Observed source and behavior | Reuse decision |
| --- | --- | --- |
| app/general/country_panel.py | Country/area dropdown from frozen WHO snapshot; preview/apply conflict controls and descriptive incidence chart | Reuse interaction design later; no page import now |
| engine/profiles/country.py | Country identity, incidence and contextual national population; explicitly leaves age distribution, TBI and medical factors unchanged | Reuse conflict/provenance concept; it is not a country demography loader |
| app/demographic_profile.py | Resolves APY config/data into age and risk summaries, including broad age bands and restoring APY defaults | Audit only; do not import APY defaults or config resolver |
| ui/static_ui.py, ui/dynamic_ui.py, ui/update.py | Legacy local cache loaders, manual ISO3 controls/uploads, five-year displays; fallback is uniform ages 0–100 | Do not reuse fallback or initialization; interface ideas only |
| ui/dynamic_age_distribution.py | Five-year expansion distributes each source band uniformly; invalid values are coerced to zero/clipped; missing distributions fall back to uniform ages, and 100+ becomes age 100 | Do not reuse validation or tail assumptions; rederive aggregation under the new interval contract |
| engine/demography_fetch.py | Live World Bank total population request plus mutable OWID age CSV; country filter searches Entity text for code; latest-year selection/fallback | Not production-ready: no timeout, stable vintage, checksum, exact identifier join or reviewed provenance |
| load_custom_age_distribution in that module | Reads counts/proportions and divides by their sum | Reimplement strict validation later; rejects neither all invalid inputs nor ambiguous/overlapping ages before normalization |
| data/population_age_latest.csv | Local single-year age cache attributed to OWID by UI comments, reference year 2021 | Candidate raw material only; source publication/vintage/units/terms and duplicate keys require review |
| data/age_dist.csv and abm/default_age_distribution.csv | Existing community/APY distributions | Protected/context-specific evidence; no national defaults |
| WHO pop.rda | Historical/future national age/sex groups 1950–2050 | Evidence candidate; reuse rights and public-source equivalent need resolution |

## Cache findings

The local cache contains 21,917 rows, 205 distinct ISO codes and ages 0–100. It has no sex column. There are 1,212 duplicate (iso_code, year, age) keys, producing 20,705 unique keys. BOL, BRN, CIV, LBY, MMR, RUS, KOR, TWN, VEN and VNM each have 202 rows, while VCT has 303 rather than 101. Alias/name duplication is a possible explanation, not a validated deduplication rule. Do not sum these rows or choose one silently.

SHA-256 of the current file bytes: e51101fa3f5054b9aeb11999d38f6b39a3b0ab4625daa00c5c7fd52e77e535d4. Its labels do not establish whether age 100 means exactly 100 or an open 100+ tail. The cache has no manifest establishing the UN WPP revision, extraction date, transformations or source licensing. National all-sex data cannot represent a subnational or selected community by merely changing its label.

The old live fetch requests World Bank SP.POP.TOTL but does not use that total to construct its returned age profile. It then normalizes OWID values and can silently switch to cached data or a global/uniform fallback. The legacy local UI loader does not select a reference year. These behaviors are incompatible with reproducible calibrated exposure histories.

## Proposed loader behavior

Build a separate offline loader under engine/country_tb in a later milestone. Require the demography contract described in [schemas.md](schemas.md), explicit geography and reference year, unique disjoint age/sex cells, finite nonnegative counts, positive total and reconciliation with target size. Keep source counts and transformed proportions separate. Do not erase a count discrepancy by automatic normalization.

Use exact identifiers, display canonical country names, preserve source aliases, and distinguish national default, subnational estimate and user-defined population. Changing the country previews demographic consequences and asks for an explicit resolution only if it would overwrite a user choice in the eventual UI.

Retain single-year resolution where trustworthy. Reconcile counts into 0–4, 5–14, 15–44, 45–64 and 65+ for display, permitting expansion. A five-year source band can aggregate into a wider band; a broad band cannot be split without a reviewed allocation rule. Preserve the oldest open interval and uncertainty about its tail. Sex-specific counts must sum to an explicitly defined all-sex total without counting aggregate cells twice.

Freeze a reviewed national demographic snapshot, preferably a specified UN WPP revision/public export or national statistical source. For a subnational/custom population, require local age counts and scope, or a visible national-structure assumption with sensitivity. Mortality/life-table provenance and historical migration are separate inputs, not inferred from this cross-sectional age distribution.
