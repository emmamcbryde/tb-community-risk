# Frozen demographic source policy

## Deterministic hierarchy

For Australia: reviewed user national/subnational/custom data, ABS, UN WPP, then documented fallback. For other countries: reviewed user data, UN WPP, then documented fallback. There is no undocumented cache fallback. Candidates must match the requested geography, source year and scenario. Source names never substitute for geography identifiers.

The resolver returns candidates considered, eligibility, priority, rejection/selection reasons, source year, age completeness, warnings and the selected immutable snapshot. Equal-priority eligible vintages require an explicit snapshot choice; neither directory order nor the latest filename resolves a tie. An explicitly restored snapshot bypasses reranking but must still be present and valid.

User demographic CSVs use the documented eleven canonical source columns in persons. Their original file, source URL/citation, vintage, terms and explicit review status are preserved. Import does not promote an unreviewed user file to reviewed. Subnational/custom geographies remain distinct from the parent national population.

## Frozen sources selected

Australia uses the [ABS June 2024 release](https://www.abs.gov.au/statistics/people/population/national-state-and-territory-population/jun-2024), Table 59, published 12 December 2024. This supplies exact 2024 estimated resident population at 30 June, single ages, female/male and published Persons totals. Its original workbook is retained.

Other countries use [UN World Population Prospects 2024](https://population.un.org/wpp), the original Medium single-age/sex CSV, not a later update. The frozen selection includes 237 countries/areas and years 2024–2030. All 215 WHO entities with 2024 records have a WPP counterpart. Historical WHO entities ANT and SCG lack a corresponding 2024 demographic entity and are not merged into successors.

WPP publication units are thousands at 1 July. Multiplication by 1000 is exact. The open 100+ age group remains open. WPP estimates end in 2023; its 2024 baseline and later selected years are projections under Medium. Profiles expose that status prominently.

ABS [copyright terms](https://www.abs.gov.au/website-privacy-copyright-and-disclaimer) provide CC BY 4.0 with listed exceptions. The UN official download catalogue identifies CC BY 3.0 IGO and the WPP citation. Source terms are recorded without suggesting agency endorsement.

## Year and scenario resolution

Default target is 2024. All installed baseline examples use exact source-year 2024, with no interpolation, extrapolation, nearest-year substitution or count adjustment. ABS is an estimate; WPP is a projection.

If exact coverage is unavailable, the default resolver fails. An optional explicitly reviewed nearest-year proxy records the actual source year in every row, target baseline year separately, year displacement, unchanged counts and the rule/rationale. Equidistant choices fail. Count adjustment is zero because this method does not estimate population growth. Interpolation and extrapolation are not implemented; their absence never triggers implicit carry-forward.

After 2024, the default scenario is medium where supplied. Alternatives remain stored and distinguishable, and an explicit scenario argument can select them. The installed WPP snapshot ends in 2030; later years fail until an explicitly imported snapshot supplies them. Australia can use WPP after its ABS 2024 coverage, with the source switch visible in candidate diagnostics; ABS estimates are not relabelled as projections.

## Sex and age representation

Canonical sex values are female, male, both and unknown. Female/male are retained separately when available. Published both-sex records remain in frozen source files and reconciliation checks, never in their population sum. Both-only data remain both-only and visibly lack sex-specific modelling. Unknown means the source genuinely did not specify sex. No 50:50 split or subtraction-derived sex is generated.

ABS requires exact male+female equality with published Persons. WPP permits an absolute one-person rounding discrepancy plus relative 1e-9 times the published cell. Every achieved difference is reported. This allows publication rounding; it does not average or repair conflicts. The full WPP 2024 country audit has maximum age-cell differences of at most one person.

Internally retain source single years and an open 100+ tail. Display bands are [0,5), [5,15), [15,45), [45,65), [65,infinity). Decimal count arithmetic gives exact sex-specific and combined reconciliation. Coarse intervals crossing a display boundary fail aggregation rather than being split uniformly. Gaps remain explicit diagnostics. Different male/female age supports cannot silently produce a combined age profile.
