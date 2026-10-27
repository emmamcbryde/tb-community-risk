# Proposed country TB versioned contracts

## Status and shared conventions

Milestone 2A's first executable incidence/demography profiles use revision 1.1.0, documented in [data_layer_contracts.md](data_layer_contracts.md) and schemas/country_tb/v1/data_layer.schema.json. The 1.0.0 material below is the archived Milestone 1 proposal. It is preserved for review and is not automatically migrated or accepted by the production data loaders.

Contract release 1.0.0 is a review proposal. [contracts.schema.json](../../schemas/country_tb/v1/contracts.schema.json) is a JSON Schema Draft 2020-12 document with seven named definitions and a root union. Its URI is urn:country-tb:contracts:1.0.0; validation is local and requires no network. [review_examples.json](../../schemas/country_tb/v1/review_examples.json) contains seven invented structural examples with placeholder hashes and failed readiness, not country data or a calibrated population.

The JSON Schema checks structure, enums, required fields, ranges and value-state nullability. Cross-field scientific checks below are normative proposed behavior; only the isolated prevalence, age reconciliation and association helpers are executable in Milestone 1. Schema validation alone never confers analytical readiness. Full cross-artifact validation is Milestone 2 work.

Each artifact requires schema_id, contract_version, artifact_id and provenance. Versions use semantic versioning: changing units, age boundaries, state meanings or infection-clock semantics requires a new major contract. Adding optional fields requires a minor version with an explicit migration; reject unknown fields in v1 rather than dropping them.

Each provenance entry contains source_id, locator (stable URL/DOI/relative archive reference), source_vintage, SHA-256, transformation_rules, review_status and terms_status. Keep acquisition date, report year, source revision and observation/reference year distinct. Reviews are not_reviewed, reviewed or rejected. A source must resolve to frozen bytes or a frozen entered-value record; a citation alone is insufficient for a reproducible numerical default.

Use a stable semantic geography_id separately from artifact identity. Proposed IDs are who:ISO3 for a national WHO area, iso3166-2:CODE for a supported subnational boundary, and custom:UUID for a defined population. Parent links and boundary vintages are mandatory where applicable. Former entities retain distinct IDs rather than being merged into present countries.

For content addressing, calculate the digest of the complete analytical artifact excluding its artifact_id, then set artifact_id to schema_id:sha256:DIGEST. The caller constructs that hash-free payload explicitly; content_hash does not remove fields. Store audit receipts and timestamps separately. All provenance and review decisions in the analytical artifact affect its identity. Parent hashes refer to the same identity payload convention. Raw data hashes always identify source bytes with a documented newline policy.

The prototype canonical serializer sorts string keys, escapes Unicode deterministically, uses compact Python JSON number formatting, disallows NaN/infinity and non-JSON types, and retains array order. It does not claim cross-language RFC 8785 equivalence. Integer 1 and floating 1.0 are distinct encodings; normalize types through a specified importer and record its version. Source text precision is preserved in frozen WHO CSV; numerical transforms use full supported precision, never display rounding. A later non-Python implementation must match canonical golden bytes or introduce a new serialization version.

Age intervals are [age_lower, age_upper_exclusive), in completed years, with null upper bound for the oldest open group. Continuous individual age is age_years. Sort distribution rows by sex, age lower bound, then upper bound; incidence by year; medical rows by factor ID, sex and age; population by stable person_id. Validate unique keys. Within each disjoint sex stratum age bands must form an ordered partition beginning at zero. Never add all-sex aggregate rows to male/female rows.

## Geography

Schema ID: country_tb.geography.v1.

Fields: geography_id, kind (national/subnational/custom), name, country_iso3, parent_country_iso3, reference_year, population_definition, boundary_vintage, plus the common envelope.

National records need a valid source country/area code and no parent-country link. Subnational/custom records need a parent country when the model relies on its national WHO series; multinational populations require an explicitly revised contract, not a guessed parent. Preserve target population definition, boundaries and effective dates. Reject unknown codes rather than performing substring matching.

## Incidence history

Schema ID: country_tb.incidence_history.v1.

Fields: geography_id, mode (national_who/local_replacement/national_scaled), quantity fixed to tb_disease_incidence, report_year, source_snapshot_id, source_data_sha256, scale_factors, rows and provenance. Each row has year, population denominator, incidence_per_100k, incident_cases, case_definition and ascertainment_notes. Each measurement carries state (value/missing), value, lower, upper and unit.

Rates use per_100000_person_year; counts use cases_per_year. Values and bounds must be finite and nonnegative, with lower <= estimate <= upper when bounds exist. A missing measurement has null numerical fields; measured zero is a value. Require both uncertainty bounds or neither, increasing unique years and explicit missing-year handling. Check count/rate/population consistency within published rounding, not exact equality of rounded WHO numbers.

National modes must bind the expected production WHO snapshot and its hash; national_scaled requires complete year-specific nonnegative rate ratios. The local source is explicitly independent or derived, never blended invisibly. Preserve the unscaled parent series and derive a new artifact. Reject quantities labelled force of infection or notifications in this disease-incidence contract unless a separate reviewed observation conversion has been recorded.

## Demography

Schema ID: country_tb.demography.v1.

Fields: geography_id, geographic_status, reference_year, population_size, source_population_size, native_age_resolution, sex_availability, allocation_rule, oldest_age_rule, rows and provenance. Rows hold age bounds, sex and count. geographic_status distinguishes national_default, subnational and user_defined.

Source counts sum to source_population_size; a documented scaling operation maps them to population_size, the represented target population. Preserve both totals and reject unexplained discrepancies. Derived age proportions are count divided by the relevant all-sex or sex-stratum total. Positive population totals and finite nonnegative counts are required; no silent global or uniform fallback.

The user-facing bands are 0–4, 5–14, 15–44, 45–64 and 65+, while the loader retains supported finer ages. The prototype aggregates only whole source age cells and rejects splitting a coarse cell. A source oldest interval cannot be split without a reviewed tail model. Source vintage, checksum, transformation and review remain in provenance after any aggregation.

## TBI profile

Schema ID: country_tb.tbi_profile.v1.

Fields: geography_id, reference_year, mode (model_estimated/observed_anchor), denominator, calibration and rows. Calibration identifies status (compatible/partial/failed/not_fitted), parameter hash, natural-history ID, explicitly selected backcast and provenance ID, historical start year and positive recent_window_years.

Rows contain age/sex, uninfected, recent, remote, total, lower_total, upper_total, uncertainty_kind, observed_total, residual, standardized_residual and derivation. All model prevalences are proportions. Require recent + remote = total and uninfected + total = 1 using a stated numerical validation tolerance for computed floating outputs (proposed 1e-12); never repair or clip mismatches. Apply bounds only to their stated estimand; uncertainty draws must preserve the state simplex, rather than drawing each state independently from marginal intervals.

Observed anchors retain their own survey age/sex bands, source weights, assay definition, dates and measurement-model provenance; the displayed model profile may use different bands only through explicit weighted aggregation. Use residual = model prediction on the observation scale minus observed value. Standardized residuals require an identified variance/covariance and are null otherwise.

Observed_anchor mode needs at least one valid observation and a calibration attempt. It never substitutes observations for the fitted displayed model values. Incompatible observations set failed/partial status with diagnostics. Missing observed data in model_estimated mode are null, not zero. Additional survey metadata and a calibration-report artifact are required before an analytical anchor loader is production-ready.

## Risk factor prevalence

Schema ID: country_tb.risk_factor_prevalence.v1.

Fields: geography_id, reference_year, denominator (whole_population/non_active_tb_population), definition_registry_id, rows and provenance. The registry is tstin3d_candidates_m1 for the eleven IDs in the [evidence table](tstin3d_evidence.md). It is an evidence inventory, not an approved effect set. Transform whole-population marginals explicitly before using them in the non-active-TB association contract; do not equate those denominators.

Each age/sex/factor row requires state, value_percent, lower_percent, upper_percent, reason, source_ids, review_status and geographic_status. Source IDs must resolve in provenance. Percentages are visibly labelled and remain 0–100 inclusive; internal proportions are percentage divided by 100 without display rounding. Optional bounds must occur together and enclose the value within 0–100.

| State | Numerical fields | Meaning and behavior |
| --- | --- | --- |
| value, including 0 | Required percentage; optional pair of bounds | Measured/entered value; zero means explicitly absent in that stratum |
| missing | All null | Unknown; blocks a reviewed included-factor baseline |
| excluded | All null, reason required | Explicit analysis decision, distinct from measured absence; retain the underlying observation in its source artifact |
| not_applicable | All null, reason required | Definition does not apply in this stratum; requires review, especially in children |

Recency, close_contact and AIDS as a separate HIV duplicate are invalid factor IDs. Every included factor must cover the target age/sex support once, without duplicate cells. Excluded/not-applicable factors do not contribute to the chosen any-factor union, but the report discloses their omission. Do not assign unknown factors as absent.

The eventual editor validates a staged edit before applying it, highlights the offending cell, keeps valid source precision and never changes missing to zero through an empty text box. Country defaults and user overrides are separate artifacts with explicit conflict resolution.

## TBI and risk factor association

Schema ID: country_tb.tbi_risk_association.v1.

Fields: geography_id, denominator fixed to non_active_tb_population, method fixed to frechet_piecewise_overlap_v1, factor_allocation_policy, age/sex rows, overall joint, within_age_covariance, age_driven_covariance and provenance. Row fields include demographic weight, requested_alpha, TBI/any-factor marginals, minimum_overlap, maximum_overlap, achieved overlap, phi, phi_status and four cells.

See [association.md](association.md) for equations. Validate the complete joint against the marginals and alpha, not merely field ranges. Require nonnegative cells and correct row/column sums. Phi is null only for degenerate marginals with an explicit status. Never report an arbitrary requested Pearson value as achieved. Overall results derive from demographic weights and overlap; there is no single overall requested alpha.

The factor-allocation policy must state whether q is derived under independent age-patterns or observed and jointly reconciled. Version and hash any pattern table, pairwise constraints or copula/IPF model. A feasible aggregate joint alone does not validate individual-factor allocation.

## Common baseline population

Schema ID: country_tb.baseline_population.v1.

Manifest fields: geography_id, reference_time (calendar years), represented population_size, n_records, parameter_draw_id, realization_seed, rng_algorithm, recent_window_years, infection_clock_semantics, parent hashes, readiness and people. The schema uses inline people for the small structural example; future columnar storage must preserve the same row contract, stable sort and logical hash, with a storage manifest and physical-file checksum.

Every person contains:

- person_id, population weight, continuous age_years and sex or explicit unknown;
- TBI state (uninfected/recent/remote/unresolved/not_applicable_active_tb), first and most recent infection time, time_since_infection_years;
- all medical-factor states and the any_risk_factor indicator;
- progression_risk policy, multiplier or unresolved representation, effect_measure, reference population and review status;
- mortality_input_id, screening_eligibility and eligibility_reason;
- clinical active-TB/prior-TB status, treatment status and infectiousness or explicit unknown.

Manifest hashes bind configuration, data, calibration, geography, demography, TBI profile, medical prevalence and association. A reference to mortality data must itself be covered by the data hash. Both static and dynamic initialization must return exactly the same baseline identity.

Cross-field checks require unique person IDs, n_records matching row count, positive weights summing to represented population, nonnegative age and time, infection events no earlier than birth or later than reference_time, first event no later than most recent event, time_since matching the last event and derived recent/remote classification. Uninfected people have null infection times. Clearance/reversion must introduce reviewed state semantics before a previously infected person can be labelled uninfected.

Active-TB rows use not_applicable_active_tb in the TBI field rather than counting active disease as latent infection or unknown. Retain known historical event times separately. TBI profile and association aggregates condition on non-active clinical rows; the total population and mortality account for everyone. Unknown clinical partition prevents readiness. A later treated/recovered state model needs an explicit reviewed transition into an infection-state definition.

Any-risk is true exactly when an included factor is present, false when all included factors are known absent/non-applicable/excluded, and null when remaining uncertainty prevents classification. The progression representation and age permissions must be compatible with the selected evidence registry. Do not assign a known multiplier to unreviewed child effects.

Readiness is separate for static and dynamic engines. Unknown clinical state, unresolved mortality, unapproved natural history or incompatible calibration blocks the relevant engine. Dynamic readiness additionally needs reviewed active disease, infectiousness and transmission initialization; those are not identified by disease incidence alone. The one-row structural example deliberately fails these checks and is not a synthetic population sample.

## Executable coverage in this milestone

Pure helpers prove deterministic JSON/hashing, percentage/state validation, age-band aggregation, Fréchet mapping, phi and demographic aggregation. Focused tests validate all seven structural examples and provenance round-trips, invalid inputs and isolation from protected apps/runners. They do not prove natural-history calibration, epidemiological validity, a complete country data loader or an executable population generator. Those remain explicit future work.
