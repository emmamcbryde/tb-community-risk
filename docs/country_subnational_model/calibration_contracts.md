# Proposed calibration contracts: 1.0.0-design.1

The [machine schema](../../schemas/country_tb/v1/calibration_design.schema.json)
and [six synthetic examples](../../schemas/country_tb/v1/calibration_design_examples.json)
are **design proposals**. They do not replace the frozen 2A snapshots, their
1.1.0 executable profiles or importer versions. Every proposal carries
status design_only and contract_version 1.0.0-design.1; production loaders must
not accept these as validated calibration artifacts. Namespace family suffix
v1 is retained, with an explicit contract revision. Promotion/migration needs
review rather than changing a status label.

Each artifact contains schema_id, contract_version, status, nonempty source
provenance, typed data, scientific profile_hash and optional display labels.
Provenance records source_id, immutable snapshot_id, SHA-256 and transformation
version; the actual future implementation must additionally retain the full
source metadata and parent artifact references. A hash consisting of zeros in
the examples is an invented source reference, not a claim of evidence.

| Family ID | Required scientific contents |
| --- | --- |
| country_tb.infection_pressure_history.v1 | Geography/sex, pressure units, history kind (supplied, conditional scenario, disease-calibrated), mapping and disease-operator IDs, scale/shape normalization, age-exposure policy, backcast hash, explicit future policy, finite nonnegative interval hazards |
| country_tb.historical_backcast_policy.v1 | Method, earliest observed time, required start, parameters and bounds, slope window/distribution where used, join rule, sensitivity IDs and review state |
| country_tb.tbi_observations.v1 | Survey or reviewed override mode, parent hash, observations with IDs, exact age bounds, sex, geography, survey year/date, original numerator/denominator and prevalence, assay/product, positivity threshold with units, sampling frame, uncertainty/design, source/review state |
| country_tb.calibration_policy.v1 | Scenario/TBI/joint mode, free parameters and fixed assumptions, objective and weights, assay operator, disease operator if applicable, parameter bounds, tolerances, identifiability checks, incompatibility action |
| country_tb.calibrated_tbi_age_profile.v1 | 2024 baseline and exact time, target state, recency clock/window, pressure/policy/demography/observation hashes, fit state, exclusions, exact-age sex/geography/scenario rows, probabilities, cumulative hazard, uncertainty, residuals and observation references |
| country_tb.uncertainty_scenario_set.v1 | Parameter draws, structural scenarios, population-variation policy, revision snapshot IDs, explicit RNG metadata or null, weights and interval interpretation |

The calibrated-profile family is deliberately distinct from the older M1
country_tb.tbi_profile.v1 proposal and from the implemented
country_tb.catalytic_diagnostic.v1 result. A future explicit migration must
establish state/recency/survival semantics before either is treated as a
calibrated profile. The diagnostic's total_ever_infected and
recently_first_infected names avoid implying retained viable infection or
most-recent exposure.

## Semantic validation to implement in 2C

The machine schema checks shape, version, units and nonnegative/probability
strings. It does not claim to validate all cross-field scientific constraints.
The following are mandatory before accepting real observations/calibration:

- Histories are unique non-overlapping half-open intervals, cover births and
  survey dates, and contain no hidden joins/interpolation. Supplied risks must
  name their duration and susceptible denominator before hazard conversion.
- A disease_calibrated history or joint objective requires a reviewed disease
  operator, natural-history parameter artifact and historical denominators.
  A null operator is allowed only for supplied pressure or conditional scenarios.
- Backcast methods require method-specific parameters: constant level; early
  trend window/slope/finite limits; historical input quantity and checksum; or
  probability distribution and truncation/feasibility policy. Unsupported
  histories fail; bounds reject rather than silently cap.
- Observation numerator cannot exceed denominator. Original unweighted counts,
  weighted prevalence and effective information must not be confused.
  Reported rounded prevalence disagreement is retained and reviewed. Missing
  counts are null, not zero. Survey dates and geography must match the relevant
  pressure/weights. Sex unknown/both does not become female/male automatically.
- Distinct assay/threshold/design observations are distinct records, even when
  their prevalence matches. Specify uncertainty level/method, covariance or
  design effect where available; indeterminate results have an explicit
  denominator/exclusion policy. TST and IGRA cannot share an identity operator
  silently.
- An override needs reviewed provenance and parent hash; it is a new anchor
  artifact. Exact incompatible targets return an explicit failure. No source
  mutation or displayed-prevalence replacement can masquerade as calibration.
- \(U+R+M=1\), \(R+M=T\), \(H\ge0\), all probabilities in [0,1].
  Clock and infection target must match the transition model. Rounding and
  aggregation tolerances are stated, not repaired by clipping.
- No categorical "recent prevalence" input; window and infection histories
  generate recency. Fine-age probabilities aggregate using appropriate
  age/sex population counts. Missing/open-tail inference has explicit coverage.
- Scenario weights are nonnegative and sum to one when probabilistic; do not
  assign probabilities to unweighted scenarios without review. RNG metadata
  specifies algorithm, seed, stream/draw ID and sampling implementation version.

## Canonical serialization and provenance

Use existing data_common canonical JSON: UTF-8, sorted object keys, compact
separators, no NaN/infinity. Normalize scientific decimal values through
numeric_text; serialize finite decimals as plain strings without redundant
zeros. Sort set-valued rows, source_files and warnings by canonical contents;
preserve scientifically meaningful sequences such as ordered transformations.
The history prototype explicitly sorts intervals before hashing.

profile_hash is SHA-256 of the scientific payload, excluding only display,
profile_hash and receipt fields. Display names belong under display, not in
source/geography IDs. Source snapshots, sex, ages, pressure, year, observations,
assay, bounds, mortality/migration assumptions, transformations, RNG policy
and parent hashes are material. Ordering or label changes do not alter
scientific identity; material transformations do. Raw byte checksums and
manifest analytical hashes remain separately visible.

Restoration resolves explicit original snapshot IDs and verifies hashes. It
never resolves "latest", reranks demographic candidates or rewrites a prior
result. Derived calibration histories and overrides get new immutable IDs.
An eventual common baseline references exactly one calibrated age-profile
hash and one configuration/scenario set; both engines consume those same
inputs and infection-history semantics.
