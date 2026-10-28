# Milestone 2C executable contracts

Executable revision **1.0.0**, implementation
`country_tb_catalytic_calibration_1.0.0`, uses canonical UTF-8 sorted-key JSON
and SHA-256. Frozen 2A revision 1.1.0 contracts and 2B design proposals remain
unchanged. Design-only artifacts require explicit migration; changing a
status label does not make them executable.

[Machine schema](../../schemas/country_tb/v1/catalytic_calibration.schema.json)
checks versions, structure, units, state ranges, semantics and exclusions.
Pure Python APIs additionally enforce cross-field history coverage, source
references, review, age/sex weights, assay informativeness, target feasibility,
state/count reconciliation and restoration from original embedded inputs.

| Family | Creator / validator |
| --- | --- |
| `country_tb.infection_pressure_history.v1` | `pressure_shape.apply_backcast / validate_history` |
| `country_tb.historical_backcast_policy.v1` | Embedded, separately hashed by `apply_backcast` |
| `country_tb.tbi_observations.v1` | `tbi_observations.observations_profile` |
| `country_tb.tbi_calibration_policy.v1` | `catalytic_calibration.calibration_policy` |
| `country_tb.calibrated_tbi_age_profile.v1` | `tbi_baseline.baseline_tbi_profile / restore_tbi_profile` |
| `country_tb.uncertainty_scenario_set.v1` | `tbi_baseline.uncertainty_scenario_set` |

Supporting families include pressure-shape proxy, assay operator, fitted
scale, catalytic age result, reviewed TBI override, preview diagnostics,
export receipt and unavailable-pressure diagnostic.

The state fields are `never_infected`, `recent`, `remote_only`, `any_tbi` and
`prior_remote_plus_recent`. `most_recent_infection_time` contains conditional
calendar-event and time-since summaries. The progression clock resets at the
most recent reinfection without any progression-engine connection.

Each profile row retains population/geography/2024 timing, shape identifier,
alpha, remote/recent hazards, all probabilities, expected lifetime events,
conditional event-time summaries for any/recent/remote-only, expected state
counts, historical policy, recent window, calibration source, scenario,
reinfection/clock policies, exclusions and provenance/scientific hashes.
Required camel-case probability aliases are exported alongside canonical
state names. Unresolved open tails retain population with null estimates.
Sex, derived-both, display-band and complete-population reconciliations are
separate from published-both checks.

Observation records require geography and level, age range, sex, survey year
and evaluation time, numerator/denominator (jointly present or jointly null),
prevalence, uncertainty, assay/threshold, sampling frame/design, source,
review, observation policy, weight source and explicit age/sex strata.
Direct latent-TBI and binary assay positivity are distinct observation modes.
Performance provenance is retained separately for sensitivity and specificity.
The immutable source receipt is checked independently from scientific hashing.

A calibration policy selects `overall_anchor`, `binomial` or
`weighted_least_squares`, one free alpha, finite bounds, tolerances, explicit
weights and compatibility policy. Supplied-scale mode uses empty observation
weights and mandatory scale provenance. Observations are evaluated at their
survey time; 2024 weights are not copied backwards without a declared survey
weighting assumption. Objective, residuals, feasible range, convergence,
boundary status and identifiability limits are retained.

The minimal workflow is: frozen incidence/demography -> reviewed proxy ->
explicit lifetime backcast -> supplied scale or reviewed observations/assay
and one-scale fit -> conditional 2024 age-sex profile -> immutable export.
The CLI is `python scripts/preview_country_tb_calibration.py`. It imports no
production app or engine and labels every demonstration diagnostic.

Scientific hashing excludes only established `display`, `receipt` and
`profile_hash` fields, and sorts known set-valued row collections. Creators
canonicalize sequence-sensitive strata/curves before sealing. Original
observation receipts retain their own canonical source-record SHA-256; input
row order can change the receipt without changing scientific identity.
Changes to shape, source vintage, backcast, assay, alpha, sex exposure,
age policy, time/window or override produce a different scientific hash.
