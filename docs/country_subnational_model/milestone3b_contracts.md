# Milestone 3B executable contracts

The new [Draft 2020-12 schema](../../schemas/country_tb/v1/progression_forward.schema.json)
uses version `1.0.0`, implementation `country_tb_progression_forward_1.0.0`.
No frozen 2A–3A artifact is revised. All model modules use the standard library.

| Required family | Constructor/module |
| --- | --- |
| `country_tb.progression_natural_history.v1` | `progression_natural_history.natural_history` |
| `country_tb.progression_effect_application.v1` | `progression_policies.effect_application` |
| `country_tb.reinfection_progression_clock.v1` | `progression_policies.reinfection_clock` |
| `country_tb.competing_mortality.v1` | `progression_policies.competing_mortality` |
| `country_tb.future_exogenous_pressure.v1` | `progression_policies.future_pressure` |
| `country_tb.active_tb_forward_operator.v1` | `progression_forward.forward_operator` |
| `country_tb.active_tb_projection_result.v1` | `progression_forward.project_active_tb` |
| `country_tb.incidence_validation_comparison.v1` | `progression_validation.incidence_comparison` |

Supporting families: `country_tb.progression_baseline_cohort.v1`,
`country_tb.active_tb_monte_carlo_validation.v1`,
`country_tb.active_tb_forward_diagnostic_set.v1`,
`country_tb.active_tb_forward_export_receipt.v1`.

Constructors enforce semantics beyond JSON structure: monotonic finite risks,
net-risk interpretation, exact age/sex applicability, explicit missing curves,
non-overlapping survival/pressure intervals, reviewed provenance, complete source
identities, first-event absorption, recency-window consistency, no inverse fitting
and no automatic future extrapolation. Source hashes do not substitute for
scientific review. National mortality additionally requires licensing review and
an immutable source hash. Schema validation is fully local/offline.

`baseline_cohort` validates an explicit partition with exact supplied ages and
infection-time support, exclusive baseline TBI states, highest-RR attribution,
baseline active-TB exclusions and unresolved reasons. `cohort_from_allocation`
also records effect age applicability: finite reviewed or pediatric scopes cannot
be silently extended over follow-up. `cohort_from_allocation`
preserves source 3A artifacts and their hashes. Catalytic time support comes from
the source history; direct user profiles need separately reviewed conditional
age/time support. A mean time or display band cannot fill this gap. Open tails
remain unresolved without reviewed ages. Model-derived and supplied source values
remain available for audit.

`forward_operator` requires all curve/effect/clock/mortality/future policies,
mode, numerical step, horizons, review metadata and diagnostic status explicitly.
No-mortality and provisional/synthetic curves cannot become reviewed outputs.
Partial reset is single-event calculation-only. Central and sensitivity policies
have separate scientific hashes. There is no implicit RR cap or country curve.

`project_active_tb` returns instantaneous baseline incidence and cumulative
first-episode projections with coverage denominators, active-TB exclusions,
survival, competing deaths, integrated disease-free person-time, decompositions,
common-horizon comparison counts and extreme-risk concentration. `complete_result_available`
means the mathematical projection is covered; `complete_country_estimate` remains
false pending independent validation and unresolved historical assumptions.

`monte_carlo_validate` consumes the same hashed baseline/operator/projection,
requires an explicit RNG and bounded sample size, records RNG-state hashes and
individual event audits, and compares TB/death/survival/person-time under declared
Monte Carlo tolerances. It is not a production ABM or a factor-membership generator.

`incidence_comparison` consumes the unchanged 2D disease-observation contract.
It returns target/predicted incidence, difference, ratio, uncertainty, source,
scenario hash and compatibility warnings. Its status is always
`validation_only_not_calibration`; fitted parameters and TBI/pressure changes are
null. Baseline incidence inversion retains its original blocked status.

All artifacts use canonical JSON and SHA-256 `profile_hash`. Scientific rows are
order-independent; event sequences retain order. New source decimals remain
unchanged; derived forward quantities retain full binary64 precision, with
explicit numerical policy. The source 3A rational quantities remain embedded.
No source bound or numerical tolerance is relabelled as a TB confidence interval.

Developer preview:

```
python -m engine.country_tb.preview_progression --output abm/output/country_tb_m3b_preview_20261029 --mc-samples 4000
```

It writes immutable hash-named JSON/CSV/Markdown with a file-checksum manifest.
Paths must stay under `abm/output/`. No production navigation, runner, service,
network download or deployment is connected.
