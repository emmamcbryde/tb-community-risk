# Milestone 2D executable contracts and APIs

Revision `1.0.0`, implementation `country_tb_baseline_composition_1.0.0`.
[Offline machine schema](../../schemas/country_tb/v1/baseline_source_composition.schema.json)
uses explicit schema-ID dispatch and contains no external-reference download.
Pure APIs additionally validate numerical, membership, coverage and source
invariants. Frozen 2A contracts and implemented 2C calibration contracts remain
unchanged; no implicit migration occurs.

| Contract family | Main API |
| --- | --- |
| `country_tb.baseline_tbi_source.v1` | `baseline_tbi_source.baseline_tbi_source` |
| `country_tb.baseline_tbi_resolution.v1` | `baseline_tbi_source.resolve_baseline_tbi` |
| `country_tb.population_subgroup.v1` | `population_subgroups.population_subgroup / validate_subgroup` |
| `country_tb.population_mixture.v1` | `population_subgroups.population_mixture` |
| `country_tb.subgroup_demography.v1` | `composition_common.population_demography / validate_demography` |
| `country_tb.active_tb_inference_targets.v1` | `active_tb_targets.active_tb_targets / validate_active_tb_targets` |
| `country_tb.direct_tbi_profile.v1` | `direct_tbi_profile.direct_tbi_profile / allocate_direct_profile` |
| `country_tb.baseline_tbi_interface_specification.v1` | `baseline_tbi_source.interface_specification` |
| `country_tb.baseline_tbi_synthetic_scenario_set.v1` | Isolated developer preview |
| `country_tb.baseline_tbi_export_receipt.v1` | Immutable preview export file hashes |

The source's scalar primary `mode` is mandatory. Incidence configuration
requires targets plus explicit nullable pressure-history/operator fields;
operator must remain null in 2D. Observed/user configuration requires exactly
one submode and its compatible inputs:

| Submode | Required scientific input |
| --- | --- |
| `overall_prevalence_anchor_v1` | Overall anchor, fixed pressure history, one-scale root policy, time/window |
| `age_specific_prevalence_anchors_v1` | Existing reviewed 2C observations and one-scale policy, fixed history, time/window |
| `direct_age_specific_profile_v1` | Direct profile, explicit recency policy and nullable comparison profile, time/window |
| `direct_recent_remote_profile_v1` | Direct profile with all exclusive states, recency policy and nullable comparison, time/window |

Source/resolution records include geography and level, subgroup ID, 2024,
source/review/provenance/hash, empirical/user-defined/model-inferred/diagnostic
classification, recency classification, status and downstream usability.
Resolution embeds unchanged inputs, calibration diagnostics, comparison
profile, coverage gaps, component counts and missing approved dependencies.
Unknown or mixed modes/configurations fail; there is no default source or
fallback from a failed calibration to direct values or incidence inference.

Subgroup demography may come from the unchanged frozen 2A API or explicit
user cells. A new-user demographic adapter supplies those same counts to the
unchanged 2C API when calibration is selected. Population and survey-sample
denominators are distinct. Empty cells are allowed, but each subgroup's total
must be positive and finite. Coarse cells are not split implicitly. Reviewed
wider-cell catalytic averaging and terminal-age weights use the existing 2C
policies. Direct profiles specify their own explicit within-band value policy.

Subgroups preserve current residence separately from all optional
birth/previous residence/exposure/arrival metadata. Every optional absence is
null. Their risk-factor field is exactly `{status: deferred, profile: null}`.
Membership declarations require common partition ID, unique category IDs,
definition, source/review/provenance and nullable complete person-ID lists.
Mixtures validate reviewed partition evidence, unique subgroup IDs, compatible
geography/year/grids, denominators, weights and optional reference age/sex
counts. See [scientific behavior](milestone2d_sources_subgroups.md) for the
partial-component rules and limitations of aggregate membership evidence.

Scientific hashes use the established canonical UTF-8 sorted-key JSON and
SHA-256. Known unordered rows, subgroup lists, membership IDs and partition
IDs are canonicalized deterministically. Original direct/demographic/active-
TB records retain independent receipt checksums and are checked against their
normalized artifacts. Scientific hashes exclude established display, receipt
and profile-hash fields, so reordering source rows changes the raw receipt
without changing scientific identity. Model comparisons and source histories
are embedded; no latest catalogue lookup or network is needed to resolve the
new contracts.

Developer command:

```
python scripts/preview_country_tb_sources.py --output abm/output/country_tb_m2d_preview
```

It produces invented six-scenario results, the pure interface specification,
summary CSV, README and file-checksum manifest in a new immutable hash-identified
folder. Alternative output paths outside `abm/output/` fail. No Streamlit,
APY, explicit-TBI, Starsim, MATLAB, progression or economics code is imported.
