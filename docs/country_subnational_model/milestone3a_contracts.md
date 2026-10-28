# Milestone 3A executable contracts

The new [offline JSON Schema](../../schemas/country_tb/v1/risk_allocation.schema.json)
uses Draft 2020-12, contract version 1.0.0 and implementation
`country_tb_risk_allocation_1.0.0`. Frozen 2A–2D files are not revised.

| Required family | Pure implementation |
| --- | --- |
| `country_tb.progression_risk_factor_catalog.v1` | `risk_factor_catalog.progression_risk_factor_catalog` |
| `country_tb.risk_factor_prevalence_profile.v1` | `risk_prevalence.risk_factor_prevalence_profile` |
| `country_tb.risk_factor_definition_mapping.v1` | `risk_prevalence.definition_mapping` |
| `country_tb.risk_factor_overlap_policy.v1` | `risk_policies.overlap_policy` |
| `country_tb.tbi_any_risk_association.v1` | `risk_policies.association_policy` |
| `country_tb.multimorbidity_progression_policy.v1` | `risk_policies.multimorbidity_policy` |
| `country_tb.baseline_population_risk_allocation.v1` | `risk_allocation.allocate_baseline_risk` |

Supporting families are `country_tb.progression_risk_factor_record.v1`,
`country_tb.population_risk_allocation_mixture.v1`,
`country_tb.risk_allocation_diagnostic_set.v1` and
`country_tb.risk_allocation_export_receipt.v1`.

Modules are under `engine/country_tb/`. Pure constructors validate reviewed
metadata, definition mapping, population alignment, non-overlap, exact pattern
marginals, feasible association, source identity and highest-RR isolation.
Schema validation alone is not evidence approval or a substitute for semantic
validation. All new schema references are local; no network resolution is used.

Every artifact uses canonical JSON and SHA-256 `profile_hash`. Expected numeric
quantities use `{value, exact: {numerator, denominator}}`; `value` is the decimal
display and the integer-string rational pair is authoritative. The decoder checks
their consistency. Phi is a deterministic decimal summary. Original source
receipts, mapping evidence and review status are retained. Receipt order does not
alter scientific identity; canonical scientific rows do.

`allocate_baseline_risk` requires an explicit resolved 2D baseline, catalogue,
prevalence, overlap, association and multimorbidity policy. The result embeds
these contracts and the unchanged baseline hash. No source-mode fallback,
unknown-factor fallback, implicit definition mapping or allocation policy is
accepted. The global association value must be supplied explicitly.

The result separates `complete_result`, component coverage, unresolved population,
`diagnostic_only`, and `ready_as_input_for_future_operator`. A complete diagnostic
stratification is not a reviewed progression input. `forward_operator` is always
null. Calculation-only sensitivities cannot enter population allocation.

`combine_risk_allocations` consumes a validated 2D mixture and exactly one
allocation per subgroup, checking baseline hashes and total population. Different
factor scopes remain explicit. Missing outputs, duplicate IDs or mismatched sources
fail. Incomplete groups make the combined complete result unavailable.

Developer-only diagnostics run with:

```
python -m engine.country_tb.preview_risk_allocation --output abm/output/country_tb_m3a_preview_20261028
```

Exports are immutable, scientific-hash-named JSON/CSV/Markdown bundles with a
SHA-256 file manifest. Paths must remain under `abm/output/`. There is no
production UI, runner registration, external download or network access.
