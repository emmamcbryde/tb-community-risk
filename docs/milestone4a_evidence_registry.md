# Milestone 4A APY Evidence Registry

Milestone 4A adds a machine-readable APY assumption evidence registry at
`data/apy_assumption_evidence_registry.csv` and a read-only readiness gate in
`engine.apy.evidence`.

The registry separates software inputs from reviewed clinical or economic
evidence. A value found in a CSV, MATLAB default or validation fixture is not
treated as clinician-ready unless it has reviewed provenance, source-year
metadata where relevant, and a reviewed status.

## Schema

Each row uses a stable `assumptionId` and records category, description,
current value, units, currency and price-year fields for costs, source
citation/location, repository path and variable, derivation method, inflation
index, inclusion/review status, provisional flag, bundling and double-counting
metadata, notes and unresolved reason.

Review statuses accepted for clinician-ready numerical assumptions are:

- `configured_reviewed`
- `model_derived_reviewed`

Explicit exclusions use `reviewed_exclusion` and require a rationale/scope in
the notes or rationale field. Repository-only statuses such as
`unreviewed_repository_input`, `unresolved`, `legacy_placeholder` and
`migrated_legacy_unverified` remain provisional.

## Readiness Gate

`assess_apy_reference_readiness(config, economics_config, registry)` returns:

- `epidemiologyReady`
- `costReady`
- `dalyReady`
- `thresholdReady`
- `overallClinicianReady`
- unresolved/provisional assumption IDs
- evidence and bundling conflicts
- compact readiness rows for UI and workbook export

The gate is additive. It does not alter event-ledger generation or the
Milestone 3 health-economic calculations.

## Repository Audit Summary

Cost variables requested in the milestone were found in `data/costs.csv`,
`abm/build_economics_preset_kwab150_v9.m` and `engine/apy/economics.py`:

- `cscreenqft = 113.48`
- `cscreentst = 116.07`
- `ctreat3HP = 165.5072`
- `ctreat4R = 123.3172`
- `ctreat3HR = 134.2272`
- `ctreat6H = 187.7508`
- `ctreat9H = 254.8544`
- `ctb = 19079.6`

These locations prove the software input values and local row descriptions.
They do not provide reviewed source citations, explicit source price years or
inflation-index values, so the rows remain provisional.

Conflicting or unmapped repository epidemiological values were also recorded:

- APY config defaults use IGRA sensitivity/specificity `0.95`/`0.98`, while
  `data/parameters.csv` contains `snqftgit = 0.7` and `spqftgit = 0.8553`.
- APY config defaults use TST sensitivity `0.80`, BCG specificity `0.55` and
  non-BCG specificity `0.97`, while `data/parameters.csv` contains
  `sntst10`, `sntst15`, `sptst10` and `sptst15`.
- The regimen library and `data/parameters.csv` contain different completion
  and efficacy values for some regimens.
- The synthetic `0.35` recent-LTBI validation fixture remains excluded from
  APY reference evidence.

The APY baseline recent-LTBI proportion remains unresolved. The Markov
recent-to-remote transition rate of `0.2/year` is recorded separately as
inherited model structure and does not validate the baseline recent fraction.

## Cost and Bundling Rules

Clinician-ready costs require numerical original value, currency, explicit or
reviewed-inferred source price year, source citation, compatible resource-use
basis, reviewed status, valid target-year conversion and no unresolved currency
conversion. Same-year costs still require provenance and reviewed status.

Foreign-currency cost inputs cannot enter economic totals without an explicit
currency-conversion record.

Essential costs for selected screening tests, preventive regimens and
active-TB care cannot disappear through a generic exclusion. Omission requires
a reviewed outside-perspective exclusion or explicit reviewed bundling into a
valid destination cost item.

## Workbook and UI

Workbook exports now include:

- `Assumption_evidence_registry`
- `APY_readiness`
- `Evidence_conflicts`
- `Bundling_conflicts`

The Economics page displays a compact read-only readiness section by category
and the main unresolved assumptions.
