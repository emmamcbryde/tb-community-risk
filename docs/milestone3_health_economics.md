# Milestone 3 APY Event-Ledger Health Economics

Milestone 3 introduces the versioned `ltbi_health_economics_results_v2`
contract. The authoritative epidemiological input is the APY event ledger
(`ltbi_screening_event_ledger_v3`), not median summary rows or workbook
formulas.

## Scope

The current static expected-value and individual-based APY analyses estimate
direct benefits, harms and costs for screened and treated individuals.
Transmission-mediated benefits are not included. The engine does not implement
PSA, early stopping rules, ICER-based clinical recommendations or an
Australian funding threshold.

## Cost Mapping

All economic totals use authoritative `costItems` only when
`conversionStatus = valid`, `convertedTargetYearCost` is present and the cost
item target currency and price year match the analysis target. Raw mixed-year
source costs cannot enter totals.

Each cost item must also carry an explicit resource-use basis that matches the
event mapping. Examples are `per_person_screened` for IGRA/TST,
`per_started_course` for preventive regimens, `per_active_tb_case` for active
TB disease care, `total_once_at_program_start` for setup costs and
`per_adr_stop` for ADR management. A conflicting or missing basis makes the
component unavailable rather than applying a cost to the wrong denominator.

Mapped components are:

- screening test: `screened * selected test cost`;
- preventive regimen: `tpt_started_total * regimen cost per started course`;
- false-positive increment: `tpt_started_false_positive * false-positive cost`;
- active TB disease care: `active_tb_cases * active-TB cost per case`;
- programme setup: intervention arm, model year 0 only;
- programme running: requires an explicit basis and allocation method;
- ADR management: optional, mapped to `tpt_adr_stop_total` only.

Missing unit costs block totals only when the corresponding event quantity is
positive in an included annual row. Zero resource use does not require a unit
cost for that resource.

Cost conversion to the common target price year remains separate from future
discounting. Future model costs stay in constant target-year prices.

## DALY Method

DALYs use explicit reviewed assumptions for active-TB disability weight,
active-TB duration, TB case-fatality risk and years of life lost per TB death.
The current method is a scalar expected DALY-per-active-TB-case calculation,
not an age-specific life-table calculation.

For each arm and year:

```text
expected TB deaths = active_tb_cases * TB case-fatality risk
YLD = active_tb_cases * active-TB disability weight * active-TB duration
YLL = expected TB deaths * years of life lost per TB death
active-TB DALYs = YLD + YLL
```

TPT and ADR-related DALY losses require reviewed inclusion assumptions, or a
reviewed exclusion. False-positive treatment may contribute harms and costs but
does not produce infection cure or active-TB prevention.

Numerical DALY parameters require `configured_reviewed` or
`model_derived_reviewed` status. `reviewed_exclusion` is valid only for explicit
exclusion decisions such as omitting TPT health loss, ADR health loss or
post-TB sequelae from the acute primary analysis.

## Discounting

The engine calculates 3% and 0% annual discounting from the same undiscounted
annual economic ledger. With the current annual event timing convention:

```text
discount factor in model year t = 1 / (1 + rate)^t
```

Model year 0 is undiscounted. Event quantities such as people screened and TB
cases are not discounted.

## Incremental Results

Comparator and intervention costs and DALYs are calculated by arm, model year
and paired replicate. Incremental cost is intervention cost minus comparator
cost. DALYs averted are comparator DALYs minus intervention DALYs.

Annual rows include `costComplete`, `dalyComplete`, missing-component lists,
`includedInEconomicAnalysis` and `economicExclusionReason`. Replicate rows
include `costPairComplete`, `dalyPairComplete` and `economicPairComplete`.
Incomplete pairs have null incremental cost, DALYs averted, ICER and NMB.

For stochastic agent-based analyses the primary aggregate ICER is:

```text
mean paired incremental cost / mean paired DALYs averted
```

Replicate ICERs are retained as diagnostic outputs only. Percentile intervals
are labelled as simulation distributions across finite-population replicates,
not confidence intervals or PSA uncertainty.

A numerical primary ICER is shown only in the increased-cost/increased-health
gain quadrant. Other quadrants are classified as dominant, dominated, cost
saving with health loss, health gain at no additional cost, cost saving with no
material health difference, increased cost with no material health difference,
no material difference, or incomplete.

Net monetary benefit is calculated only when an exact reviewed threshold value,
currency, reference year and source are supplied. The proportion with positive
NMB is a finite-replicate, fixed-parameter simulation summary, not a PSA
probability of cost-effectiveness.

The threshold currency and reference year must match the analysis target
currency and target price year unless a future explicit threshold-normalisation
record is added. Deterministic expected-value analyses report deterministic NMB
when available; they do not report probability-positive-NMB.

## Economic Horizon

`economicHorizonYears` defaults to the epidemiological follow-up horizon. Annual
event rows outside the economic horizon are retained for audit with their event
year and quantities, but are excluded from economic totals and carry an
exclusion reason.

## Provisional Results

Economic outputs inherit unresolved or provisional epidemiological, cost, DALY
and threshold status. In particular, unresolved APY recent-LTBI assumptions,
development compatibility mode, unresolved cost conversions, unresolved DALY
inputs and unresolved GDP-per-capita thresholds prevent clinician-ready
cost-effectiveness conclusions.

Legacy `costs.*`, `unitCosts` and `costEffectiveness` fields remain as
compatibility mirrors where older callers require them. They are derived from
the event-ledger calculation and do not override authoritative annual or
replicate economic results.
