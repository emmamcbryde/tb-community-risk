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

Mapped components are:

- screening test: `screened * selected test cost`;
- preventive regimen: `tpt_started_total * regimen cost per started course`;
- false-positive increment: `tpt_started_false_positive * false-positive cost`;
- active TB disease care: `active_tb_cases * active-TB cost per case`;
- programme setup: intervention arm, model year 0 only;
- programme running: requires an explicit basis and allocation method;
- ADR management: optional, mapped to `tpt_adr_stop_total` only.

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

For stochastic agent-based analyses the primary aggregate ICER is:

```text
mean paired incremental cost / mean paired DALYs averted
```

Replicate ICERs are retained as diagnostic outputs only. Percentile intervals
are labelled as simulation distributions across finite-population replicates,
not confidence intervals or PSA uncertainty.

Net monetary benefit is calculated only when an exact reviewed threshold value,
currency, reference year and source are supplied. The proportion with positive
NMB is a finite-replicate, fixed-parameter simulation summary, not a PSA
probability of cost-effectiveness.

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
