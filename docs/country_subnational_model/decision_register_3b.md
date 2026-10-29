# Milestone 3B decision register

## Approved mathematical and software policies

- Time since most recent infection drives central progression. Reinfection resets
  the central clock without multiplying RR or adding another person.
- Highest applicable RR centrally, with provisional RR-as-hazard interpretation
  explicitly named; cumulative-hazard probabilities, never `m*F(t)`.
- First active-TB episode after baseline is absorbing. Death competes explicitly;
  baseline active TB is excluded from the ordinary future-progression risk set.
- External age/sex/subgroup mortality, or labelled no-mortality diagnostics; no
  life-expectancy substitute or cross-country mortality fallback.
- Exogenous future infection scenarios are explicit, separate from retrospective
  WHO pressure, and have no feedback from generated infectious cases.
- Source TBI/risk marginals are preserved at initialization. Forward depletion
  changes the surviving population naturally; it does not recalibrate baseline.
- Actual recency and sensitivity progression clocks are separate. Prior remote
  audit remains a subset, not another mutually exclusive state or case.
- Positive quadrature, analytic competing kernels, explicit renewal-grid
  refinement and a continuous event-time Monte Carlo validator with caller RNG.
- Incidence comparisons are `validation_only_not_calibration`. Baseline TBI
  inversion remains `unavailable_pending_progression_and_risk_forward_operator`.

## Requires scientific review

- Adult curve: the requested unqualified five-knot candidate is rejected; the
  two-year knot was not verified. Four verified points are diagnostic surrogates
  pending cause-specific-hazard, population and age transport review.
- Pediatric 0–4 and 5–14 curves remain external-required. No adult or invented
  production child curve fills the gap.
- Long tail: no extrapolation centrally; final-slope, declining and fixed-remote
  alternatives remain explicit sensitivities with no automatic approval.
- RR-as-hazard approximation, contemporary effects, treatment-modified effects,
  extreme multipliers, high-risk dominance and age/sex applicability over follow-up.
- Full-reset biology versus partial reset, no reset and reinfection protection.
  Disease IRRs and total TBI alone do not identify these mechanisms separately.
- Mortality sources, licensing and immutable snapshots; period survival transport
  into future cohort survival, risk-specific mortality and cause-of-death scope.
- Prior-TB history, historical TB depletion, selective mortality, recurrence and
  the relation between first-post-baseline events and WHO new/recurrent episodes.
- Future infection pressure and migration/imported disease assumptions.

Partial reset is single-event calculation-only, not a population default.
No-reset/protection population sensitivities require diagnostic operation. No
arbitrary cap is justified or applied. No diagnostic mathematical result is
labelled a complete country estimate.

## Deferred

- Inverse active-incidence calibration and approval of incidence-based TBI inference.
- Production UI, deployed applications and production runner integration.
- Screening/treatment/interventions, economics, DALYs/QALYs and endogenous transmission.
- Historical migration, recurrence/treatment states and an individual factor
  membership population generator.

The recommended following milestone is **independent reconciliation of predicted
versus observed incidence, identifiability analysis, and explicit approval or
rejection of incidence-based TBI inversion**. It must assess observation processes,
coverage and source assumptions independently of this implementation. It has not
been started automatically.
