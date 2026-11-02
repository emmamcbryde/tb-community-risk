# Milestone 3E review decisions

Date 2026-11-02. These are software/interface decisions, not scientific approval
of country parameters or external validation.

## Implemented for internal review

- Separate Streamlit entry, explicit seven-view navigation, no root legacy pages.
- Frozen 2024 source selection, source/estimate/projection/coverage captions and
  retained female/male data. Custom/subnational proxies require explicit choice.
- Operational user-TBI route: total-only, catalytic derivation by explicit opt-in,
  direct age/sex/states, CSV, subgroup sources and population-weighted mixtures.
- Visible missingness, open ages, unresolved subgroups and source-preserving values.
- Definition-aware risk prevalence, highest reviewed RR centrally, legacy effects
  only by explicit sensitivity selection, no adult-to-child effect fallback.
- Normalized within-stratum association and achieved phi with feasible bounds;
  no changes to TBI, risk, recency or population marginals.
- Deterministic instantaneous sensitivity and optional bounded one-year 3B
  forward projection. Unknown clocks/curves/effects remain unresolved, not zero.
- Material-input hashing, stale-result/export protection, lazy full-precision
  JSON/CSV, and source-bound provenance/version records.
- 15-item scientific gate plus separate native numerical clearance for a guarded
  uploaded diagnostic configuration. No UI checkbox or direct route bypass.

## Not approved

- Production incidence-to-TBI inference, automatic WHO fitting or fictional
  counterfactuals presented as empirical country estimates.
- The four-point adult surrogate as a central national curve, pediatric scaling,
  an unreviewed remote tail, legacy RR 74 as a universal hazard, or an arbitrary cap.
- National data automatically representing a subgroup, Australian mortality for
  another country, unlogged assay/migration/recurrence adjustments, or use of an
  old cached result after material edits.
- Scientific sign-off inferred from AppTest, unit tests, an inverse root or an
  empirical residual. 3D remains internal synthesis/model checking.

## Reviewer questions and remaining decisions

Confirm usability and labels for source/model/sensitivity/unresolved quantities;
verify CSV definitions, age/sex coverage and subgroup partitions; assess whether
the explicit catalytic/hybrid assumptions are sufficiently visible; review
factor definitions and conditional effect measures; decide adult, pediatric,
remote-tail, mortality and assay inputs; reconcile resident/first-episode targets
with migration, importation, prior TB and recurrence; approve appropriate exclusion
and numerical thresholds only with source-bound evidence.

Cold WPP validation is a known first-load cost, explicitly invoked and cached.
The optional prospective calculation has a preflight node budget and is not
eagerly expanded during page rendering. These are software limits, not scientific
acceptance thresholds.

## Deferred

External scientific review, matched data acquisition and held-out validation;
production UI and country inversion approval; stochastic UI; interventions,
screening/preventive-treatment strategies, economics, DALYs/QALYs/ICERs, and
endogenous transmission. No Starsim/Summer or MATLAB work is started.

Recommended next milestone: reviewer recheck of this isolated interface and its
assumptions, then targeted scientific input decisions/data acquisition. A dynamic
model remains necessary where incidence and infection pressure are endogenous.
