# Milestone 2C: one-scale reinfection-capable catalytic calibration

This isolated Python pathway is for scientific review and planning/sequencing,
not for denying care. It is not a production country estimator. No reviewed
empirical TBI anchor exists in the current review package. All demonstrated
scales, anchors and assay-performance parameters are diagnostic assumptions.
The [recovery/validation report](milestone2c_validation.md), [contracts](milestone2c_contracts.md),
[examples](milestone2c_examples.md) and [decision register](decision_register_2b.md)
record implementation, evidence limits and remaining decisions.

## Infection events and baseline states

The central model is a non-homogeneous Poisson event process:

\[
\lambda(t,a,s,g)=\alpha q(t,a,s,g),\quad \alpha\geq0.
\]

Pressure applies equally before and after any infection. Previous infection
provides no acquisition protection. Repeated events are possible; there is no
clearance or reversion. Any past event establishes the model's `any_tbi` state.
The most recent event determines recent versus remote-only classification.
A recent reinfection overrides remote-only classification; earlier remote
exposure remains auditable. The provisional progression-clock policy resets
the clock at the most recent reinfection. It has no active-TB progression link.

For exact baseline time t0 within 2024, attained age a and window w, split
lifetime integrated pressure at max(0,a-w) into H_remote and H_recent.

```
never_infected          = exp(-(H_remote + H_recent))
recent                  = 1 - exp(-H_recent)
remote_only             = exp(-H_recent) * (1 - exp(-H_remote))
any_tbi                 = recent + remote_only
prior_remote_plus_recent = (1-exp(-H_remote)) * (1-exp(-H_recent))
expected_lifetime_infection_events = H_remote + H_recent
```

The first three states are mutually exclusive and sum to one. The
`prior_remote_plus_recent` subset is already within `recent` and is never
added to total TBI or population reconciliation. It indicates at least one
remote event and at least one recent event; it is not every possible repeated
event (multiple events can also all fall in one interval).

The old Milestone 2B diagnostic remains unchanged as historical software
validation. Its state semantics are superseded for this pathway and it is
never imported by the new calibration/profile implementation.

## Conditional most-recent-event distribution

For [c,b], conditional on at least one event:

\[
f_{last}(u)=\lambda(u)\exp[-\int_u^b\lambda(v)dv]/(1-\exp[-\int_c^b\lambda(v)dv]).
\]

`most_recent_event.py` supplies conditional CDF, inverse CDF, expected time
since the most recent event, event-calendar-time and time-since quantiles
(5%, median, 95%), and sampling requiring an explicit caller RNG. Separate
summaries are returned for any TBI, recent, and remote-only. Conditioning on
remote-only includes the requirement of no subsequent recent event.
Zero-probability groups have null summaries. Very rare groups retain a log
probability and a numerically defined conditional law even when their
unconditional probability underflows the stored numeric range.

Piecewise exponential density integrals and moments are analytic. Constant
eligibility segments use an analytic inverse; linear birth-eligibility
segments use bounded inversion with 1e-24-year tolerance. Uniform within-cell
birth eligibility is integrated explicitly. It is never represented by a band
midpoint. Reviewed discrete open-tail ages mix unconditional densities before
conditioning. Small-hazard series avoid subtractive cancellation; scaled
log masses preserve rare-group conditionals at very large hazards.

The fixed Decimal context uses precision 80 and half-even rounding. Stored
state probabilities round at 40 decimal places, with an absolute 1e-40 bound;
complement/difference construction makes the exclusive states reconcile
exactly. Thus an extremely small valid stored probability can round to zero;
this is distinct from replacing missing inputs or clipping invalid targets.
Expected event counts and conditional time distributions retain their own
precision. Caller Decimal settings and RNG state do not alter deterministic
scientific results.

## Pressure shape and history

The central descriptor is a centred, unweighted local-linear OLS fit to
positive log WHO incidence, using a fixed five-calendar-year window. Endpoint
windows truncate visibly to three observations. Smoothed levels are
exponentiated and divided by the smoothed 2024 reference level. Annual q is
constant on [year,year+1), dimensionless, and labelled
`infection_pressure_shape_proxy`. Scale alpha has units per person-year.
No force-of-infection estimate is inferred from incidence. There is no Styblo
conversion, fitted progression parameter, disease-forward operator or
incidence likelihood.

Raw WHO incidence and lower/upper bounds remain embedded unchanged. The
output preserves raw and smoothed series, level/log residuals, local slopes,
window membership, edge flags and specific 2020-2022 diagnostics. Coherent WHO
lower/central/upper trajectories are separately smoothed and normalized.
Normalizing each curve removes amplitude information: WHO upper/lower curves
need not order TBI profiles and must not be called TBI confidence intervals.

Missing, zero, negative or nonfinite selected log-inputs fail with offending
years. Absent calendar years also fail with the missing-year list. No epsilon,
interpolation or replacement is invented. A reviewed supplied shape/history
is a separate supported alternative, including explicitly supplied zero
intervals with a positive normalization reference. PRK has no usable frozen
WHO series and remains blocked for incidence-informed reconstruction. There
is no post-2024 extrapolation; the final annual interval ends at 2025, allowing
retrospective integration to a specified time within 2024.

The central pre-observation policy holds the earliest available normalized
shape constant back to the required birth time. Sensitivities support bounded
backward log-linear continuation (annual calendar-grid levels evaluated at
bin starts, with supplied slope and finite bounds) and reviewed supplied
historical intervals. Bounds reject violations instead of capping values.
Supplied history must cover the reconstructed period exactly and declare a
continuous or reviewed-discontinuous join. Histories never silently truncate
at 2000. Backcast policy, parameters, first observed year, reconstructed
period, bounds, warnings, provenance and scientific hash are embedded.

## Observations, assays and calibration

Observations retain geography and level, survey year and exact evaluation
time, age range, sex, actual numerator/denominator or null, prevalence,
uncertainty, assay/product/threshold, sampling frame and design, age/sex
weights and source, review status, observation policy and provenance hash.
Original source records and their byte-canonical checksums remain in audit
receipts. Missing values are not zero and prevalence-only observations never
acquire invented counts. Counts inconsistent with published weighted/rounded
prevalence require an explicit source-preserving rationale.

Direct latent-TBI observations require a reviewed direct operator. Binary
assays use p_positive = Se*p_TBI + (1-Sp)*(1-p_TBI). Se and Sp require supplied
provenance and review, with Se+Sp>1; they stay fixed during scale fitting.
TST, IGRA and TBST remain distinct hash-identified operators with distinct
thresholds. Real-test perfect identity is rejected. Synthetic identity is
explicitly confined to invented validation fixtures.

One non-negative scale is shared across ages and sexes. Supported methods:

- Supplied scale with provenance and finite policy bounds, without fitting.
- One overall anchor: deterministic bounded bisection, matched to its actual
  survey time and reviewed population weights.
- Multiple count observations: independent binomial likelihood using actual
  numerator/denominator and compatible simple-random/independent sampling
  design. Unit observation weights prevent reweighting counts implicitly.
- Multiple prevalence observations: explicitly selected weighted least
  squares with supplied uncertainty and declared positive weights. No counts
  are inferred. Complex survey likelihood/covariance models remain deferred.

Bounded multi-observation optimization scans linear and logarithmic scale
knots, refines detected local basins with golden-section search and compares
boundaries. This deterministic basin search is a diagnostic optimizer, not a
proof of a unique global optimum. Exact synthetic recovery is tested. Results
retain objective (or log likelihood excluding the constant binomial term),
per-observation contributions, observed/achieved values, residuals, feasible
ranges, numerical convergence, rank/derivative, boundary status and warnings.
A converged result whose residuals fail declared compatibility criteria is
`incompatible_anchors` and cannot generate an accepted model profile.

Impossible individual or aggregate targets fail without clipping, including
assay floors, structural upper limits and finite-scale bounds. A true unit
latent-TBI target is unattainable at finite pressure. Coherent uncertainty
endpoint inversion is supported for one overall anchor; multiple-observation
uncertainty uses separate reviewed scenarios rather than a moment-based
substitute for the selected objective.

Total TBI identifies one scale only conditional on pressure shape, historical
backcast, assay observation model, no-clearance assumption, reinfection-event
process, demographic weights and migration assumptions. It does not
independently identify reinfection protection, clearance, progression-clock
effects or infection-time detail. Infection-time summaries are shape-
conditional reconstructions. No assay parameter is estimated with alpha.
National anchors are not automatically transported to subnational areas.

## Ages, cross-sectional interpretation and reconciliation

Frozen 2024 demography weights people alive in 2024. Single-year cells use an
explicit uniform attained-age distribution integrated analytically, with
pressure/birth/window boundaries split exactly. A wider finite cell needs an
explicit reviewed uniform policy; no display-band midpoint enters a
calculation. This within-year distribution remains an assumption for review.
For the country diagnostic the exact convention is t0=2024.5, using annual
2024 pressure retrospectively with ABS June/WPP July stocks.

An open age cell is never treated as its lower-bound exact age. It needs a
reviewed within-cell age distribution with provenance; otherwise all its
population stays visibly unresolved/unestimated. Population and proportion
excluded from age-specific estimation are reported. Full-population point
prevalence is null when age coverage is incomplete; covered prevalence and
conservative missing-age bounds are separately labelled.

Females and males remain separate internally. The same initial shape applies
to both unless reviewed sex-specific exposure curves are supplied; population
size never changes acquisition hazard. Derived both-sex totals and separate
derived-both display bands reconcile counts. Published both-sex rows remain
reconciliation-only, including genuine source-rounding discrepancies. Five
display bands remain 0-4, 5-14, 15-44, 45-64 and 65+.

Survival is assumed independent of infection. Historical deaths are not
simulated and living demographic counts receive no second survival factor.
Infection-selective mortality is not handled. No migration histories are
invented: national and subnational residents may have acquired infection
elsewhere. This limitation is especially material for Australia and prevents
interpreting a national incidence backcast as every resident's actual exposure.

## Overrides, uncertainty and isolation

`reviewed_tbi_override` creates a new, source-preserving artifact referencing
the unchanged parent profile and reviewed latent-TBI observation artifact.
Age/sex overlaps and incompatible geography/year fail. Total-only override
rows have null recent, remote-only, prior-remote-plus-recent and infection-time
composition. They are not replaced with model-derived composition. A reviewed
observation can alternatively be selected explicitly as a new calibration
anchor; that creates a distinct conditional model profile and hash.

Scenario sets preserve source identities and parameter choices separately:
WHO bound trajectories (source-data sensitivity), historical policies
(structural scenarios), assay/anchor/scale parameters, and stochastic
validation (explicit RNG, no stochastic uncertainty in deterministic output).
No scenario spread is labelled a confidence interval. Immutable hash-identified
exports stay under `abm/output/`; saved profiles rebuild from embedded original
inputs, without consulting a latest source catalogue.

No risk-factor allocation, associations, interventions, screening/treatment,
economics, DALYs/QALYs, endogenous transmission, Starsim, APY/explicit-TBI
runner, MATLAB, deployed application or production navigation is connected.
No frozen Milestone 2A contract or artifact is changed. No merge, tag or
deployment is part of this milestone.
