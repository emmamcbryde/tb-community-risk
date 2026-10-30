# Diagnostic inverse contract and identifiability

Runtime is pure/offline. Production availability is always false. Existing
public source dispatch and2A–3B contracts remain unchanged.

## Mapping

For fixed q, exact age a and baseline2024, define
`Q(tau)=integral(2024-tau,2024)q(t)dt`. Last-event time-since-infection has
unconditional density `alpha*q(2024-tau)*exp[-alpha*Q(tau)]`,0<=tau<=a, plus
never-infected mass `exp[-alpha*Q(a)]`. Repeated events remain possible. Recent,
remote-only and prior-remote-plus-recent retain the2C formulas; the last is an
audit subset, never added again.

Response kernels use the unchanged3B instantaneous-hazard/trajectory functions.
Exact exponential interval masses integrate piecewise-linear responses in tau.
Conditional interval means use stable small/large-hazard formulas. Bins split
at pressure/hazard knots, recent-window boundaries and knots minus follow-up.
Kernel-width refinement is independent of the alpha search. Only exact supplied
ages or reviewed within-cell supports enter; no demographic band midpoint or
invented terminal age. Unresolved100+ stays in denominators.

Instantaneous incidence is `100000*sum N*m*h0(tau)*mass / sum N`. Mortality and
future pressure have **exactly zero effect on this instantaneous measure**.
Prospective incidence uses first-episode cases divided by disease-free
person-time, including depletion/competing death. Observation measure and period
must match. Future exogenous pressure is fixed separately from historical alpha;
WHO q is never extrapolated automatically. Central reinfection resets the clock.

Risk factors are conditionally independent within exact age/sex/subgroup, with
at most11 factors/2^11 patterns. Only the highest RR applies. Scalar scans can
compress patterns exactly by maximum RR and union membership; detailed factor
audits retain full patterns. Frechet union association preserves the conditional
factor pattern distribution, all marginals and recent/remote composition.
Only alpha fits; effects, mortality, prevalence and other policies remain fixed.

## Structural identification

The bounded nonnegative bracket is scanned with dense linear, independent
half-grid and log-spaced near-zero points. Finite-difference sign changes refine
stationary points. All detected crossing, boundary and tangent roots are
reported. No monotonicity is assumed. Grid-resolved identification is not a
global analytical proof against arbitrarily narrow unsampled oscillations.

An actual non-monotonic operator fixture has q=1,age20,h(tau)=0.05 at5–10years,
zero elsewhere. `I(alpha)=5000*(exp(-5alpha)-exp(-10alpha))`, maximum1250 at
`log(2)/5`. Target800 gives roots `-log(.8)/5` and `-log(.2)/5`, approximately
0.0446287103 and0.3218875825. Both are reported and neither adopted. Any
non-monotonic search interval fails single-target identification. Flat/saturated
maps, unresolved kernels, no root, infeasible target and failed convergence also
fail without clipping.

Root finding uses deterministic bounded bisection. Explicit numerical tolerances:
incidence1e-7,alpha1e-10,flat derivative1e-7. Derivative step and half-step,
inverse sensitivity `1/I'`, relative condition `abs[I/(alpha*I')]`, residual,
boundary and attainable-range status are reported. Relative condition at0 is
unavailable. Independent narrower-bracket recovery and halved kernel width check
convergence; future-pressure time steps are also halved. An excessive refinement
budget or undefined refined kernel returns a structured refusal. Integration
allowance is1e-5 times max(1,rate), not a scientific
uncertainty. Numerical bracket/grid evidence accompanies each fit.

## Gate and approvals

Review receipts bind to the exact scenario hash and target population. Hashes
prove identity, not scientific validity. Required checks: compatible observation;
matching measure/period/age/sex/geography; reviewed material curves and explicit
remote tail; reviewed risk application; mortality and immutable licensed snapshot
for empirical configurations; addressed migration/importation/recurrence/prior
TB/depletion/future pressure; complete provenance; approved unresolved threshold.
After fitting, identification, attainable target, concentration and numerical
checks must pass. Unreviewed risk effects fail irrespective of concentration.

**No default unresolved threshold.**0%,0.1%,1% are candidate thresholds presented
for approval. Until answered, actual configurations record pending and block.
A small allowed share cannot approve missing material children/tails. Fictional
test protocols declare0% only under `synthetic_test_only`; these cannot authorize
empirical populations or represent user/biomedical approval. Production stays
false even after a successful diagnostic fit.

This solver also requires a defined forward prediction for the target's entire
risk set. A positive exclusion ceiling cannot supply an unknown100+ hazard or
an observation model for a covered-only denominator. Such inputs return a
structured gate refusal; no threshold is silently used to drop people.

| Status | Meaning |
| --- | --- |
| blocked_scientific_gate | No root search or forward fit invoked |
| rejected_nonidentifiable_or_infeasible | No unique admissible root |
| blocked_postfit_gate | Numerical/concentration checks failed |
| diagnostic_only_not_approved_for_production | Fictional or externally reviewed diagnostic conditions satisfied; no production authorization |

Actual Viet Nam/Australia configurations fail several gates. A country forward
scan can report **unadopted crossing brackets**, never an empirical fitted scale.
A separate fictional closed population with Viet Nam-shaped pressure explores
mathematical roots at a synthetic target numerically equal to182. Its demography,
risk laws and clinical surrogate are declared mathematical assumptions, not
approval of the actual country target or national estimates.

## Practical identification and user TBI

Exports retain scenario hashes, fitted scale where permitted, age/sex TBI,
recent/remote, predicted incidence, derivative/condition and failure information.
Source bounds, parameter/scenario variation and integration error are separate;
none is labelled a TBI confidence interval. Progression, tail, association,
factor omission, backcast, source/target, mortality and future-pressure scenarios
show conditional scale instability. Missing pediatric or migration information
remains missing. Instantaneous mortality sensitivity is zero by construction;
prospective mortality/future effects are evaluated separately.

The2D overall/age-specific anchors, direct age/sex/recent/remote and subgroup
pathways remain available. Validation stores an unchanged supplied profile and
reports incidence disagreement; it never forces user TBI to match incidence.
Incomplete subgroups and100+ prevent apparently complete results unless the
necessary reviewed age/time support resolves their forward profile.

## Contracts

All six use canonical JSON, SHA-256,version1.0.0,source/review/provenance and
`production_available:false`. The [machine schema](../../schemas/country_tb/v1/incidence_inversion.schema.json)
checks shape; constructors enforce numeric/semantic rules.

| Identifier (`country_tb.` prefix and `.v1` suffix) | Content |
| --- | --- |
| active_tb_incidence_observation | Original target, person-time/period, episode/acquisition/residence, uncertainty and reviewed adjustments |
| incidence_inversion_gate | Bound scientific receipts, explicit threshold approval, failed criteria/post-fit requirements |
| incidence_to_tbi_diagnostic_calibration | One-scale fit or structured refusal, full gate and scenario/target hashes |
| incidence_identifiability_report | Grid, extrema/all roots, range, derivatives/conditions, refinement and numerical limitations |
| progression_evidence_review | Adult/child/tail/mortality decisions and11 factor qualifications, including explicit unknowns |
| incidence_reconciliation_report | Preserved user profile, disagreement and exclusions; validation only |
