# Milestone 3B forward natural history

This isolated model is for planning and sequencing, not denying care. It projects
the first active-TB episode after a supplied 2024 baseline. It does not calculate
active-TB prevalence, interventions, economics or transmission, and it does not
infer TBI from disease incidence. Read the [evidence decisions](milestone3b_evidence.md)
before interpreting any numerical output.

## Retained 2C–3A semantics

All retained implementation, schema, test and developer-preview files remain
unchanged. Infection is a repeated-event Poisson process; the most recent event
defines recent versus remote-only infection. Prior-remote-plus-recent remains a
subset, never another mutually exclusive state. Source modes, user-TBI inputs,
one-scale calibration, subgroup count weighting and unresolved 100+ populations
remain intact. Baseline TBI/risk marginals and highest-RR allocation are retained.
Adult effects are not automatically assigned to children.

`active_tb_incidence_inference_v1` still returns
`unavailable_pending_progression_and_risk_forward_operator`. A forward operator
existing is not approval of an inverse observation/calibration model.

## Curves and equations

The mathematical object supplies instantaneous hazard `h0(tau; ai, ac, sex)` and
cumulative hazard `H0`. The explicit age policy selects the current attained-age
curve, retaining age at most recent infection and the full age trajectory.
Curve groups are 0–4, 5–14 and 15+. Finer supplied ages are retained. Sex-specific
curves require explicit reviewed inputs; otherwise females and males share a
baseline curve while their counts, exposures and mortality remain separate.

For a reviewed **net** cumulative disease risk, without competing death:

```
H(t) = -log(1-F(t))
segment hazard = [H(t_next)-H(t_previous)]/(t_next-t_previous)
DeltaH = integral h0(tau+u; age_at_event, age_current+u, sex) du
P(TB in interval | disease-free at its start) = 1-exp(-m*DeltaH)
```

Age-group and hazard-knot boundaries are split explicitly. Forward conditional
increments can be calculated for an adult now without inventing a missing
historical pediatric hazard: this does not resolve the historical depletion issue.
`cumulative_hazard()` over a lifetime crossing an unavailable child curve fails.

The four-knot diagnostic surrogate has the following transformed hazards. The
unverified 5.0% two-year knot is not inserted. The interpolation is linear in
cumulative hazard, not in probability; it is not a reproduction of an original
continuous posterior curve.

| Time interval (years since infection) | Constant baseline hazard/year |
| --- | ---: |
| 0–1 | 0.0387408283164305622 |
| 1–5 | 0.00738450310921596722 |
| 5–10 | 0.00128894108852840934 |
| 10–25 | 0.000504779768726248618 |

All source risks, uncertainty bounds and high-precision transformed knots/segments
remain in the contract. Four-point interpolation must not be mistaken for evidence
at every intervening year. The source risk-to-hazard conversion is a diagnostic
surrogate pending cause-specific-hazard review. Observed cumulative incidence
with competing deaths is not silently accepted as a cause-specific hazard.

Central tail: no extrapolation. Any required time beyond 25 years leaves that
cohort unresolved. Sensitivities explicitly select the final segment hazard,
a declared fixed remote hazard, or an exponential decline toward a declared
plateau. For the latter, after knot K:

```
h(t) = plateau + (h_final-plateau)*exp(-k*(t-K))
H(b)-H(a) = plateau*(b-a)
             + (h_final-plateau)*exp(-k*(a-K))*(1-exp(-k*(b-a)))/k
```

No tail parameters are production defaults. Child curves remain external-required.
Missing curves, infection-time distributions, risk effects or mortality coverage
remain visible and block a complete mathematical result for the affected horizon.

## RR application and extreme risk

The central effect policy is `highest_applicable_reviewed_rr_v1`, interpreted via
`provisional_rr_as_hazard_multiplier_v1`. Multiple factors are neither added nor
multiplied. RR is not relabelled as an observed hazard ratio. `m*F(t)` is never
used. The approximation bounds probabilities but may not reproduce the original
effect measure, and requires scientific review and sensitivity analysis.

The 3A adapter aggregates over factor-count categories only where the same highest
RR and any-risk status imply identical forward hazards. Highest-contributor
probabilities use the original exact rational ratios. All ties are retained as
nonexclusive attribution weights; summing tied factor counts can exceed unique
cases. That sum is never added to the total. Full factor provenance stays in the
source allocation; the Monte Carlo module validates event times, not a full new
individual risk-factor population generator.

Each effective RR has a declared follow-up age scope. The adapter conservatively
intersects the included positive factors' reviewed scopes. Crossing a finite
boundary leaves that horizon unresolved; it does not retain an inapplicable RR
or switch to a legacy effect. Pediatric effects need explicit reviewed scopes.
Factor membership/effects remain fixed only within that declared scope.

No arbitrary RR cap is imposed. The result reports case counts by effective RR
and the top 1%, 5%, 10% of the resolved disease-free population ranked by baseline
expected individual risk. A declared diagnostic warning triggers when the top 1%
accounts for more than half of expected cases. Rapid saturation under extreme
candidate multipliers is a warning about scientific transport, not evidence that
those probabilities are clinically plausible. No capped sensitivity is adopted
without scientific justification.

## Competing death and the absorbing risk set

With progression hazard b(t) and external background mortality d(t):

```
Y(t) = exp(-integral [b(u)+d(u)] du)
TB_CIF(T) = integral Y(u)*b(u) du
Death_CIF(T) = integral Y(u)*d(u) du
TB_CIF + Death_CIF + Y = 1
person-years(T) = population * integral Y(u) du
incidence rate = 100000 * expected cases / disease-free person-years
```

Death and the first TB event remove future risk. Competing deaths are deaths
**before** that first TB event; this is not a model of death after TB. Baseline
active TB is excluded explicitly from each risk set and reported separately.
Assuming zero baseline active TB is permitted only as a named diagnostic policy.

Externally supplied one-year survival becomes `d=-log(survival)`, with explicit
age/sex coverage and optional subgroup-specific rows taking precedence over
general rows. Hazards are stationary across follow-up under this contract.
Reviewed national life tables require licensing and snapshot provenance. No crude
life expectancy or Australian fallback is used. No-mortality is diagnostic only.

All horizons describe the same absorbing process from baseline. Cumulative cases
are never summed across horizons. Incidence per person-year uses the integrated
risk-set denominator; 20-year cumulative cases are not annual incidence. Coverage
can shrink at longer horizons because a curve/tail becomes unavailable. In that
case covered totals describe different populations and may decrease; use the
reported common-horizon cohort population/cases for longitudinal comparisons.

The catalytic baseline does not remove pre-2024 TB episodes, treated previous TB
can remain TBI positive, prior-TB history is unknown and infection-selective
historical mortality is absent. The operator is first **after baseline**, not a
validated first-ever disease model. A synthetic historical-removal sensitivity
quantifies this gap without altering the source baseline or proposing a country
correction. Recurrence, relapse, failure, imported disease, migration and
endogenous transmission are excluded.

## Prospective infections and clock policies

Baseline-cohort mode prohibits nonzero future pressure. Exogenous-pressure mode
requires a reviewed explicit scenario: zero, constant at an explicitly supplied
2024 rate, user intervals, or bounded annual staircase increase/decline. It
cannot extrapolate the WHO shape automatically and never feeds predicted disease
back into infection pressure.

Central policy `most_recent_reinfection_resets_progression_clock_v1` resets the
clock to zero for every admitted event, without multiplying RR. First infection
and reinfection both remain possible. Event counts and cases following first
infection versus reinfection are separate. Prior remote baseline exposure is
retained in source audit fields and never adds another person or case.

No-reset and partial protection against infection events are explicitly selected
diagnostic sensitivities. Actual last-event time still defines recency under
no-reset; the effective progression clock is reported separately. Partial reset
is provided as a **single-reinfection calculation-only** function; it is rejected
as a population operator policy. The synthetic susceptibility value is not an
estimate derived from a disease IRR. Total TBI alone identifies none of these
clock/protection alternatives independently.

## Deterministic integration and stochastic validation

Supplied discrete age/time supports are used unchanged. For a catalytic source,
the adapter reconstructs the joint attained-age/last-event density from the
original pressure history, integrating it by positive eight-point Gauss-Legendre
quadrature. The original within-single-year age policy is respected; no display
midpoint or mean infection time replaces a distribution. Conditional quadrature
weights are normalized within the original 3A state/RR totals. Direct total or
recency profiles alone do not identify infection times: they require separately
reviewed conditional support. Open tails need reviewed age support.

The baseline-only solver uses analytic piecewise competing-hazard kernels;
declining tails use composite positive quadrature with an explicit work bound.
Future pressure uses Strang splitting: progression/death for half a step,
Poisson events/reset, then progression/death for half a step. Zero/one/multiple
infection-event probabilities within the infection step are retained. Event
locations on this deterministic grid are a numerical approximation, not new
source uncertainty. Explicit step refinement is required before interpretation.

The separate validator samples continuous event times by integrated competing
hazard inversion, including infection/reinfection and death; it does not round
events to the deterministic grid. It requires a caller RNG, records its state
hash, retains individual event audits and stops at the first TB event or death.
Tests use predeclared five-Monte-Carlo-SE tolerances plus a numerical-error
allowance; disease-free person-time is also checked. No global RNG is used.

Forward calculations retain full binary64 precision with stable `expm1` kernels
and `fsum` totals; source quantities and transformed curve contracts retain their
original decimal precision. Numerical tolerances and quadrature/grid error are
explicit. They are not confidence intervals. Parameter, source, scenario and
stochastic uncertainty remain separate. There is no production UI.

## Validation-only incidence comparison

The comparison retains a supplied target and its original uncertainty/provenance,
predicted rate, difference, ratio, unresolved population and exclusions. Units and
prediction measure are explicit. Instantaneous baseline rate differs from a
calendar-year average, and cohort person-time differs from a whole-population
surveillance denominator. WHO includes new/recurrent disease and additional
processes excluded here. Notifications and prevalent/screen-detected TB remain
distinct observations. The status is always `validation_only_not_calibration`;
there are no fitted parameters or changes to TBI/pressure.

See [contracts](milestone3b_contracts.md), [diagnostics](milestone3b_examples.md),
[validation](milestone3b_validation.md) and [decisions](decision_register_3b.md).
