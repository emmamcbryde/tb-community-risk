# Milestone 3B diagnostic results

**Invented TBI, exposure and conditional age/time support. Not Australian or
Viet Nam incidence estimates.** Adult progression uses a source-audited diagnostic
surrogate, not an approved country curve. Except where explicitly identified,
future pressure and mortality are zero. Presentation values below are rounded;
the export retains full calculated precision.

## Constant-hazard analytic check

With baseline hazard 0.02/year and no mortality, risk is `1-exp(-0.02*T)`.
The deterministic solver reproduces 1/5/10/20-year risks of 1.98013%, 9.51626%,
18.12692% and 32.967995%. Risk-set depletion leaves disease-free probability
`exp(-0.02*T)`; the incidence rate per disease-free person-year remains 0.02.
Twenty-year cumulative risk is not an annual rate.

## Recent, remote and extreme multipliers

Under the diagnostic four-knot curve, an infection 0.1 years before baseline
has much greater one-year progression risk than an infection ten years before
baseline. For infection exactly at baseline, hazard-multiplier risks are:

| Effective multiplier | One-year risk |
| ---: | ---: |
| 1 | 3.8% |
| 3.6 | 13.0178257% |
| 25.3 | 62.4742550% |
| 74 | 94.3120558% |

These bounded probabilities demonstrate rapid saturation and the scientific
limitations of transporting candidate RRs to hazards. They are not validated
clinical predictions. No cap is imposed. A separate scenario puts five infected
RR-74 people in a population of 1,000, with everyone else uninfected and no future
pressure. Those 0.5% generate all expected cases. The top-1%, top-5% and top-10%
case fractions and the small-group dominance warning are reported.

## Reinfection and tail sensitivities

For 1,000 initially remote-infected adults, tau=10, RR=1 and explicitly supplied
exogenous infection pressure 0.05/year:

| Policy | Five-year expected cases | Twenty-year expected cases |
| --- | ---: | ---: |
| Central full reset | 13.045136 | 64.425929 |
| No-reset sensitivity | 2.520716 | 9.426981 |
| Synthetic susceptibility 0.2 sensitivity | 4.707831 | 21.369046 |

These scenarios use an explicit fixed-remote tail sensitivity of 0.00038/year.
The susceptibility input is invented, not an estimate of infection protection
from a disease IRR. Actual last-event recency remains distinct from the no-reset
progression clock. A partial reset retaining half the old clock is demonstrated
only as a conditional single-event calculation, not a population default.

All four tail policies are exported. Under no extrapolation, a tau-10 stratum
has no complete twenty-year result because follow-up would extend beyond tau=25.
The fixed tail, final-slope continuation, and decline-to-0.0001 with decay 0.1/year
are explicitly selected alternatives, not evidence-supported central defaults.

Numerical step refinement from 0.125 to 0.0625 years changes risk by less than
0.00000016 at every exported horizon; the predeclared tolerance is 0.0001.
This is numerical error assessment, not an uncertainty interval.

## Competing mortality

For 1,000 synthetic 70-year-olds with tau=30, RR=1 and the fixed remote tail,
twenty-year expected cases decrease from 7.571193 without competing death to
3.780296 with a synthetic mortality hazard of 0.08/year. The latter produces
795.851784 competing deaths. The cause-specific rate remains 38 per 100,000
disease-free person-years in this constant-hazard example; cumulative risk and
case counts fall because mortality removes future person-time.

The stochastic diagnostic uses continuous event-time sampling with reinfection
and explicit synthetic mortality, seed 31829 and 4,000 samples. It checks first
TB, death, disease-free survival and disease-free person-time against the matching
deterministic projection. Tolerances are predeclared as five Monte Carlo SEs plus
0.0001 probability numerical allowance (scaled by horizon for person-time).
All four checks pass. This validates the implementation, not source assumptions.

## Invented Australian subgroup mixture

The retained 3A mixture contains 8,000 Australian-born, 1,000 recent-arrival and
1,000 long-term-migrant residents, with different invented TBI/risk inputs.
This forward demonstration explicitly supplies synthetic conditional age/time
support and synthetic age-specific mortality; it does not infer these from a
display band or a summary mean. Infected remote support for older cells lies
beyond the last evidence knot, so the numeric table uses the explicit fixed-tail
sensitivity.

| Follow-up | Expected cases, 10,000 baseline people | Competing deaths |
| ---: | ---: | ---: |
| 1 year | 4.983872 | 49.858810 |
| 5 years | 11.549682 | 337.190785 |
| 10 years | 14.545588 | 2,033.052109 |
| 20 years | 17.871219 | 4,027.792708 |

Subgroup/sex/RR/state totals reconcile by counts. Under the central no-extrapolation
policy, 527.550385 people (5.2755%) remain unresolved at horizons 1–10, increasing
to 689.695918 (6.8970%) at twenty years. The full result is unavailable. Covered
sets change with horizon, so common-horizon counts are supplied separately.

## Diagnostic Viet Nam comparison

A separate invented population of 1,000, with ages 30 and 60, TBI 20%, recent 5%
and remote-only 15%, uses synthetic compatible heavy-smoking prevalence 1% in
females and 20% in males. It is not a national demographic or survey sample.
Under the diagnostic curve and fixed-tail sensitivity, instantaneous baseline
incidence is 242.411388 per 100,000 disease-free person-years. The frozen WHO
2025 snapshot's 2024 estimate is 182; difference 60.411388 and ratio 1.331931.
The original WHO uncertainty bounds and row provenance remain in the comparison.
No fitting occurs and numerical proximity carries no validation of country risk.

One-year cumulative cases are 1.451467; twenty-year cumulative cases are 5.046120.
The no-mortality label, first-episode limitation, synthetic population mismatch,
calendar-time mismatch, WHO recurrence and other exclusions are explicit.

## Children and historical conditioning

The child diagnostic has 250 people aged 3, 250 aged 10 and 500 adults. Both child
groups remain unresolved, so 500/1,000 (50%) block complete incidence at every
horizon. Resolved outputs for synthetic child curves exist only in tests.

The historical sensitivity starts with 9,500 recent RR-1 and 500 remote RR-74
adults. Hypothetically removing everyone with pre-baseline TB according to the
same surrogate hazards removes 680.264223 people and changes one-year expected
cases from 234.920504 to 212.505883. This is an illustrative depletion calculation,
with a different eligible denominator; it is not a country correction, survival
model, treatment model or recurrent-TB model. The source baseline remains intact.
