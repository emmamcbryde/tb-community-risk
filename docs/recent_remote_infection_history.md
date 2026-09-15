# Recent-versus-remote LTBI infection-history scenarios

This note documents the explicit infection-history diagnostic pathway used by
the APY Python model. It is separate from the SA Health MATLAB-v9-compatible
reference anchor and is not validated for decision-making or health-economic
reporting.

## What the data identify

- IGRA and TST identify immunological evidence of MTB infection; they do not
  identify when infection was acquired.
- The current working LTBI prevalence target is `47 / 624`, or approximately
  7.53% of the baseline population.
- The LTBI age odds ratio identifies an age relationship in prevalent LTBI.
  It is not evidence that infection is recent.
- Prevalence plus an age odds ratio cannot identify a calendar-time infection
  trend. A rising, steady or falling transmission-history trajectory is a
  scenario assumption.

## Mathematical model

For a person of current age `a`, infection pressure is represented over
lookback time `tau`, where `tau = 0` is the present and larger values are
further in the past:

```text
lambda(a, tau) = s * exp(-g * tau) * (a - tau + 0.5)^gamma
```

`s` is the fitted hazard scale, `gamma` is the fitted age-shape parameter and
`g` is the fixed scenario trend:

- Rising: `g = +0.01 per year`
- Steady: `g = 0 per year`
- Falling: `g = -0.01 per year`

Positive `g` means infection pressure is lower in the past and higher toward
the present. Negative `g` means infection pressure was higher historically and
has declined toward the present.

The probability of prevalent LTBI by age is:

```text
P(infected by age a) = 1 - exp[- integral(lambda(a, tau), tau = 0..a)]
```

Infection acquired within the preceding two years is distinguished from
infection acquired earlier than that window. This two-year timing definition is
not equivalent to the early higher-progression-risk state used by the current
natural-history equations, which has a five-year mean residence time before
transition to the later lower-risk progression state.

## Calibration

For each fixed trajectory, the model calibrates:

- hazard scale `s`;
- age-shape parameter `gamma`.

The targets are:

- overall baseline LTBI prevalence, `47 / 624`;
- the configured LTBI age odds ratio.

The trajectory slope itself is not fitted because the available data do not
identify it. The preset slopes are conservative scenario assumptions, not
observed APY transmission trends.

## Relationship to the SA Health compatibility reference

The SA Health working-reference economic package continues to use the
MATLAB-v9-compatible stochastic epidemiological anchor. That anchor preserves
the earlier software semantics and must not be interpreted as a measured
recent-LTBI fraction.

The explicit infection-history pathway is currently an experimental diagnostic.
It is intended to show how recent-versus-remote LTBI composition can be derived
from a transparent historical infection-pressure assumption. It must not be
used as the SA Health report reference or as report-facing health-economic
evidence until the progression calibration has been scientifically corrected.

## Deterministic and stochastic reconciliation

The frozen SA Health report reference is the stochastic MATLAB-v9-compatible
Python anchor with 2,000 repetitions and seed 1. Its comparator active-TB mean
is approximately 37.305 cases over 20 years.

A read-only diagnostic at commit `6171d13` found that the deterministic
compatibility expected-value path gives approximately 30.634 comparator cases
for the same headline scenario. This deterministic value is therefore not a
drop-in replacement for the frozen stochastic reference. The current difference
is retained as a documented approximation/reconciliation gap pending later
scientific and numerical review; no report-facing result should substitute the
deterministic compatibility value for the stochastic reference mean.

## Limitations

- Baseline or near-baseline active TB is not automatically evidence of future
  progression from remote LTBI.
- The inherited `10/770` active-TB calibration target is unresolved and is not
  validated as future two-year progression from baseline LTBI.
- Prevalent or near-baseline active TB is not yet separated from future
  incident, preventable TB in the experimental pathway.
- Current recent and remote progression hazards are software-calibrated
  diagnostic quantities, not reviewed natural-history estimates.
- Disease-risk odds ratios are still multiplied as hazard multipliers; extreme
  joint risk-factor combinations can imply near-certain progression.
- The trajectory slopes are scenario assumptions pending better APY-specific
  epidemiological history.
- Contact and risk-factor information does not by itself identify infection
  timing.
- Wider sensitivity analysis should vary the trend magnitude and the
  recent-infection window.
