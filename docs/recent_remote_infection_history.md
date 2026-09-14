# Recent-versus-remote LTBI infection-history scenarios

This note documents the explicit recent/remote LTBI scientific scenario used by
the APY Python model. It is separate from the SA Health MATLAB-v9-compatible
reference anchor.

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

Recent infection is initially defined as infection acquired within the previous
two years. Remote infection is infection acquired earlier than that window.

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

The explicit infection-history pathway is a scientific scenario or sensitivity
analysis. It is intended to show how results change when recent versus remote
LTBI is derived from a transparent historical infection-pressure assumption.

## Limitations

- Baseline or near-baseline active TB is not automatically evidence of future
  progression from remote LTBI.
- The trajectory slopes are scenario assumptions pending better APY-specific
  epidemiological history.
- Contact and risk-factor information does not by itself identify infection
  timing.
- Wider sensitivity analysis should vary the trend magnitude and the
  recent-infection window.
