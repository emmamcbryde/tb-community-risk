# Incidence-trend estimation: methods and interpretation

Status:
* log-linear and penalised smooth trends: **implemented but provisional**, pending
  scientific review;
* state-space trend: **designed but not implemented**.

Code: `engine/who_incidence/trend.py`.

**Scope.** The trend summarises *estimated TB disease incidence*. It is not a
transmission or infection-pressure trend. It is descriptive and is not used by
any model in this version.

## Log-linear recent trend (primary)

`log(I_t) = a + b (t - t_mean) + e_t`, fitted by ordinary least squares to the WHO
point estimates in the fitting period.

* **Annual percentage change (APC):** `100 (exp(b) - 1)`.
* **Fitting period:** by default the 10 years ending at the latest WHO year.
  Alternatives are 5 years, 15 years, the full series, or custom first and last
  years. Missing years inside the period are reported and never imputed.
* **Minimum data:** at least 5 years are recommended. With 3-4 years the result is
  flagged as a short series; with fewer than 3 no trend is estimated.
* **Regression interval:** the t-based 95% interval for `b`, transformed to the APC
  scale.
* **Diagnostics:** R-squared and residual SD on the log scale, Durbin-Watson,
  largest residual, and the p-value of a quadratic term (curvature test).

## Segmented disruption option

`log(I_t) = a + b (t - t_mean) + c D_t + e_t`, where `D_t = 1` in the nominated
disruption years (by default 2020-2022). This is a temporary level shift. The APC
comes from `b`, and `100 (exp(c) - 1)` is reported as a level shift. It is **not**
a causal estimate of the effect of COVID-19.

## Penalised smooth trend (sensitivity)

This is a Whittaker-Henderson smoother, a discrete P-spline, on log incidence over
the annual grid:

`minimise sum w_t (y_t - s_t)^2 + lambda sum (Delta^2 s_t)^2`

* **Weights:** 0 for missing or excluded years, 1 otherwise.
* **Smoothing parameter:** `lambda` is chosen by generalised cross-validation on
  the grid 10^0 to 10^6 (quarter-log steps), never below 10.
* **Reported with every result:** the chosen `lambda` and the effective degrees of
  freedom.
* **Recent slope:** the mean annual change of `s` over the last 3 years of the
  period, transformed to an APC.
* **Boundary caution:** smooth slopes are least stable at the end of a series, so
  every smooth result carries a boundary warning.
* **No extrapolation:** the smooth is not extended beyond the observed years.

## Uncertainty propagation from WHO bounds

1. On the log scale, `sigma_low = (log I - log lo) / 1.96` and
   `sigma_high = (log hi - log I) / 1.96`. The WHO bounds are treated as the 2.5%
   and 97.5% points of a split-normal distribution. The asymmetry is kept.
2. Draw 1,000 series with a fixed recorded seed (default 20250101). Years are
   independent by default, which gives the most conservative slope uncertainty. A
   `common_shift` option treats years as fully correlated.
3. Refit each draw with the same method, settings and `lambda`.
4. The 2.5th-97.5th percentiles of the refitted APC form the **propagated
   uncertainty interval**.

What the interval covers: only the published WHO bound uncertainty. It excludes:
* WHO methodological uncertainty;
* future structural change;
* model structural uncertainty;
* infection-pressure uncertainty;
* all other parameter uncertainty.

It is not a confidence interval. If any year in the period lacks bounds, or has a
lower bound of zero, propagation is not available.

## COVID-era handling

| Option | Effect |
| --- | --- |
| Include all WHO estimates (default) | All years are used. |
| Exclude 2020-2022 | Those years get zero weight or are dropped. |
| Segmented disruption adjustment | Temporary level shift (log-linear only). |
| Exclude selected years | Only the chosen years are dropped. |

When the period includes 2020-2022, the log-linear APC is also computed with those
years excluded. The trend is flagged as disruption-sensitive if the two differ by
more than 1 percentage point per year. Changes during 2020-2022 are not assumed to
be artefacts.

## Diagnostic classification (first matching rule)

| Classification | Criterion |
| --- | --- |
| Unsuitable for automated trend inference | Fewer than 3 years, zero estimates, or inconsistent bounds |
| Unstable because of short series | Fewer than 5 years (log-linear) or fewer than 8 (smooth) |
| Possible series discontinuity | A year-on-year change of more than 50% |
| Disruption-sensitive | As above |
| Non-linear trend - review smooth estimate | Curvature p < 0.05 with at least 6 years, or residual SD > 0.10 on the log scale |
| Bounds incomplete | Any year lacks usable bounds |
| Smooth estimate - boundary slope | Every smooth result |
| Adequate log-linear fit | None of the above |

The estimation-method field is not published in the public CSV, so a method change
cannot be confirmed; only possible discontinuities are flagged.

## Qualitative summary

The summary uses the propagated interval, or the regression interval if bounds are
unusable:

* **Estimated disease incidence increasing:** the whole interval is above 0.
* **Estimated disease incidence decreasing:** the whole interval is below 0.
* **No clear recent change:** the interval lies within ±1% per year. This is not
  evidence that the epidemic is constant.
* **Trend uncertain:** otherwise.

A summary is always shown with the APC, the interval, the period, the method and
the fit classification.

## Examples (WHO 2025 round, 10-year log-linear, COVID years included)

| Country | APC | Propagated interval | Classification | Summary |
| --- | --- | --- | --- | --- |
| South Africa | -9.3% | -13.4% to -4.3% | Adequate log-linear fit | Decreasing |
| Indonesia | +2.3% | +0.7% to +3.5% | Non-linear | Increasing |
| Philippines | +1.7% | -3.8% to +8.0% | Adequate log-linear fit | Uncertain |
| Australia | +0.1% | -2.0% to +2.3% | Adequate log-linear fit | Uncertain |
