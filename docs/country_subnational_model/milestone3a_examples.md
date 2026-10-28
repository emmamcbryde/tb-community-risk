# Milestone 3A invented diagnostic examples

**All values are invented. These are not Australian or Viet Nam estimates.**
Legacy adult RRs are enabled explicitly for diagnostics. Demography uses exact
synthetic ages 20 and 60, equally divided by sex and age, with explicit zero-count
100+ cells; these ages are not band midpoints. Positive-tail handling has separate
tests. Baseline anchors use the retained one-scale catalytic model. Small solver
residuals are retained, not altered by risk allocation; exports retain full values.

## Invented Australian subgroup mixture

| Subgroup | Population | TBI anchor | Diabetes | >1 pack/day smoking | Haemodialysis | c | Any-risk union | Joint TBI/any risk |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Australian-born | 8,000 | 4% | 8% | 10% | 1% | 0 | 18.028% | 0.72112% |
| Recent arrivals | 1,000 | 30% | 4% | 20% | 0.5% | 0.5 | 23.584% | 13.7410817% |
| Long-term migrants | 1,000 | 18% | 15% | 3% | 2% | -0.5 | 19.199% | 1.72791% |

Population-weighted TBI is approximately 8% (800 expected people), any risk
18.7007% (1,870.07 people), joint TBI/any risk 2.12379517% (212.379517 people),
and multimorbidity 0.9416% (94.16 people). Expected factor count is 0.1965 per
person. No RR or association control is averaged before allocation.

## Synthetic Viet Nam heavy smoking by sex

The compatible, invented >1-pack/day prevalence is 1% in females and 20% in males.
With equal sex populations, the union is 10.5%; TBI is approximately 20% and joint
prevalence under within-cell independence is 2.1%. General current smoking is
demonstrated separately as a definition mismatch: all 1,000 people stay unresolved
for risk allocation. No current-smoking survey value is adopted or transported.

## Multimorbidity and RR combinations

For diabetes 20%, heavy smoking 10% and haemodialysis 5%, independent overlap gives
union 31.6%, multimorbidity 3.3%, and expected factor count 0.35. Factor counts are
68.4% with none, 28.3% with one, 3.2% with two and 0.1% with all three. Each
original marginal is unchanged after coupling. The highest-RR distribution is
68.4% at RR 1, 7.6% at RR 3, 19% at RR 3.6, and 5% at RR 25.3. For the all-three
pattern the central RR is 25.3; calculation-only multiplicative and additive-excess
sensitivities are 273.24 and 29.9. Neither can enter the central allocation API.

## Feasible control versus achieved phi

For `p=0.3`, `q=0.2`, joint bounds are 0 and 0.2; independence is 0.06.

| Requested normalized c | Joint prevalence | Achieved phi |
| ---: | ---: | ---: |
| -1 | 0 | -0.32732684 |
| 0 | 0.06 | 0 |
| 0.5 | 0.13 | 0.38188131 |
| 1 | 0.2 | 0.76376262 |

c is not an achieved correlation. Every scenario preserves TBI 30% and risk union
20%. State-specific and factor tests preserve recent, remote-only, prior-remote
subset and factor marginals exactly. Aggregate phi can be nonzero with c=0
everywhere because age/sex/subgroup structure induces association; covariance
components are reported separately.

## Incomplete and pediatric examples

Selecting incidence inference for the 8,000 Australian-born people leaves them
unresolved and the complete 10,000-person result unavailable. Known risk union
remains a separately labelled component; full-population TBI and joint TBI/risk
are unavailable. There is no source-mode fallback.

The children example uses exact ages 3, 10 and 20, population 1,200. Adult inputs
do not cover the 800 children and adult effects are not applied. Those 800 stay
unresolved rather than risk-free. Tests permit an explicitly reviewed synthetic
pediatric effect without altering adult evidence.

Positive-tail tests retain 100 of 1,000 people unresolved when the catalytic
source has no reviewed 100+ distribution. Directly supplied TBI states can resolve
the same open cell without inventing an exact age.
