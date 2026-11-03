# Scientific definition and evidence status

## Four distinct quantities

1. Infection model: probability of true TB infection conditional on available
   evidence. A positive test is an observation, not certainty of infection.
2. Progression model: probability of first active TB conditional on infection
   and remaining disease-free at the observation date.
3. Absolute risk: infection probability times conditional progression risk.
4. Ranking: ordering that absolute risk over a common horizon, with missing
   information and uncertainty visible.

For time since most recent infection tau and horizon t, define a cumulative
baseline hazard H0. Its future increment is
`DeltaH = H0(tau+t) - H0(tau)`. Without competing death:

```
P(TB by t | infected, evidence) = -expm1(-m * DeltaH)
P(TB by t | evidence) = P(TBI | evidence) * P(TB by t | infected, evidence)
```

Knots contain cumulative **hazards**, not risks. Linear interpolation in cumulative
hazard gives constant non-negative segment hazards. Curves must start at (0,0),
with strictly increasing durations and nondecreasing cumulative hazards.
Integration sums segment rate times duration; unsupported tails return unavailable.
No extrapolation never means zero progression. Age-specific curves are selected
using current exact age; this milestone uses that baseline-age applicability
through the horizon and does not implement age-at-infection effects implicitly.
A supplied pediatric curve must be independent of the adult curve.

`recent_infection_progression` explicitly selects the most-recent-infection clock
and recency window. Recent exposure, conversion and close contact do not each
multiply risk again. The person-supplied duration must concern the most recent
infection, including reinfection where known. Historical infection events are
not reconstructed and future events are not generated in this milestone.

## Observation model

A supplied latent-TBI probability may be used directly, with its source retained.
Alternatively, fixed assay sensitivity Se, specificity Sp and a supplied pretest
probability p give, for a positive assay:

```
posterior = Se*p / (Se*p + (1-Sp)*(1-p))
```

For a negative assay replace Se with 1-Se and 1-Sp with Sp. Log likelihoods avoid
underflow. Impossible evidence under the selected likelihood is unavailable.
Assay accuracy is not fitted from the same record. TST threshold, BCG-specific
accuracy, IGRA product and any use of two conditionally independent assays must
be selected explicitly. Unknown, indeterminate and borderline outcomes are not
negative. Product mismatch and contradictory TST size/interpretation are rejected.
There is no automatic BCG/NTM adjustment or country-incidence-to-infection rule.

## Multimorbidity

The working default is `highest_applicable_rr = max(1, RR_i)`. Tied contributors
are retained. Alternatives are explicitly selected:
`additive_excess_rr = 1 + sum(RR_i - 1)` and
`multiplicative_rr = product(RR_i)`. No multiplier cap or automatic policy change
is implemented. Effects must be at least one for these progression risk-factor
policies. Protective effects need a separate reviewed contract.

Relative risks are not automatically hazard ratios. Applying them as multipliers
of cumulative hazard is an explicit approximation recorded in warnings. Multiple
extreme effects trigger a separate warning; the configurable warning threshold
is an engineering diagnostic, not a clinical boundary or numerical cap. Numeric
multiplier overflow returns unavailable rather than choosing a different policy.
Adult effects are not applied outside their declared ages. Context constraints,
such as untreated HIV, must match or remain unavailable. Known unmodelled factors
require an effect or an explicit source-bearing scope exclusion.

## Optional competing death

No-mortality mode must be explicitly selected and is labelled as an omission.
An externally supplied constant death hazard mu is also implemented. Over a
piecewise-constant TB hazard h and interval d, disease increment is
`S*h/(h+mu) * (1-exp(-(h+mu)*d))`, death increment replaces h by mu, and
`S <- S*exp(-(h+mu)*d)`. Disease and death remove subsequent first-event risk.
This alternative does not equal simply multiplying latent probability by a
mortality-free cumulative risk. No life expectancy, country mortality, migration
or recurrence assumptions are inferred from demographic fields.

## Relationship to historical tools

Menzies et al. (2008), *Int J Tuberc Lung Dis* 12:498-505, PMID 18419884,
separates test interpretation from disease risk in an adult algorithm. Its
published parameter table is historical documentation, not validation of today's
conditional hazards. The reference baseline cohort was healthy TST-positive US
military recruits followed for four years. This library does not adopt that
baseline, the Styblo conversion or lifetime-to-age-80 multiplication.
[Primary paper](https://www.tstin3d.com/thinking3D_2008.pdf).

The current TSTin3D page labels itself version 3.0 and notes changes since 2008.
The About page describes IGRA-only positive PPV 0.98, a 0.1% baseline annual
risk, and distinct recent-contact assumptions. These are documented claims,
not verified software equivalence or portable defaults. The [About page](https://tstin3d.com/en/about.html)
and [calculator](https://tstin3d.com/en/calc.html) were accessed 2026-11-03.
The [disclaimer](https://tstin3d.com/en/disclaimer.html) limits use and redistribution.
No website source code, executable behaviour or copied parameter set is included.
A commercial or clinical derivative would require independent terms/licensing
review. The paper and live site's numerical descriptions must not be assumed
identical.

[Canadian Standards 2007, p65](https://publications.gc.ca/collections/collection_2011/aspc-phac/HP40-18-2007-eng.pdf)
reports ranges for several effects, including transplantation 20-74 and renal
failure requiring dialysis 10-25. A selected historic upper value is not a
universal conditional progression hazard ratio. [Current Canadian diagnostic guidance](https://www.canada.ca/en/public-health/services/diseases/tuberculosis/health-professionals/canadian-tuberculosis-standards/diagnosis-tuberculosis-infection.html)
also presents absolute disease rates in specific test-positive groups. Neither
source is silently converted into model inputs.

The [evidence table](evidence.csv) distinguishes a value verified **in a source**
from an effect reviewed for this model. Every clinical component remains
`external_parameter_required`. All example values are invented independently;
none should be used for clinical decisions. No external validation or exact
TSTin3D reproduction is claimed.
