# Milestone 2D: baseline-TBI source and population subgroups

This is an isolated pure Python model/data-contract layer for planning and
sequencing, not denying care. It does not expose production controls or run a
projection, intervention, progression, transmission or economic model.
See [contracts](milestone2d_contracts.md), [synthetic results and validation](milestone2d_validation.md)
and [decisions](decision_register_2d.md).

## Two explicit sources; no fallback

The top-level question is **How should baseline TBI prevalence be determined?**
Exactly one scalar mode is required. Missing/unknown/multiple modes, mixed
configuration payloads and attempts to supply an executable forward operator
fail. A blocked or incompatible source is never replaced by another source.

| Mode | User-facing option | Behavior in 2D |
| --- | --- | --- |
| `active_tb_incidence_inference_v1` | **Estimate from active TB incidence** | Validated inputs; structured unavailable result; no TBI estimate |
| `observed_or_user_tbi_prevalence_v1` | **Use observed or user-specified TBI prevalence** | Reviewed overall/age anchors or direct profile inputs |

The first description is: ?Infer the baseline TBI profile using active TB
incidence, demographic structure, infection-pressure history, progression
natural history and population risk factors.? This describes the future
method; its operational flag is false.

The second description is: ?Use a reviewed TBI prevalence estimate or profile
for this population or subgroup. Where appropriate, the model can use the
historical infection-pressure shape to distribute an overall prevalence across
age and recent/remote infection states.?

Each selection records mode, residence geography/level, subgroup, 2024 baseline,
source, review, provenance/hash, profile classification, recency-composition
classification, configuration and scientific hash. A selection is pending
resolution and not simulation-ready. The separate resolved artifact records
actual component coverage and usability. `profile_kind` distinguishes
empirical, user-defined, model-inferred and diagnostic profiles. Derived
anchor profiles are model-inferred (or diagnostic); their input evidence kind
is retained separately. Source review labels are supplied metadata, not a
software certification of scientific validity.

## Why active-TB incidence inference is blocked

The deliberate status is:

`unavailable_pending_progression_and_risk_forward_operator`

No approved disease-forward operator exists to connect exposure/infection
histories to observed active disease. Required assumptions remain unresolved:
recent and remote progression, age-specific progression, risk-factor
prevalence and effects, competing mortality, migration/imported disease,
incidence ascertainment, and dynamic transmission where applicable.

The active-TB target contract retains target type; a rate or count; units;
positive denominator and denominator unit; observation period; calendar year;
age range; sex; geography/level; uncertainty; source and provenance; WHO-estimate
versus notification/source class; ascertainment status; review and any explicit
observation-model reference. Counts without denominator or observation period
fail. Rate and count inputs are mutually exclusive. A prevalent-disease rate
has person units, not person-year incidence units. Counts may be estimated
case counts; no binomial interpretation is imposed.

The following observation meanings remain distinct:

- incident disease;
- notifications;
- baseline prevalent disease;
- screen-detected prevalent disease;
- imported active disease;
- mixed or unclear observations.

Only a reviewed incident-disease target with reviewed ascertainment can be
marked eligible for a future calibration. Notifications remain ineligible as
raw notifications even if an observation-model reference is supplied; they
are not silently converted to incidence. Prevalent, screen-detected and
imported disease are not locally generated incident progression. Eligibility
never enables computation in 2D: the result contains null TBI probabilities,
null calibration/model profile, the full unresolved denominator and all eight
missing operator components.

No SA Health compatibility pathway, APY default, provisional multiplier,
inherited OR multiplication, Styblo conversion, fixed annual progression
rate, explicit-TBI experiment or Starsim demonstration is used. Tests replace
catalytic computation with a function that raises and verify that selecting
incidence mode still returns its blocked result.

## Functional observed/user-input pathway

All prevalence inputs use proportions in [0,1], not percentages or percent
strings. Invalid and infeasible values fail without clipping. Missing
numerators/denominators remain jointly null; zero observed counts are distinct
from missing counts. Original records and metadata have audit checksums.

### Overall prevalence anchor

Supply population/subgroup, population denominator, prevalence, uncertainty
(or explicit unspecified uncertainty), actual sample numerator/denominator if
known, 2024 year, source, review, provenance and a reviewed observation operator.
The denominator of the population and the denominator of a survey sample are
separate fields. The helper uses 2024 population age/sex cells as weights and
requires `apply_anchor_to_model_covered_population_v1`. This explicitly applies
the overall anchor to the age-covered population when an open tail cannot be
estimated. It does not claim to have fitted the complete population.

The retained 2C catalytic model fits one non-negative pressure scale with
bounded root finding. Age-specific TBI, recent, remote-only,
prior-remote-plus-recent and most-recent-event summaries are derived
conditionally on that overall value, pressure proxy, historical backcast,
reinfection assumptions, age structure and any explicit exposure-proxy
application. Actual counts are preserved, never manufactured from prevalence.
This convenience helper does not copy 2024 demographic weights backwards to
an earlier survey year; an earlier observation needs the age-anchor pathway
with explicit reviewed survey weights and actual evaluation time.

A direct latent-TBI operator makes the covered latent prevalence match the
supplied value. If a binary assay operator is explicitly supplied, its target
is observed positivity and the achieved positivity, rather than latent TBI,
matches that value. Se/Sp, assay/product and threshold remain fixed and
provenance-bearing. TST and IGRA are not interchangeable.

### Age-specific anchors

Reviewed single-year or age-band observations can be both-sex, female, male,
or separate female/male records. They use the existing 2C age/sex weight and
assay contracts, with an explicit subgroup ID and source population match.
Actual-count data support binomial likelihood; prevalence/uncertainty data
support explicitly selected WLS. A single suitable observation can use the
bounded-root method. There is one common scale across ages and sexes within
each subgroup, not separate scales that force each observation to match.
Different subgroups may each have their own reviewed source and scale.

Observed/achieved values, residuals, per-observation objective contributions,
feasible ranges, boundary/convergence and incompatibility diagnostics are
retained. An incompatible fit returns its calibration diagnostics, null
resolved TBI values and an unresolved population, with no fallback.

### Direct age/sex profile

A direct profile supplies latent TBI rather than untransformed test positivity.
It is separate from any model-derived comparison profile, which is embedded
unchanged when available. If no comparison exists, its absence is explicit;
the model does not invent one through blocked incidence inference.

A row declares one of three age-allocation policies:

- `exact_single_year_values_v1`: a value for the finite demographic cell
  [age, age+1), not an invented exact birthday;
- `piecewise_constant_within_band_v1`: explicitly use the supplied prevalence
  throughout the named band, including an open band if explicitly supplied;
- `supplied_within_band_distribution_v1`: explicit within-band age/sex values,
  without interpolation. If a parent-band total is also supplied, its
  population-weighted child values must reproduce that total.

A both-sex probability explicitly applies that same probability to the
population's separate sex cells. It does not manufacture sex counts.
Overlapping both/sex-specific rows fail. Separate female/male values are
supported. Every positive population cell must be covered for a complete
profile. Gaps remain visible with unresolved counts; records outside the
population scope and profile boundaries cutting demographic cells fail.
Demographic counts are never split using prevalence values as age weights.
Original numeric notation and all supplied source values remain unchanged in
the source receipt; calculation rows use canonical decimal values.

### Direct recent/remote states

The user can supply all three mutually exclusive probabilities:

\[
p_{recent}+p_{remote\_only}+p_{never\_infected}=1.
\]

Any separately supplied total must equal recent plus remote-only. An optional
prior-remote-plus-recent probability must be no larger than recent. It is a
subset and never added again to total. These direct states are retained even
if a model comparison or conditional-allocation policy is present. Their
composition is labelled observed for empirical inputs and user-defined for
user-defined inputs. Infection-time detail is not inferred from three state
probabilities alone.

### Total only: recency unavailable unless explicitly selected

The default is `recency_unavailable_v1`. A supplied total determines never
infected as its complement; recent, remote-only and the prior-remote subset
remain null. Having a model comparison does not activate allocation.

Only `model_conditional_recency_allocation_v1` opts into the hybrid:

\[
p_{user,recent}=p_{user,TBI}\frac{p_{model,recent}}{p_{model,TBI}},\qquad
p_{user,remote\_only}=p_{user,TBI}-p_{user,recent}.
\]

The auditable prior-remote subset uses the same explicit conditional scaling.
It remains included within recent. This is labelled **hybrid model-derived
allocation**, not observed recency. The comparison must match population,
age/sex cells, residence, baseline time and window. Positive user TBI cannot
use a zero-probability or unresolved model conditional distribution. Missing
model estimates in empty demographic cells remain unavailable but do not
create missing people. Directly supplied recency is never overwritten.
No infection-time summaries are invented for the override, including the
hybrid. The unchanged comparison retains its own conditional time summaries.

## Subgroups, residence and membership

Subgroups have stable IDs, labels/descriptions, current residence geography
and level, 2024 denominator, age/sex population cells and margins, selected TBI
source and resolved profile, exposure-history source when supplied, review,
provenance/hash and a null risk-factor placeholder. Positive subgroup
population must equal its cells, including the complete open-ended tail.
Both-sex population rows cannot be mixed with female/male counts; sex-specific
age grids must align before age margins can be derived without rebinning.
The declared supplied-cell partition is the population scope, including
restricted occupational or service populations, not an inferred complete
national population.

Birth geography, previous residence, exposure-history geography, arrival
date/year and duration in the current setting are distinct nullable fields.
Missing values stay explicitly null. Date/year contradictions and invalid
durations fail. Current residence is never automatically substituted for
historical exposure. No migration or infection-selective mortality model is
created. The 2C model still weights people alive in 2024 and assumes
infection-independent survival.

If a selected pressure history has a different geography from current
residence, use requires a reviewed `reviewed_subgroup_exposure_proxy_v1`
application with source, provenance and rationale. The original history,
source geography, curves and backcast remain embedded unchanged in the
application provenance. A new supplied scenario history is created for the
receiving subgroup, preserving all lifetime intervals and re-normalizing at
2024. This is an explicit proxy transport assumption, not a reconstruction
of an individual's route, arrival-time exposure or migration process.

Aggregate labels such as migrant and occupational group can overlap. A
mixture therefore requires a reviewed disjoint partition: common partition
ID, unique category IDs, definitions and source evidence, plus an exact list
of selected subgroup IDs. Where complete enumerated person IDs are supplied,
their counts must match the subgroup denominator and duplicated membership
within/across groups fails. Without IDs, disjointness is a reviewed source
declaration, not something software can prove from aggregate labels.
Duplicating a subgroup for editing requires a new stable ID and a reviewed
membership decision before inclusion in a mixture.

## Mixture arithmetic and partial coverage

For compatible non-overlapping groups:

\[
N=\sum_gN_g,\qquad
C_x=\sum_g\sum_{a,s}N_{g,a,s}p_{g,a,s,x},\qquad p_x=C_x/N.
\]

Counts determine weights; simple unweighted subgroup means are never used.
Supplied weights must agree with counts (declared 1e-12 rounding tolerance),
and are retained separately from count-derived weights. Counts are not
adjusted to make rounded weights exact. Unique IDs, membership, positive
populations, year, geography, source/provenance and all age/sex sums are
validated. Subnational containment needs explicit reviewed geography evidence.
Age/sex grids must align; no implicit age rebinning or sex split occurs. An
optional reference age/sex population distribution is reconciled exactly.

Outputs retain each subgroup, combined age/sex cells, total population,
source composition and per-component estimated/unresolved populations,
estimated counts, complete counts and covered/full prevalences. Unavailable
complete values are null. `estimated_count=0` with zero covered population
means no known contribution, not observed absence of infection.

If a subgroup selects blocked incidence inference or otherwise lacks a
usable profile, all its people stay in the denominator. The complete mixture
result is unavailable. Covered contributions are labelled partial; the
unresolved subgroup and population remain explicit. A total-only profile can
make total-TBI coverage complete while recent/remote coverage remains
incomplete. A prior-remote subset can remain unavailable even with complete
exclusive-state coverage. Each component has its own coverage fields.

`complete_mixture_result` refers to complete exclusive baseline-state
coverage. It is never an operational or empirical claim for diagnostic
inputs. `usable_for_downstream_simulation` requires complete total/exclusive
states, reviewed non-diagnostic inputs and observed/user mode. Its scope is
only future baseline categorical-state initialization: duration-dependent
progression, interventions and clinical use remain unapproved. Total-only
profiles are still usable for complete total-TBI summaries when covered,
without pretending that recency is known. All demonstration profiles remain
diagnostic and have simulation usability false.

## Open-ended ages and uncertainty

No exact age or terminal age is assigned to 100+. The tail always contributes
to subgroup/mixture denominators. Catalytic estimates require reviewed
within-cell ages and weights; otherwise that population is unestimated.
Direct user-specified 100+ probabilities can resolve the cell without claiming
an exact attained age. No full-population total-TBI prevalence is reported
while any positive tail/subgroup population remains unresolved.

The layer is deterministic. Original observation uncertainty and existing 2C
parameter-scenario information are preserved. It does not turn subgroup
bounds into a mixture confidence interval or invent covariance, sampling
counts, stochastic histories or empirical anchors. Review of those
assumptions remains separate from population arithmetic.

## Intended interface; specification only

For observed/user inputs, expose these future subchoices:

- Overall TBI prevalence;
- Age-specific TBI prevalence;
- Upload age/sex-specific profile;
- Enter recent and remote infection directly.

For multiple populations, specify Add subgroup, Duplicate subgroup, Remove
subgroup, Edit subgroup, Show subgroup results separately, and Show
population-weighted combined result. The incidence choice must show its
blocked status and dependencies. Total-only uploads must default to recency
unavailable; optional conditional allocation requires explicit selection and
its hybrid label. Exposure geography and current residence need separate
fields. The pure `interface_specification()` artifact describes these controls
and marks incidence inference non-operational. No production navigation or
application imports it.
