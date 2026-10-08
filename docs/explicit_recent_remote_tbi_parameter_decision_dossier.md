# Explicit recent/remote TBI parameter decision dossier

Status: decision-support documentation only. This dossier does not approve
production defaults and does not change model code, runner behavior, Streamlit
UI, event ledgers, economics, DALYs, MATLAB, frozen artifacts or the dynamic
model.

Branch context: `feature/explicit-recent-remote-tbi-v2`.

Access date for live source checks: 2026-10-08.

## Executive recommendation

The evidence is strongest for a declining progression risk with time since
infection, not a biologically sharp recent/remote boundary. For initial runner
integration, the most defensible decision set is:

1. Retain the five-year recent-infection classification for the explicit TBI
   prevalence targets, but add two-year and continuous-decline sensitivity
   analyses before production use.
2. Use externally supplied natural-history hazards for the first production
   integration, and treat active-TB observations as validation unless a
   genuinely prospective incident target with fixed ascertainment is approved.
3. Use no risk-factor progression multipliers in the central production policy
   until individual factors have been reviewed as hazard multipliers.
4. Add competing mortality through externally supplied age/sex survival curves
   or death hazards. Use ABS life tables for Australian examples and a generic
   WHO life-table interface for international examples.
5. Do not treat an observed active-TB case count plus population denominator as
   sufficient to identify progression hazards.

## Evidence sources used

| source | model use | access date |
| --- | --- | --- |
| WHO consolidated guidelines on tuberculosis, module 1: prevention, tuberculosis preventive treatment, second edition, 2024, https://www.who.int/publications/i/item/9789240096196 | TB infection burden, lifetime progression framing, high-risk groups | 2026-10-08 |
| Behr, Edelstein and Ramakrishnan 2018, "Revisiting the timetable of tuberculosis", https://pubmed.ncbi.nlm.nih.gov/30139910/ | interpretive review of timing after infection | 2026-10-08 |
| Menzies et al. 2018, progression assumptions in TB transmission models, https://pmc.ncbi.nlm.nih.gov/articles/PMC6070419/ | model assumptions and caution on natural-history calibration | 2026-10-08 |
| Time-since-infection natural-history synthesis for the United States, https://pmc.ncbi.nlm.nih.gov/articles/PMC7707158/ | numerical central progression anchor points by time since infection | 2026-10-08 |
| Trauer et al., "Risk of Active Tuberculosis in the Five Years Following Infection", https://pmc.ncbi.nlm.nih.gov/articles/PMC7190060/ | higher-progression sensitivity anchored to five-year risk after infection | 2026-10-08 |
| Late-reactivation systematic review, https://pmc.ncbi.nlm.nih.gov/articles/PMC6850269/ | remote/late progression upper context | 2026-10-08 |
| Jeon and Murray 2008, diabetes mellitus and active TB systematic review, https://doi.org/10.1371/journal.pmed.0050152 | diabetes risk-factor evidence | 2026-10-08 |
| Bates et al. 2007, tobacco smoke and tuberculosis systematic review, https://pubmed.ncbi.nlm.nih.gov/17325294/ | smoking risk-factor evidence | 2026-10-08 |
| Lonnroth et al. 2008, alcohol use and tuberculosis systematic review, https://pmc.ncbi.nlm.nih.gov/articles/PMC2533327/ | harmful alcohol risk-factor evidence | 2026-10-08 |
| Martinez et al. 2020, child close-exposure individual participant meta-analysis, https://pmc.ncbi.nlm.nih.gov/articles/PMC7289654/ | contact and child-age risk evidence | 2026-10-08 |
| CKD without kidney failure TB-risk systematic review/meta-analysis, https://pmc.ncbi.nlm.nih.gov/articles/PMC10573716/ | CKD risk-factor evidence | 2026-10-08 |
| CKD/dialysis TB-incidence systematic review/meta-analysis, https://pubmed.ncbi.nlm.nih.gov/35609860/ | dialysis and advanced kidney disease evidence | 2026-10-08 |
| Chronic airway disease and active TB systematic review, https://pmc.ncbi.nlm.nih.gov/articles/PMC9070518/ | COPD/chronic airway disease evidence | 2026-10-08 |
| ABS life expectancy/life tables, latest available series checked as 2022-2024, https://www.abs.gov.au/statistics/people/population/life-expectancy/2022-2024 | Australian all-cause mortality source candidate | 2026-10-08 |
| WHO Global Health Observatory life tables, https://www.who.int/data/gho/data/themes/mortality-and-global-health-estimates/ghe-life-expectancy-and-healthy-life-expectancy | international mortality interface candidate | 2026-10-08 |

## Progression after infection

Values below are reported only where the cited source supports them. Derived
segment hazards are shown later in the candidate parameterisations.

| study/reference | setting and population | infection definition and follow-up | infected/events | risk by time since infection | measure reported | competing mortality and reinfection | limitations/applicability |
| --- | --- | --- | --- | --- | --- | --- | --- |
| WHO TPT guidance, 2024 | global programmatic guidance | TB infection; lifetime framing | not a cohort estimate | lifetime progression commonly framed as about 5-10% among infected people | guidance summary, not a hazard | not a reusable competing-risk model; reinfection not parameterised for this model | useful for broad plausibility, not for fitting early/remote hazards |
| Time-since-infection natural-history synthesis | United States model synthesis, adult newly infected without other progression risk factors | recent infection with modelled risk over 25 years | synthesis/model-derived, not a single observed cohort numerator | cumulative risk anchor points: 3.8% at 1 year, 5.0% at 2 years, 6.6% at 5 years, 7.2% at 10 years, 7.9% at 25 years; annual rates decline from about 38 to 0.38 per 1,000 person-years between the first and 25th year | model-derived cumulative incidence and rates | competing mortality handling is source-model specific; not adopted as a mortality input here; reinfection not resolved for this branch | best numerical anchor for a central time-since-infection shape, but not directly observed hazards |
| Trauer et al. five-year risk after infection | infected contacts, with higher child risks highlighted | infection/conversion after contact; about 1,650 days | reported cohort included 613 infected contacts and 67 TB events | five-year cumulative hazard reported as 11.5%, with imputation-adjusted estimate 14.5%; risk concentrated early, especially in young children | cumulative hazard/risk after infection | does not provide a production reinfection-reset rule for this model | useful for higher-progression sensitivity and child/contact sensitivity, not a universal adult default |
| Late-reactivation systematic review | untreated cohorts eligible for late/reactivation assessment | later TB after initial high-risk period | varies by included cohorts | rates decline over time and are approximately 200 per 100,000 person-years or lower by the fifth year in eligible untreated cohorts; evidence beyond 10 years is limited | incidence rates from review synthesis | no reusable branch-specific competing mortality curve; reinfection/exogenous infection difficult to separate | useful remote-rate context; not enough to define all late hazards alone |
| Behr et al. timetable review | interpretive review across historical and epidemiologic evidence | infection-to-disease timetable | not a single cohort denominator | emphasizes that much disease occurs soon after infection and that "remote reactivation" may be over-attributed | interpretive synthesis | reinfection and misclassified timing are central concerns, not resolved parameters | supports caution against a flat five-year biological risk window |

### Risk window interpretation

The model must keep three concepts separate:

- Recent infection classification: a person had at least one infection in the
  recent window used for TBI prevalence targets.
- Duration of elevated progression risk: the biological period during which
  risk remains materially higher after infection.
- Mathematical approximation: the constant, piecewise-constant or continuous
  hazard shape used for computation.

The evidence supports a continuous decline more strongly than a hard boundary.
A two-year window captures the highest-risk core and aligns with common
contact-investigation language. A five-year window captures most of the
near-term cumulative excess risk and matches the current explicit recent/remote
prevalence target definition, but it is less biologically precise if treated as
a flat hazard. The recommended structural choice is to retain five years for
central recent-infection classification and require two-year and
continuous-decline sensitivity analyses before production reliance.

## Candidate progression parameterisations

These are decision candidates, not implemented defaults. Percentages are
cumulative progression probabilities without competing mortality.

### Candidate definitions

| candidate | recent-window duration | early formulation | remote hazard | early:remote relationship | source and derivation | evidence quality |
| --- | --- | --- | --- | --- | --- | --- |
| Conservative | five-year classification | constant early hazard `0.0034/year` through five years | `0.00038/year` | about 8.9 if collapsed to two hazards | remote value uses the low late point from the time-since synthesis; early value is a deliberately lower sensitivity that implies about 1.7% five-year risk from infection at baseline | sensitivity only; not directly reported as a natural-history hazard |
| Central | five-year classification, but biological risk represented as piecewise decline | `0-1y: 0.03874083/year`; `1-2y: 0.01255247/year`; `2-5y: 0.00566185/year`; `5-10y: 0.00128894/year`; `10-25y: 0.00050478/year`; `25+y: 0.00038/year` | if forced into the existing two-hazard approximation, use a review-only representative `0.00076/year`; otherwise retain the piecewise shape | collapsed 0-5 average hazard is `0.01365577/year`; ratio to `0.00076/year` is about 18.0 | segment hazards derived from the reported cumulative risks of 3.8%, 5.0%, 6.6%, 7.2% and 7.9% at 1, 2, 5, 10 and 25 years using `lambda = -log(S_b/S_a)/(b-a)` | best-supported working choice, but still model-derived |
| Higher-progression sensitivity | five-year classification, with a front-loaded high-risk period | `0-1y: 0.060/year`; `1-5y: 0.02416345/year`, chosen to match 14.5% by five years | `0.002/year` | front-loaded, not a single stable ratio | five-year cumulative target anchored to the imputation-adjusted Trauer estimate; remote value uses the late-reactivation review's approximate fifth-year upper context | sensitivity only; not a reviewed default |

Transformation formulas:

```text
cumulative probability P(T <= t) = 1 - exp[-H(t)]
piecewise segment hazard lambda_ab = -log(S_b / S_a) / (b - a)
collapsed five-year hazard lambda_0_5 = -log(1 - P_5) / 5
```

### Implied cumulative progression risks

| candidate/state | 1y | 2y | 5y | 10y | 20y |
| --- | ---: | ---: | ---: | ---: | ---: |
| Conservative, infection at baseline | 0.339% | 0.678% | 1.686% | 1.872% | 2.244% |
| Conservative, infection 2.5y before baseline | 0.339% | 0.678% | 0.941% | 1.129% | 1.504% |
| Conservative, infection 4.9y before baseline | 0.068% | 0.106% | 0.220% | 0.409% | 0.787% |
| Conservative, remote-only infection | 0.038% | 0.076% | 0.190% | 0.379% | 0.757% |
| Central, infection at baseline | 3.800% | 5.000% | 6.600% | 7.200% | 7.667% |
| Central, infection 2.5y before baseline | 0.565% | 1.126% | 1.723% | 2.162% | 2.655% |
| Central, infection 4.9y before baseline | 0.172% | 0.301% | 0.686% | 0.944% | 1.443% |
| Central, remote-only infection using representative remote hazard | 0.076% | 0.152% | 0.379% | 0.757% | 1.509% |
| Higher sensitivity, infection at baseline | 5.824% | 8.072% | 14.500% | 15.351% | 17.027% |
| Higher sensitivity, infection 2.5y before baseline | 2.387% | 4.718% | 6.332% | 7.264% | 9.100% |
| Higher sensitivity, infection 4.9y before baseline | 0.421% | 0.620% | 1.214% | 2.197% | 4.134% |
| Higher sensitivity, remote-only infection | 0.200% | 0.399% | 0.995% | 1.980% | 3.921% |

Interpretation: the central evidence-supported shape is not a flat early
hazard. If the first runner integration only supports a two-phase hazard,
that should be labelled a lossy approximation and compared with the piecewise
shape above.

## Calibration policy recommendation

| policy | production recommendation | why |
| --- | --- | --- |
| A: externally supplied hazards | recommended initial production mechanism | avoids fitting unidentifiable progression parameters from weak observation targets; supports transparent sensitivity choices |
| B: fixed early-to-remote ratio with fitted scale | allowed only for a reviewed eligible prospective incident target with fixed ascertainment `q` | one parameter can be fit, but `R`, `q`, target eligibility and denominators must be externally reviewed |
| C: validation only | recommended for retrospective notifications, prevalence rows, screening prevalence and mixed observations | prevents inappropriate conversion of incompatible data into progression hazards |
| D: jointly estimated early and remote hazards | not recommended for initial production | one aggregate target cannot identify two hazards; multiple targets need structural and practical identifiability review |

An observed case count and population denominator are not sufficient. A fitting
target also needs prospective timing from baseline, numerator composition,
source denominator, baseline active-TB exclusions, at-risk population,
ascertainment, observation duration, population composition and an explicit
observation model.

Changing denominators must trigger a recalculation or recalibration of derived
quantities. It must not overwrite the original observed numerator,
denominator or observation period.

## Risk-factor evidence table

| factor | evidence and reported measure | pathway represented | defensible as progression hazard multiplier? | overlap and status |
| --- | --- | --- | --- | --- |
| Age | WHO TPT guidance and child-contact IPD evidence support higher risk in young child contacts; adult biological progression effects are less cleanly separable | age affects accumulated infection probability, progression biology, mortality and comorbidity prevalence | not as one generic multiplier without age-specific model review | reviewed candidate for age sensitivity; central acquisition uses age/time alive and mortality uses life tables |
| Diabetes | Jeon and Murray reported cohort-study RR 3.11 (95% CI 2.27-4.26) for active TB among people with diabetes | active-TB disease incidence, usually without known infection timing | possible reviewed candidate, but RR is not automatically a hazard ratio | overlaps with age, kidney disease and other comorbidities; sensitivity only until approved |
| Chronic kidney disease | CKD reviews show increased TB risk, with estimates dependent on CKD stage and setting | active-TB incidence among CKD populations | possible candidate by CKD stage only after source-specific review | overlaps diabetes, dialysis, age and immunosuppression; unresolved/sensitivity |
| Dialysis | CKD/dialysis review evidence reports high TB incidence, including pooled incidence around 3,718 per 100,000 in CKD populations and higher estimates in dialysis groups | active-TB incidence in advanced renal disease/dialysis | possible high-risk sensitivity candidate, not a generic multiplier | overlaps CKD and diabetes; sensitivity only |
| Smoking | Bates et al. reported associations with infection, pulmonary TB disease and mortality; summarized disease RRs are about 2.3-2.7 in review summaries | mixed acquisition, disease and mortality pathways | not yet defensible as a pure progression multiplier | overlaps alcohol, socioeconomic exposure and lung disease; sensitivity only |
| Harmful alcohol use | Lonnroth et al. reported increased TB risk for heavy alcohol use/alcohol use disorder, with definitions varying across studies | active-TB disease incidence, often confounded | not yet defensible as a pure progression multiplier | overlaps smoking, drug exposure, housing and comorbidity; sensitivity only |
| Drug exposure | no sufficiently specific source was approved for the inherited factor semantics in this branch | unclear; may reflect exposure, social risk, immunologic risk or ascertainment | no | unresolved/unsuitable until source and pathway are identified |
| Close TB contact | WHO and contact studies identify recent exposure/contact as a key acquisition and near-term disease marker | acquisition/recent infection, and possibly shared environment | no as a generic progression hazard multiplier | would double count recent infection if reused carelessly; unsuitable as central progression multiplier |
| Chronic lung or airway disease | COPD/chronic airway disease reviews report heterogeneous active-TB associations, including HRs roughly 1.44-3.14 across some high-income cohorts | active-TB incidence among people with lung disease | possible sensitivity candidate after outcome and confounding review | overlaps smoking and age; sensitivity only |
| Other inherited factors, including `MJ` | inherited code contains factor labels without a reviewed branch-specific source | unclear | no | unsuitable/unresolved |

No OR, RR or incidence-rate ratio should be relabelled as a hazard ratio
without a documented transformation and assumptions. Adjusted estimates must
not be multiplied with other adjusted estimates as if all were independent.

## Joint risk-factor policy

The recommended initial production policy is:

```text
none
```

All progression multipliers should equal one in the central integration.
Reviewed subsets can be added later as explicit sensitivity policies once
their pathway, effect measure, adjustment set and overlap are approved.

Reasons:

- Inherited diagnostic multiplication can reach `2916`.
- A synthetic legacy diagnostic attributed about 83% of expected cases to the
  highest-risk 1% of population weight.
- Calibrating a lower baseline hazard can hide an excessive multiplier
  concentration but does not make the relative-risk structure correct.
- Diabetes, CKD/dialysis, smoking, alcohol, close contact and lung disease can
  overlap biologically and socially.
- Odds ratios are non-collapsible and cannot be treated as independent hazard
  multipliers by default.
- A ceiling on joint multipliers would be a safety/review rule, not biological
  evidence.

If the user wants an early limited sensitivity subset, the most defensible
starting candidates are diabetes, dialysis/advanced CKD and chronic airway
disease, each introduced one at a time or with explicit overlap rules. That is
not recommended as the central policy yet.

## Reinfection

The current isolated assignment design can record `recent` infection with
prior remote exposure and can use the most recent event as the effective
infection clock. The evidence review does not establish that reinfection
universally resets the high-progression-risk clock to the same level as first
infection. Reinfection can contribute to disease in higher-transmission
settings, and prior infection may modify subsequent risk, but the needed
parameter is not resolved for this model.

Recommendation: retain the ability to flag prior-remote-plus-recent infection,
but treat "recent reinfection resets the clock" as a sensitivity decision
requiring user approval before runner integration. Do not present it as an
empirically settled default.

## Age effects

Age should enter the model through separate mechanisms:

- Acquisition calibration: age affects exposure opportunity through time alive
  and the explicit recent/remote windows.
- Progression: do not add a generic age multiplier in the central policy until
  age-specific progression evidence is selected; child-contact and older-age
  sensitivity analyses may be justified later.
- Mortality: age and sex should enter through an external all-cause mortality
  curve or hazard.
- Risk-factor assignment: age-related comorbidity prevalence can affect future
  risk flags, but risk-factor multipliers are not central defaults.

The same observed age association must not be applied simultaneously to
acquisition, progression, mortality and comorbidity without a source-specific
justification.

## Competing mortality

For Australian working examples, use ABS life tables as the preferred
all-cause mortality source. The latest checked ABS life expectancy/life-table
release was 2022-2024. The integration should retain provenance, reference
years, sex resolution, age resolution and the transformation from death
probabilities to annual hazards or survival curves.

For international analyses, use a generic interface that can accept WHO GHO
life-table outputs or a user-supplied survival curve. Do not download or bundle
mortality tables unless redistribution rights and provenance are clear.

Implementation requirements for future integration:

- distinguish period from cohort life tables;
- transform annual death probability `q_x` to hazard as
  `mu_x = -log(1 - q_x)` when a piecewise-constant hazard is needed;
- start survival at one and validate that survival is finite, non-increasing
  and in `[0, 1]`;
- handle the final open age category explicitly, for example by carrying the
  final hazard only when that rule is documented;
- record `competingMortality = not_modelled` when mortality is omitted.

## Active-TB observation decision table

| observation type | fitting permitted? | calibration policy | denominator retained | effective at-risk population | additional information required | validation only? |
| --- | --- | --- | --- | --- | --- | --- |
| Baseline screening prevalence | no | baseline active-TB sequencing or validation | source population denominator | none for prospective TBI progression | prevalence ascertainment, numerator meaning, baseline state sequencing | yes for progression hazards |
| Prospective incident cases after baseline | yes, if all eligibility criteria pass | Policy A validation, or Policy B only with fixed `R` and fixed `q` | source population denominator | non-baseline-active-TB, non-prevalent population with known TBI composition | numerator excludes baseline disease, observation horizon, ascertainment, composition, mortality mode | no, if eligible |
| Retrospective notifications | no | Policy C validation unless retrospective population reconstruction is built | source denominator for the historical period | not the present baseline cohort | migration, turnover, prior treatment, deaths, changing ascertainment, historical infection pressure | yes |
| Mixed prevalent and incident numerator | no | validation or split after source review | source denominator | ambiguous until split | numerator decomposition | yes |
| Changed population denominator | fitting can proceed only after rerun with a new source row and compatible composition | same policy as the observation type after eligibility review | both old and new source denominators retained | recalculated from baseline active-TB exclusions and population composition | whether numerator changed, why denominator changed, source period, ascertainment | depends on eligibility |
| Incomplete ascertainment | yes only if `q` is fixed externally | Policy A or B with explicit `q` | source denominator | at-risk population multiplied by explicit ascertainment in expected cases | `q`, source, review status and uncertainty | no, if otherwise eligible |

Motivating changing-denominator case: if stakeholders change the population
denominator and predicted active-TB detections change, the calibration and
expected-case calculation must rerun using the new population composition and
at-risk count. The original observed numerator, source denominator and
observation period remain visible. If the observed numerator also changes, it
is a new observation row rather than a rewrite of the old row.

## Decisions required before runner integration

1. Natural-history shape.
   Recommended option: central piecewise decline using the reported 1, 2, 5,
   10 and 25 year anchors, with a labelled two-phase approximation only if
   the runner cannot yet support piecewise hazards.
   Alternatives: conservative low progression; higher-progression sensitivity.
   Consequence: determines expected active-TB cases before interventions.
   Sensitivity: yes.

2. Recent-window classification.
   Recommended option: retain five years as the central classification window,
   with two-year and continuous-decline sensitivity analyses.
   Alternatives: switch central classification to two years; abandon a
   discrete classification and use continuous time since infection.
   Consequence: changes state assignment, remaining early-risk duration and
   interpretation of user targets.
   Sensitivity: yes, but it affects calibration targets and cache keys.

3. Progression-calibration policy.
   Recommended option: externally supplied hazards for production plus
   validation-only active-TB comparisons until an eligible prospective target
   is approved.
   Alternatives: Policy B with externally fixed ratio and fixed ascertainment;
   Policy D only after identifiability review.
   Consequence: determines whether active-TB data can move natural-history
   parameters.
   Sensitivity: yes.

4. Ascertainment.
   Recommended option: require explicit `q`, source and review status; no
   silent assumption of 100% ascertainment.
   Alternatives: validation-only when `q` is unknown.
   Consequence: prevents confounding ascertainment with progression scale.
   Sensitivity: yes.

5. Risk-factor progression multipliers.
   Recommended option: `none` for the central production policy.
   Alternatives: reviewed subset sensitivity; legacy OR-as-hazard diagnostic
   only.
   Consequence: controls concentration of expected cases in high-risk strata.
   Sensitivity: yes.

6. Reinfection reset.
   Recommended option: keep prior-remote-plus-recent auditable and run the
   reset assumption only after explicit approval or as sensitivity.
   Alternatives: most-recent-event reset as central; no separate prior-remote
   effect.
   Consequence: affects progression clocks for people with both remote and
   recent infection.
   Sensitivity: yes.

7. Competing mortality source.
   Recommended option: ABS life tables for Australian working examples; WHO
   GHO or user-supplied life tables for international examples.
   Alternatives: no mortality with explicit limitation; bespoke cohort
   survival if supplied.
   Consequence: affects long-horizon expected active-TB cases, especially in
   older groups.
   Sensitivity: yes.

8. Active-TB observation use.
   Recommended option: fit only prospective incident targets with source
   denominator, numerator composition, at-risk population, horizon,
   ascertainment and mortality mode defined.
   Alternatives: validation-only for all active-TB observations.
   Consequence: determines whether local active-TB observations influence
   progression parameters.
   Sensitivity: yes.

Until these decisions are approved, runner integration remains blocked.
