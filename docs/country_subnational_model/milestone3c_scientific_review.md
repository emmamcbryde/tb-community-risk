# Milestone 3C independent scientific review

**Country incidence inversion is not approved.** Numerical inversion is useful
in fully specified fictional closed populations. A root does not approve the
assumptions. Public status remains
`unavailable_pending_progression_and_risk_forward_operator`.

## Independent 3B audit

Starting commit: `213cacda333c5dbc334e6982891aa357bc81ab46`. The audit read the
actual natural-history, policy, baseline, forward and Monte Carlo code. No
retained 2A–3B source/schema/test/preview file is changed.

| Item | Finding |
| --- | --- |
| Hazard conversion | Decimal `-log(1-F)`; cumulative hazard, not probability, interpolated |
| Segment hazards | Differences of cumulative hazard divided by interval length; exact knots recovered |
| Multiplier | `1-exp(-m*dH)`; no multiplication of cumulative probability |
| Reinfection | Repeated events; most recent infection resets central clock; no extra RR or person |
| Recency | Actual last event determines state; prior remote plus recent is a subset of recent |
| Depletion | First TB and competing death absorbing; baseline active disease excluded |
| Future pressure | Explicit exogenous scenario, Strang renewal approximation and refinement; no WHO extrapolation |
| Monte Carlo | Explicit caller RNG; continuous integrated-hazard event-time inversion and first-event accounting |
| Units | Instantaneous rate, prospective cases/person-time and cumulative cases distinguished |
| Isolation | Public source dispatcher never imports inversion; canonical JSON/SHA-256 and protected-path checks |

Fresh reproduction with4,000 MC draws exactly matched export hash
`585539a943d7c1aee57d520aa907734abdce2e0112ac46db3711391aaf8188c4`.
Independent Decimal arithmetic recovered segment hazards (/year):
0.03874082831643056,0.007384503109215967,0.0012889410885284093,
0.0005047797687262486. All requested diagnostics reproduced:

| Diagnostic | Value |
| --- | ---: |
| Australian fictional mixture, one-year cases | 4.9838718713968015 |
| Same, cumulative20-year cases | 17.871218907806064 |
| Viet Nam diagnostic baseline rate/100,000PY | 242.41138771559562 |
| Frozen WHO target | 182 |
| RR74 one-year risk after new infection | 94.31205575309726% |
| Older adult20-year cases/1,000, no mortality | 7.57119302386863 |
| Same, invented mortality | 3.780295971824996 |

No numerical discrepancy was found. `no_extrapolation_v1` raises an
external-required condition, never a zero tail hazard. Numerical agreement
does not establish clinical validity.

## Adult curve: sensitivity only

[Menzies NA et al. *Time Since Infection and Risks of Future Disease for
Individuals with Mycobacterium tuberculosis Infection in the United States*.
Epidemiology2021;32:70–78, doi:10.1097/EDE.0000000000001271](https://pmc.ncbi.nlm.nih.gov/articles/PMC7707158/)
Figure1B gives modelled cumulative risks3.8%,6.6%,7.2%,7.9% at1,5,10,25years.
The two-year5% knot remains unverified and excluded. Reported95% intervals are
3.2–4.5%,5.8–7.7%,6.4–8.3%,7.0–8.9%. The focal person is infected at20 with no
additional progression factors. These are US evidence-synthesis model outputs,
not cumulative risks observed in a representative national infection cohort.

Historical inputs include young UK test converters (Sutherland), household
contacts (Ferebee/Mount), and Netherlands molecularly inferred timing (Borgdorff).
Age distributions, infection ascertainment and settings differ; infection time
is not observed continuously. A common baseline disease exclusion, preventive-
treatment and first-versus-recurrent estimand matching this operator has not
been established here. These are review gaps, not assumed exclusions. The
paper explicitly uses US background mortality for remaining-life outputs and
HIV excess mortality. The four finite-horizon cumulative outputs have not been
independently established as portable cause-specific net risks. Transforming
them to hazard is mathematically reproducible but scientifically conditional.
**Status: sensitivity only; no national central adult curve approved.**

## Remote tail: no central extrapolation approved

[Dale KD, Karmakar M, Snow KJ, Menzies D, Trauer JM, Denholm JT. *Quantifying
the rates of late reactivation tuberculosis: a systematic review*. Lancet
Infect Dis2021;21:e303–e317, doi:10.1016/S1473-3099(20)30728-3](https://doi.org/10.1016/S1473-3099(20)30728-3)
finds sparse untreated known-acquisition follow-up and declining late risk. It
does not identify a universal hazard beyond25years. The Menzies rate near25years
(about0.00038/year) also differs from the average10–25-year segment0.00050478.
Neither establishes a permanent constant tail or zero hazard.

| Sensitivity | Hazard just after25years | Infinite-time net cumulative risk, without death |
| --- | ---: | --- |
| Constant final segment | 0.0005047797687/year | 1 asymptotically |
| Declining to0.0001, decay0.1/year | 0.0005047797687 initially | 1 asymptotically |
| Declining to zero, decay0.1/year | 0.0005047797687 initially | `1-exp[-H25-h_final/0.1]`, about8.36% |
| Supplied fixed remote hazard, illustrative | 0.00038/year | 1 asymptotically |
| Bounded scenarios | Separately declared lower/upper histories | Conditional on each; no clipping/CI claim |
| No estimate | Unavailable | Unavailable |

The net infinite-time limit is **not remaining-lifetime risk**. Finite-life risk
requires age-specific survival and integration of `S_TB*S_death*h_TB`. No
arbitrary terminal age is introduced. Exports report current/50-year hazards,
net limits, older-age incidence and conditional fitted scale/TBI for each tail.

In the reproduced10,000-person Australian toy population, no extrapolation
leaves527.550385 people(5.2755%) unresolved at1/5/10years,689.695918(6.8970%)
at20years. Complete cases are unavailable. Over a finite first-episode horizon,
the unknown contribution is between zero and the unresolved count; zero is a
bound, not an estimate. Instantaneous contribution has no finite upper bound
without a hazard bound. A separate1,000-person,tau10 example gives20-year cases
10.044806 constant-final,9.617831 declining-positive-plateau,9.426981 fixed0.00038,
and unavailable with no extrapolation.

## Pediatric gate

[Martinez L et al. *The risk of tuberculosis in children after close exposure:
a systematic review and individual-participant meta-analysis*. Lancet2020;
395:973–984, doi:10.1016/S0140-6736(20)30166-5](https://pmc.ncbi.nlm.nih.gov/articles/PMC7289654/)
contains137,647 exposed children under19,1,299 prevalent and999 incident cases.
The under-five two-year19% combines prevalent and incident disease;83% of
under-five cases were diagnosed within90days. Contact ascertainment is not
known infection time; baseline TST/IGRA positivity does not date infection.
Contact intensity, BCG, treatment selection, age at infection, ascertainment and
competing death prevent direct adoption as a postinfection progression curve.

[da Costa FBP et al. *Mycobacterium tuberculosis infection and tuberculosis disease
in the first decade of life: a South African birth cohort study*. Lancet Child
Adolesc Health, doi:10.1016/S2352-4642(24)00256-6](https://doi.org/10.1016/S2352-4642(24)00256-6)
adds prospective Drakenstein repeated-TST data. Ongoing high exposure, BCG context,
conversion intervals and selective preventive treatment do not identify complete
untreated age-at-infection/time-since-infection curves, especially for5–14years.
**Both pediatric groups remain external-required.** No adult scaling or invented
clinical defaults. Frozen2024 child/open-tail population shares are quantified
in the diagnostic report; missing child incidence is not recovered by rescaling
an adult denominator. The earlier child toy example leaves500/1,000 unresolved.

## Risk-factor qualification

All11 candidates remain traceable to [Menzies D et al. IJTLD2008;12:498–505,
Table2](https://www.tstin3d.com/thinking3D_2008.pdf). The new evidence contract
retains exact definitions, effects, reference chains, estimand fields, explicit
unknowns, independent comparisons and decisions. None is newly approved as a
conditional-TBI hazard multiplier. Adult-only applicability remains unchanged.

| Candidate | Legacy RR | Independent qualification/primary comparison |
| --- | ---: | --- |
| Granuloma | 2 | Horwitz/Grzybowski chain not fully reappraised; imaging, prior disease/treatment and active-disease exclusion unresolved |
| Smoking>one pack/day | 3 | [Leung2004](https://hub.hku.hk/handle/10722/86866), elderly current-vs-never pulmonary-TB HR2.87; different intensity definition, no conditional-TBI denominator |
| Weight<90%ideal | 3 | Comstock1975 chain incomplete; BMI threshold not equivalent; reverse causation/baseline TB unresolved |
| Diabetes | 3.6 | [Leung2008](https://pubmed.ncbi.nlm.nih.gov/18400769/), older Hong Kong cohort; control-dependent disease association, not universal TBI HR |
| TNF inhibitor | 4 | [Dixon2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2927681/), agent-specific incidence ratios, rheumatoid arthritis selection and screening/TPT era |
| Fibronodular CXR | 6 | Prior disease, lesion classification and preventive-treatment cohorts need extraction; first-episode hazard not established |
| Head/neck carcinoma | 16 | [Simonsen2017](https://pubmed.ncbi.nlm.nih.gov/28366685/), Danish matched cancer cohort,aerodigestive aHR8.12; different site/treatment/reference mixture |
| Haemodialysis | 25.3 | [Dobler2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3246480/),14,506 Australian dialysis patients,mean3.8years,aRR7.8 vs general population; not conditional-TBI HR |
| Silicosis | 30 | [Cowie1994](https://pubmed.ncbi.nlm.nih.gov/7952577/),1,153 older male miners without prior TB,7years,RR2.8 with severity gradient; exposure/LTBI denominator differ |
| HIV | 50 | [Selwyn1989](https://pubmed.ncbi.nlm.nih.gov/2915665/),pre-modern-ART injection-drug-user cohort; CD4,ART,anergy,TPT and reference group matter |
| Transplant/immunosuppression | 74 | [Torre-Cisneros2009 RESITRA](https://pubmed.ncbi.nlm.nih.gov/19445585/),4,388 recipients,reported general-population RR26.6; organ/regimen/prophylaxis differ, not TBI HR74 |

These comparisons are not replacement defaults. Full primary-study extraction
of the old table is incomplete. Unknown reference group, adjustment, latent
definition, age, treatment and horizon remain null with reasons; this supports
withholding approval, not claiming a complete contemporary review.

The94.31% one-year RR74 result is a scientific warning, not merely numerical
boundedness. Explicit alternatives: legacy RR-as-hazard sensitivity; matching
factor-specific absolute progression hazards; contemporary conditional calibrated
HRs; a capped sensitivity only with independent evidence; exclusion pending
review. No arbitrary cap is supplied. Exclusion blocks inference rather than
turning exposed people into healthy RR1 people. Explicit omission scenarios are
mathematical sensitivities, not evidence that a factor is absent.

Concentration ranks age/time/RR cells and reports top1%,5%,10% case shares, with
each effective RR's population/rate. Top1% producing>50% of cases is a declared
diagnostic warning, not a clinical cutoff. Unreviewed effect use blocks the gate
regardless. Definition mismatch fails; union association preserves TBI, factor
and recency marginals. Recent infection is never an ordinary factor.

## Mortality and observations

Existing `data/mortality.csv` has calendar-year/age/probability without a reviewed
country/sex provenance contract; `sae_mortality.csv` and `tb_mortality.csv` concern
other outcomes. None is imported. Frozen ABS/WPP population counts are not survival.
Candidate sources: [ABS2022–24](https://www.abs.gov.au/statistics/people/population/life-expectancy/2022-2024),
[WHO methods](https://www.who.int/publications/m/item/who-methods-and-data-sources-for-life-tables-2000-2021),
[UN WPP2024 methodology](https://www.un.org/development/desa/pd/node/4326).
No snapshot is adopted. Minimum: country/subgroup, age/sex coverage through the
horizon, conditional survival, vintage, death-cause scope, licensing review,
immutable offline source and SHA-256, explicit period/cohort assumptions.
All-cause survival needs a review of TB deaths to avoid double counting.
No Australian mortality is transferred to Viet Nam; no survival is inferred
from age structure. No-mortality is diagnostic and overstates long-horizon risk.

The operator predicts first active-TB episodes after baseline among represented
disease-free residents. [WHO routine definitions](https://data.who.int/dashboards/tuberculosis/routine-tb-surveillance-indicators)
include new/recurrent episodes and distinguish re-registration categories.
Treatment-failure notifications must not simply be added as new incidence.
Aggregate estimates do not identify the resident exposure history. Notifications
need ascertainment models; prevalent/screen-detected TB is not incident TB;
imported active disease is not locally generated progression. The observation
contract preserves the original and any explicit reviewed adjustment; it never
invents an adjustment or clips a target.

## Migration, prior TB and dynamics

The [Australian Immunisation Handbook](https://immunisationhandbook.health.gov.au/contents/vaccine-preventable-diseases/tuberculosis)
reports most Australian TB among overseas-born people. Birthplace is not proof
of imported active disease or infection location. National incidence alone
cannot reconstruct subgroup exposure/arrival/treatment histories. Prefer user
TBI and explicit subgroup forward validation; no empirical subgroup values are
invented here.

History metadata includes exposure geography, arrival-year policy, imported
active contribution, prior-TB prevalence, recurrence contribution and historical
depletion policy. Required missing values remain missing. No migration/prior-TB
model is invented. The3B fictional depletion sensitivity removes680.264223 from
10,000 and changes one-year cases234.920504 to212.505883. A3C20% imported/recurrent
target adjustment is explicitly hypothetical. Actual contributions and historical
TB mortality/selective survival remain unidentified.

In high transmission, incidence and pressure are jointly generated. Fixed
pressure cannot represent feedback or intervention effects. Inversion could
only initialize/inform a dynamic model, not replace Starsim or another engine.
