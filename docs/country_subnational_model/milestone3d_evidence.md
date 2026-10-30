# Milestone 3D scientific evidence consolidation

Search and extraction date: **2026-10-30**. Read the prospectively recorded
[protocol](milestone3d_protocol.md), [reviewer package](milestone3d_review.md)
and [empirical results](milestone3d_empirical.md) together. This is a focused
internal synthesis, not an exhaustive systematic review or external validation.
Numerical extraction/consistency review does not approve transport to countries.
The executable registry retains included, excluded, overlapping and data-required
sources; it never updates the frozen 3A factor catalogue.

## Adult natural history

[Menzies NA, Swartwood N, Testa C, et al. Time Since Infection and Risks of Future
Disease for Individuals with Mycobacterium tuberculosis Infection in the United
States. Epidemiology. 2021;32:70–78.](https://doi.org/10.1097/EDE.0000000000001271)
provides model outputs, with historical cohort inputs, rather than a new national
infection-conversion cohort. The focal example is infection at age20 without
other risk factors. Source inputs include the UK MRC15-year, US placebo10-year
and Netherlands contact15-year follow-up evidence. Contact ascertainment and
test conversion do not locate biological infection exactly. Baseline screening,
early/prevalent disease, source treatment and survivor selection matter.
Remaining-lifetime estimates include US mortality and age assumptions; their
transport to a country-specific cause-specific hazard is not established.

Four published focal risks support a **sensitivity-only surrogate**:
1year3.8%,5years6.6%,10years7.2%,25years7.9%. The earlier2year5.0% value remains
unverified and excluded. Below,2year risk is interpolation, not a new observation.
Neither exact ages at infection other than the focal example nor very early
within-year hazard are independently identified by these four points.

| Proposed diagnostic set | F(1) | F(5) | F(10) | F(25) | Diagnostic fixed tail/year |
| --- | ---: | ---: | ---: | ---: | ---: |
| Conservative | .032 | .058 | .064 | .070 | .00032 |
| Central candidate | .038 | .066 | .072 | .079 | .00038 |
| Higher progression | .045 | .077 | .083 | .089 | .00045 |

Lower/upper knots are published **marginal** uncertainty endpoints. Connecting
them is a modelling scenario, not a simultaneous confidence band. In particular,
the higher set's10–25year segment is lower than the central segment: ordered
cumulative risks do not order every instantaneous hazard or remote-age outcome.
The tail rates use the source's approximate year25 rate and uncertainty only as
an illustrative extension; constancy beyond25 is not evidence-supported approval.
All three sets have `selected_by_default=false`, pediatric curves null, central
adult/tail null, and scientific approval pending.

The transformation is H(t)=-log(1-F(t)). Between knots cumulative hazard is
linear, so each segment hazard is [H(t2)-H(t1)]/(t2-t1). This mathematical
transformation alone does not turn a competing-risk cumulative incidence into
a transportable net hazard.

| Set | H(1) | H(5) | H(10) | H(25) | h0–1 | h1–5 | h5–10 | h10–25 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Conservative | .03252319 | .05975000 | .06613980 | .07257069 | .03252319 | .00680670 | .00127796 | .00042873 |
| Central candidate | .03874083 | .06827884 | .07472355 | .08229524 | .03874083 | .00738450 | .00128894 | .00050478 |
| Higher progression | .04604394 | .08012604 | .08664781 | .09321238 | .04604394 | .00852053 | .00130435 | .00043764 |

[Sloot R et al. Risk of tuberculosis after recent exposure: a10-year follow-up
study of contacts in Amsterdam. Am J Respir Crit Care Med.2014.](https://pubmed.ncbi.nlm.nih.gov/25265362/)
reports five-year incident disease2.4%(95%CI1.2–4.7) in untreated infected
contacts. Its contact-timing, baseline disease exclusion and treatment selection
differ from the focal model's infection-time estimand. It is an important
lower-risk challenge, not another knot to splice into that curve. No single
reviewed adult central curve emerges from this synthesis.

## Pediatric gate

[Martinez L et al. The risk of tuberculosis in children after close exposure:
a systematic review and individual-participant meta-analysis. Lancet.
2020;395:973–984.](https://doi.org/10.1016/S0140-6736(20)30166-5)
included137647 exposed children/adolescents,1299 prevalent and999 incident
outcomes. The often-quoted under5 two-year19%(8.4–37.4) includes prevalent and
incident disease. Approximately83% of under5 disease was diagnosed within90days
of baseline evaluation. This is not an untreated first-episode progression
curve starting at a known infection event in a disease-free risk set.
Infection definitions, contact intensity, BCG, preventive treatment, follow-up,
ascertainment and competing death differ among cohorts.

[da Costa FBP, Nicol MP, Botha M, et al. Mycobacterium tuberculosis infection
and tuberculosis disease in the first decade of life: a South African birth
cohort study. Lancet Child Adolesc Health.2024.](https://doi.org/10.1016/S2352-4642(24)00256-6)
adds prospective repeated-test evidence under ongoing community exposure.
BCG context, selectively offered preventive treatment and repeated infection
must be retained. The protective treatment aHR.23(.12–.47) is not an untreated
natural-history parameter. Follow-up to10years cannot supply the entire5–14
curve, and test-conversion intervals are not exact infection dates.

**Decision:**0–4 and5–14 remain external-required. Retain finer infant,1–4,5–9
and10–14 evidence groups when acquired; do not multiply adult hazards to fill
them. At the frozen2024 population, unresolved pediatric denominators are
4,794,286/27,204,809 in Australia(17.6229%) and
23,453,385/100,987,682 in Viet Nam(23.2240%). These alone preclude complete
country progression/inversion, even before adult approval is considered.
The corresponding100+ populations6,588 and7,139 remain separately visible.

## Late reactivation and remote tail

[Dale KD et al. Quantifying the rates of late reactivation tuberculosis:
a systematic review. Lancet Infect Dis.2021.](https://doi.org/10.1016/S1473-3099(20)30728-3)
supports persistent late disease but provides sparse untreated evidence with
known acquisition time beyond25years. Unknown-exposure cohorts can show
rates below100/100000PY beyond10years; reinfection, age, selective survival,
treatment and infection ascertainment prevent a universal remote hazard.
Lack of evidence is not evidence of zero progression.

The following are transparent **mathematical sensitivities**, using central
H(25)=.08229524. F(45) is cumulative risk45years after infection, not an invented
terminal age. No-mortality infinite-horizon limits are not human lifetime risks.

| Tail after25 | Hazard just after25/year | F(45) | No-mortality infinite limit | Central decision |
| --- | ---: | ---: | --- | --- |
| Constant final segment | .00050478 |8.8251%|100%|Sensitivity only |
| Fixed external lower/mid/upper illustration | .00032/.00038/.00045 |8.4876/8.5973/8.7252%|100% for each|No reviewed external national rate supplied |
| Decline to zero, half-life10years | Starts .00050478; at26 .00047098 |8.4017%|8.5683%|Half-life is an explicit assumption, not a source estimate |
| Decline to positive plateau .0001, half-life10years | Starts .00050478; at26 .00047767 |8.4857%|100%|Plateau/half-life unapproved |
| Age-dependent or subgroup-specific external hazard |External-required|Unavailable|Requires survival and reviewed rates|Data required |
| No estimate |Unresolved|Unavailable|Unavailable|Retained public scientific policy |

For a declining tail h(x)=p+(h25-p)exp(-kx), its incremental hazard overd years
is pd+(h25-p)(1-exp(-kd))/k. k=log(2)/10 above is labelled an assumption.
Uncapped constant hazards imply eventual saturation only in a biologically
unrealistic death-free infinite horizon. A mortality-weighted remaining-lifetime
risk requires a reviewed country/sex/subgroup survival function.

The retained3C **fictional** Viet Nam-shaped comparison shows conditional alpha
.01782389–.01869405 and TBI64.6287–66.2328% across its explicit positive tails;
no-tail gives no estimate. These are not national estimates or confidence
intervals. On the independently rechecked fictional population, no-tail retains
unresolved exposed population and refuses an incidence total; it never assigns
zero beyond25. Tail uncertainty particularly affects older people. A defensible
central late tail remains unresolved.

## Progression factor qualification

The exact legacy definitions and candidate RRs remain unchanged in3A. Every
row below is **not approved as a central conditional-TBI hazard multiplier**.
Newer disease associations may include acquisition, detection and confounding.
Different effect measures/reference groups must not be substituted silently.
Missing TBI definitions, adjustments or follow-up are explicit in the registry.

| Factor: exact legacy exposure; RR | Contemporary evidence and estimand | 3D status |
| --- | --- | --- |
| CXR granuloma;2 |[CDC2020 LTBI guide,p18](https://stacks.cdc.gov/view/cdc/118719/cdc_118719_DS1.pdf): fully calcified discrete granulomas do not increase progression risk; fibrotic old-TB lesions differ |Incompatible definition |
| Smoking>one pack/day;3 |[J Glob Health2026;16:04079](https://doi.org/10.7189/jogh.16.04079):29 Chinese studies, any-smoking HR2.39(1.28–4.45),OR1.77(1.29–2.43), separate measures; not heavy-smoking TBI-conditional hazard |Incompatible definition |
| Weight<90%ideal;3 |[Sinha et al. Clin Infect Dis.2024;79:233](https://pmc.ncbi.nlm.nih.gov/articles/PMC11259213/): prospective Indian household-contact BMI evidence; BMI thresholds are not ideal-weight percentages |Incompatible definition |
| Diabetes;3.6 |[Franco JVA et al. Cochrane2024,CD016013](https://doi.org/10.1002/14651858.CD016013.pub2):48cohorts, median5years; adjustedHR1.90(1.51–2.40),RR1.60(1.42–1.80),OR1.61(1.27–2.04); mostly adult general populations |Sensitivity only |
| TNF inhibitor treatment;4 |[Kwon OC et al. Rheumatology2025](https://doi.org/10.1093/rheumatology/keaf038):43643 Korean ankylosing-spondylitis adults,2010–2021; TB13.96/10000PY without biologic vs14.28 etanercept,27.79 adalimumab,33.62 infliximab; screened/treatment era, agent-specific |Sensitivity only |
| Fibronodular CXR disease;6 |[Zhu P et al. Nat Commun.2025;16:9866](https://doi.org/10.1038/s41467-025-64834-9):183808 older Chinese adults,2020–2022; stable/uncertain/active-suggestive imaging differ. HR16.7 for active-suggestive imaging is not a scar effect |Insufficient evidence for matching estimand |
| Head/neck carcinoma;16 |[Lai et al.2017 Taiwan matched cohort](https://pmc.ncbi.nlm.nih.gov/articles/PMC5671858/): pulmonary TB4.70 vs1.64/1000PY; adjustedHR2.90(2.11–3.99), not TBI-conditioned or universal across treatment eras |Sensitivity only |
| Haemodialysis;25.3 |[Dobler CC, McDonald SP, Marks GB. PLoS One.2011;6:e29563](https://doi.org/10.1371/journal.pone.0029563):14506 Australian patients,37TB, mean3.8years;66.8 vs5.7/100000PY; adjustedRR7.8(3.3–18.7), adjusts age/sex/Indigenous status/birth-country TB incidence |Diagnostic legacy; no approved new multiplier |
| Silicosis;30 |[Ehrlich R et al. BMC Public Health.2021](https://doi.org/10.1186/s12889-021-10711-1): eight studies pooled diseaseRR4.01(2.88–5.58); severity, HIV and ongoing exposure differ; silica without silicosis is another definition |Sensitivity only |
| HIV;50 |[Gatechompol S et al. BMC Infect Dis.2022](https://doi.org/10.1186/s12879-022-07332-3):2636 Thai adults onART,24229PY, median7.6years;4.66/1000PY; competing-death CIF afterART, not since infection; CD4/regimen/TPT required |Insufficient evidence for common effect |
| Transplantation with immunosuppression;74 |[TBnet. J Infect.2026;92:106668](https://doi.org/10.1016/j.jinf.2025.106668):5805 European adult recipients, median51;68.1/100000PY,SIR6.2(3.9–9.3);33.8% infection-tested,10.3% TPT; organ/regimen/acquisition heterogeneity |Sensitivity only |

Diabetes HR prediction interval.83–4.39 and follow-up-dependent estimates
challenge universal transport. The2023 nondialysis CKD review is not evidence
for haemodialysis, and LTBI prevalence meta-analyses are not progression effects.
HIV ART-time risks cannot become infection-time curves. Transplant SIRs compare
with background populations, not with untreated infected adults.

RR74 applied to the1year3.8% risk via hazard gives
1-exp[-74(-log(.962))]=94.3121%. This is a **scientific warning**. Its numerical
boundedness does not establish plausibility. Highest-RR combination remains
central mathematics, but legacy RR-as-hazard use is diagnostic only. Appropriate
future alternatives are factor-specific absolute risks or calibrated conditional
hazard ratios with disease-free baseline, TBI ascertainment, treatment and
competing death. None is invented here; no universal cap is introduced.

## Mortality sources

| Source | Coverage/resolution/vintage | Reproducibility, licence and open ages | Decision |
| --- | --- | --- | --- |
|[ABS Life expectancy2022–2024,Table9](https://www.abs.gov.au/statistics/people/population/life-expectancy/2022-2024)|Australia, sex-specific ages0–100, period life table|CC BY4; workbook checksum and exact cell locations frozen. Derived1-qx only0–99; raw100 retained, not a100+ age distribution|Small immutable validation snapshot; not a default |
|[UN WPP2024](https://population.un.org/wpp/)|International countries including Viet Nam; sex, complete/abridged tables, estimates/projections by object|[Portal CC BY3IGO](https://population.un.org/dataportalapi/index.html); object/version/estimate-vs-projection and terminal-age extraction must be checked before freezing|Preferred international candidate, snapshot required |
|[WHO life tables2000–2021 methods](https://www.who.int/publications/m/item/who-methods-and-data-sources-for-life-tables-2000-2021)|Country/sex life tables; modelling and country data vary|Release-specific licence, age completeness, period/interval semantics and checksum needed|Candidate, no downloaded values adopted |
|User-supplied subgroup table|Depends on population, treatment and residence|Require explicit ages, sex, intervals, survival not life expectancy, source/licence/hash and open-age policy|Useful subgroup alternative; absent until supplied |

ABS source SHA-256:
`b8ed69c08132a848689b4d0a4611aa876bbabef38e464fde5da244662c2e583e`.
Snapshot contains202 raw sex/age rows and200 finite one-year survival intervals.
All-cause mortality may include TB deaths; whether/how to remove them for the
operator's competing-death cause needs review. No Australian mortality is
transferred to Viet Nam. No Vietnamese/global or subgroup survival is adopted.
Instantaneous illustrative rates do not require mortality integration; long
horizons and remaining-lifetime outputs remain incomplete without it.

## Episodes, importation and prior disease

[WHO routine surveillance definitions](https://data.who.int/dashboards/tuberculosis/routine-tb-surveillance-indicators)
and the [2025 report](https://www.who.int/teams/global-programme-on-tuberculosis-and-lung-health/tb-reports/global-tuberculosis-report-2025/tb-diagnosis-and-treatment/2-3-tb-treatment-coverage-outcomes)
distinguish new/recurrent episodes from other re-registrations after failure
or loss to follow-up. “Relapse” terminology changed to recurrent; an aggregate
incidence estimate does not identify first resident episodes generated by
local exposure. Imported disease and infection acquired before arrival require
explicit residence/exposure and episode information.

| Candidate observation policy | Requirements | Appropriate use now |
| --- | --- | --- |
| First-episode-compatible target |Disease-free baseline, first episodes, represented resident exposure, matching age/sex/time/person-years and ascertainment|Fully specified synthetic closed population; empirical only when reviewed data supplied |
| Total-episode target with recurrence model |Separate prior-TB, relapse/reinfection, treatment/failure and survival processes|Deferred; not supplied for Viet Nam/Australia |
| Reviewed adjusted target |Source-preserving component data and reviewed observation transformation; no invented fractions|Potential future country/subgroup route; no adjustment approved here |
| WHO total incidence validation only |Retain source rate/bounds and episode/acquisition mismatch|Viet Nam contextual comparisons; Australia national inversion remains unsuitable |

No imported/recurrent fraction is subtracted. Australian migrant groups cannot
have their historical TBI reconstructed from national incidence alone. A
positive test may persist after treatment and cannot identify untreated viable
infection or prior-TB history. Missing arrival/exposure, imported active disease,
prior-TB prevalence and recurrence metadata remain explicit missing inputs.
