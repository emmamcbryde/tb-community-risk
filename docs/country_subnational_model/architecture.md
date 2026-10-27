# Country and subnational TB model architecture

## Purpose and milestone boundary

Build one versioned synthetic population for a country, subnational region, or defined community. Both future projection engines must consume that same population and its calibration draws. Geography, data sources, age structure, infection history and medical factors are explicit inputs rather than inherited APY assumptions.

Milestone 1 delivers specifications, audits, proposed schemas and pure validation mathematics. It does not connect screening, preventive treatment, economics, DALYs, ICERs, either projection engine, production pages or deployed entry points. MATLAB, Starsim, frozen reference artifacts and protected workstreams remain unchanged.

## Data flow

    Frozen WHO 2025 incidence + explicit local incidence replacement or scaling
                              |
    Frozen demography + historical exposure residence + backcast scenario
                              |
    Incidence-informed infection trend + natural-history assumptions
                              |
    Calibration against disease incidence and optional observed TBI by age
                              |
    Age-specific uninfected / recent / remote profile and diagnostics
                              |
    Age-specific medical-factor marginals + residual feasible association
                              |
    One immutable synthetic baseline population with hashes and draw identifier
                              |
             Static adapter                 Dynamic adapter

Keep source observations, analytical transformations, inferred latent quantities, synthetic realization and projection results in separate artifacts. Every transformation links its parents by content hash. A change to a calibration input invalidates the derived profile and population; it never modifies the original data.

## Geography and incidence history

Use WHO ISO3 identifiers for countries and areas, preserving historical entities and effective dates. Subnational geography has its own identifier and parent ISO3; custom communities need a stable custom identifier and a documented population definition, including residence, migration and sampling frame.

Incidence has three explicit modes: national WHO series, local replacement series, and national series multiplied by a documented nonnegative local rate ratio. Store the original series and the multiplier separately; allow year-specific ratios only with complete coverage and provenance. Multiplying population size changes expected case counts, not a per-100,000 rate. Local notifications are not automatically disease incidence: record ascertainment and any conversion separately.

WHO national incidence describes residents of all ages and sexes unless explicitly stratified. It cannot be fitted as though it described only screening-eligible adults or a risk-enriched clinic. Match the calibration denominator, geography, period, case definition and HIV composition. Flag transport of a national trend to a local population as an assumption, including uncertainty in that transport.

The selected frozen snapshot is documented in the [WHO audit](who_2025_audit.md). Production must require an explicit snapshot ID and report round and refuse test fixtures, unavailable estimates and silent fallback. Maintainer imports may use a network; analysis runs must not.

## Infection history and calibration

Incidence is TB disease, not force of infection. It reflects recent infection, remote infection, reinfection, age-dependent progression, medical factors, HIV, detection, treatment and population composition. WHO estimates are estimated disease burden rather than raw notifications.

Use a cohort infection-history model with infection pressure in infections per susceptible person-year. Let a smoothed, dimensionless incidence trend inform its temporal shape; fit an absolute infection-pressure scale through a forward natural-history observation model. Observed TBI by age, when available, anchors accumulated exposure. The [identifiability note](incidence_to_infection.md) specifies equations, likelihoods, reinfection and limitations.

Exposure history must cover the birth years of the oldest supported cohorts, including before 2000. Require an explicit choice among constant earliest level, continuation of the earliest fitted trend, user historical series, or reviewed country-specific scenario. No backcast is selected automatically. Store its time coverage, parameters, boundary continuity, source, reviewer and hash. Compare all defensible backcasts, trend windows, recent windows, progression kernels, migration histories and loss-of-infection assumptions.

Do not estimate historical medical-factor prevalence from contemporary inputs by stealth. Record whether historical risk profiles are reconstructed, held constant in a declared scenario, or represented by a separately calibrated effective progression kernel.

## Demography and age resolution

Use single-year ages internally where supported, with fractional age or birth date for infants and infection timing. If only five-year bands are available, retain that native information and require an explicit reviewed disaggregation rule. The open-ended oldest group needs a tail model or a declared finite support with sensitivity analysis; do not silently truncate at 80 or 89.

The minimum user-facing bands are 0–4, 5–14, 15–44, 45–64 and 65+. They separate young children, older children, working-age adults, older working-age adults and older adults. They are reporting bands, not homogeneous biological classes. Permit five-year or single-year expansion and sex strata when available. Never split a coarse source band across these boundaries without a demographic allocation assumption.

Age influences both cumulative exposure and progression after infection. Progression depends on age at infection, current age and time since infection. For 0–4, distinguish infants from older preschool children; for 5–14, retain individual ages and the approach to adolescence; from 15 use an adult curve retaining age variation. The 15-year boundary is a modelling convention requiring review for adolescents, not proof that adult evidence applies at every age from 15.

No adult TSTin3D medical algorithm is applied below 15 without a reviewed pediatric parameter set. Pediatric HIV, immunosuppression, transplantation, renal disease, diabetes and nutritional status need child-specific definitions and evidence. Heavy smoking, occupational silicosis and head/neck carcinoma are generally inappropriate adult defaults in children; use not_applicable only when justified, never as a convenient substitute for missing data. Radiographic findings in children require separate interpretation. The [evidence table](tstin3d_evidence.md) records factor-level restrictions.

Demography must include geography, scope, reference year, target size, native age and sex distribution, source vintage, checksum, transformations and review status. Historical cohort demography and mortality are additional calibration requirements; a contemporary cross-section alone does not reveal past residence or selective survival.

## TBI profile and observed anchor

Offer two eventual modes with these exact labels:

1. Model-estimated TBI profile
2. Observed TBI profile used as calibration anchor

For each internal age and displayed band, expose uninfected, recent, remote and total TBI prevalence; uncertainty and its meaning; observation source and derivation; observed value and uncertainty; fitted value; raw and standardized residual; and compatibility status. Recent plus remote equals total; uninfected plus total equals one within the defined non-active-TB denominator. Active disease is a separate clinical state, not a fourth part of this infection-only partition.

An observed profile supplies calibration targets, never a display-only replacement. Refit the infection-pressure scale first; free further parameters only with an approved identification argument and priors. Preserve the original estimates and display both profile versions. Account for test definition, sensitivity, specificity, reversion, sampling design, survey year and exclusion of active disease. Test positivity is not automatically viable infection prevalence.

For incompatible targets, return failed or partial calibration with residuals, parameter-bound hits and affected age bands. Do not clip probabilities, silently alter observations, or add unrestricted age parameters to force agreement. Offer explicit review of migration, assay mapping, backcast, exposure heterogeneity and natural history. An exploratory non-fitting scenario must retain its non-fitting status and cannot masquerade as a calibrated population.

Without an anchor, label prevalence model-estimated and uncertain. Distinguish conditional parameter uncertainty, WHO uncertainty, between-scenario structural uncertainty and Monte Carlo variation; do not turn an arbitrary sensitivity range into a confidence interval.

## Medical-factor prevalence and association

Users enter percentages for each eligible non-recency factor by age band, optionally by sex. Validate 0–100 inclusive without rounding internally. Missing, measured zero, excluded and not applicable are distinct. Every row has source, source vintage, geographic default status, review status and optional uncertainty limits. Unknown included factors block a reviewed baseline until resolved or explicitly excluded with a reason. Recency and contact RR 15 are absent from the medical-factor prevalence input.

Age-specific medical-factor and TBI profiles automatically create age-mediated population association. Add a separate within-age normalized association between total TBI and any included non-recency factor. Its values −1, 0 and +1 mean minimum feasible overlap, independence and maximum feasible overlap. It is not an unrestricted Pearson correlation. Use the [Fréchet mapping](association.md), show requested and achieved quantities, and reject infeasible joints without clipping.

Recommend a parsimonious initial factor allocation: independent factor marginals by age define the unconditioned factor-pattern distribution and its union prevalence; couple TBI to the resulting any-factor indicator. Conditional on any factor being present, use the same nonempty pattern distribution in each TBI stratum. This preserves factor marginals and aggregate association exactly in expectation, while adding no factor-specific association parameters. It does not imply independence of individual factors conditional on age and TBI. Compare observed multimorbidity if available.

Alternatives include conditional independence given age and TBI (requiring factor prevalences in both TBI strata and union reconciliation), validated pairwise constraints, Gaussian copulas or latent variables, and iterative proportional fitting. Their constraints and limitations are in the association note. Do not introduce a complex correlation model without data capable of identifying it.

## Progression risk representation

Store the eleven requested candidate RRs in an evidence registry; do not activate them as contemporary hazard ratios. The [evidence review](tstin3d_evidence.md) separates reported legacy estimates, clinical applicability, acquisition/progression interpretation and modern review needs.

The provisional central combination policy is the highest applicable reviewed medical RR, with one if none apply. This limits unsupported compounding but may underestimate synergistic risk; it is a transparent scenario, not a guaranteed conservative estimate of disease burden. Compare additive excess risk, justified capped multiplication and explicitly defined combined risks. Uncapped multiplication is a stress test only after scientific approval. A fitted multivariable model is preferable when transported and validated in the target population.

RR, HR and OR are different quantities. A cumulative-risk RR over a stated horizon may be mapped to a hazard representation only through a documented baseline survival model. If RR times baseline risk exceeds one, reject that interpretation rather than clipping. If a scenario treats a legacy RR as a hazard multiplier, label the approximation and vary it. Baseline hazards and medical effects must be calibrated consistently: do not multiply population-average progression by factors already embedded in its reference population without standardization.

Derived infection recency enters the baseline progression kernel; do not multiply it again by recent-infection RR 15. Contact can change acquisition or inform likely exposure timing, but must not add another contact RR 15 on top of derived recency. The legacy product reaching 2,916 is not permitted.

## Shared baseline population

The baseline is a hash-addressed, immutable artifact with a manifest and stable person IDs. Include geography and provenance, continuous age, sex or explicit unknown, infection state, first and most recent infection time where represented, time since infection, recent window and derived class, factor indicators and any-factor flag, approved risk representation, mortality input reference, screening eligibility status and reason, population weights, active/prior TB clinical state, and configuration/data/calibration hashes. Keep parameter-draw and population-realization seeds separate.

TBI status describes a modelled infection state. Preserve the observation model separately so immune sensitization is not silently equated with persistent viable organisms. A history model without clearance is an explicit approximation; clearance and reinfection require reviewed scenarios.

Generate the factor and infection distribution jointly with calibration if medical factors alter survival or progression. Do not independently generate a second cohort for a dynamic engine. For each uncertainty draw, both engines receive identical person IDs, weights, initial states and history. Adapters may add engine-specific fields but must preserve and report the parent baseline hash.

The infection profile and TBI/factor association are conditional on the non-active-TB population. If factor marginals describe the whole population, transform them explicitly to that denominator before coupling. Retain active disease as separate clinical rows with TBI state not_applicable_active_tb, keeping any known historical infection times. The whole-population manifest includes their weights; infection-profile weights exclude them. An unresolved clinical partition blocks analytical readiness.

Active TB at baseline, past TB, treatment history and infectiousness require additional observations/assumptions; incidence does not identify prevalent active disease. Mark these fields unresolved until a disease-state initializer is reviewed. A population with unresolved fields may support baseline review but must fail dynamic readiness. Do not seed Starsim's demonstration active prevalence or reclassify TBI individuals to create active disease. Static readiness also requires explicit handling of baseline active/prior disease.

Mortality needs age/sex rates, units, period, source and separation of background from TB-specific deaths. Eligibility is an annotation for future screening scope, with unknown permitted during this milestone; it must not remove people from the epidemiological denominator or signify denial of care.

## Projection boundaries and interface

Both future engines implement validate_baseline, initialize, project and export_results. Inputs are a baseline manifest and people table, approved natural-history configuration, horizon/time convention, parameter-draw ID and seed streams. Validation checks contract support, all hashes, geography/denominators, required clinical states, weights and readiness. Initialization returns the exact parent population hash and any explicitly documented adapter fields.

Results share person-time, incident TB events with origins, deaths, infection/reinfection events where represented, population balance, diagnostics, seed/config/data hashes and engine name/version. Static and dynamic outcomes must be distinguishable in exports; matching initial conditions do not imply matching future processes.

| Engine | Infection pressure | Interpretation and required disclosures |
| --- | --- | --- |
| Static/direct effects | No endogenous feedback; future pressure is zero or an explicitly selected exogenous schedule | Direct outcomes among modelled people. If exogenous infection is enabled, record reinfection support and approved progression after new exposure. No claimed indirect protection. |
| Dynamic transmission | Active disease contributes to infection pressure through a specified mixing/transmission model | Reinfection, changes in future transmission and indirect effects. Needs infectious disease-state initialization, mixing, detection/treatment histories and calibration; no double addition of baseline historical pressure as endogenous plus exogenous pressure. |

In the dynamic adapter, historical infection pressure is an inferred initialization history, not a permanently imposed future force. Any imported external force must be separately identified, with a documented border/external-population interpretation.

## Application and filesystem proposal

    country_tb_ui/                  eventual pages and presentation only
    engine/country_tb/              isolated contracts, calibration and baseline generation
    data/who_incidence/             existing frozen country snapshot
    data/demography/                future frozen reviewed demographic snapshots
    schemas/country_tb/v1/          proposed versioned contracts and review example
    docs/country_subnational_model/ architecture and evidence
    tests/country_tb/               focused independent validation
    abm/output/country_tb/          future generated artifacts, ignored by Git

Only the prototype, schemas, tests and documentation are added now. Ordinary labels do not mention APY or SA Health. Eventual pages: Geography and population; TB incidence history; TBI prevalence by age; Risk-factor profile; Projection settings; Intervention; Run analysis; Results; Health economics; Evidence and technical information. The latter intervention/economic pages remain future placeholders in this architecture, with no implementation in Milestone 1.

Each input page will show source, vintage, review and uncertainty next to the value. Geography changes preview conflicts rather than silently replacing user input. The TBI page will distinguish modelled and observed profiles and show fit failures. The factor editor will distinguish zero and missing visibly and report aggregate association by age and overall.

## Acceptance and scientific review

Software validation proves contract/math behavior, not scientific validity. Require evidence review of child and adult progression, historical backcasts, assay mapping, reinfection/clearance, local trend transport, medical factors and mortality before a policy-ready baseline. See the [staged plan](implementation_plan.md) for explicit gates and Milestone 2 scope.
