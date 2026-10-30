# Milestone 3D evidence and validation protocol

Protocol recorded before the 3D evidence search on 2026-10-30. This is an
internal scientific synthesis and diagnostic validation, not external sign-off.
Production incidence inference remains blocked. Amendments must state their
reason and must not depend on whether a source makes the model fit better.

## Questions and searches

1. Which untreated, disease-free populations identify progression by exact age
   and time since infection, rather than time since contact ascertainment?
2. What late progression, competing mortality and effect estimands can be
   transported to the proposed resident first-episode operator?
3. Which observed TBI datasets support age/sex/subgroup validation, and which
   require assay, survey-design, time or migration adjustments?
4. Do the 3C prospective failures reflect quadrature, time resolution, objective
   shape or an invalid scenario? Diagnose before adjusting any tolerance.

Search date:2026-10-30. Sources: WHO guidance/technical definitions and GHO;
PubMed/Europe PMC; original journal articles and supplements; national health
surveys; ABS, UN Population Division and public national data custodians.
Priority: authoritative guidance; systematic reviews/IPD meta-analyses;
prospective cohorts; national surveys; gap-filling primary studies. Reviews map
evidence and overlapping cohorts; they do not create additional independent
observations when their primary cohorts are also used.

Search families (record actual queries in the source ledger):

- tuberculosis progression time since infection conversion adult late reactivation;
- child tuberculosis contact cohort age infection BCG preventive therapy;
- latent tuberculosis prevalence age sex Vietnam Viet Nam IGRA TST survey;
- Australia migrant refugee arrival latent tuberculosis IGRA prevalence age;
- NHANES2011–2012 tuberculosis infection age foreign born IGRA;
- each of granuloma, fibronodular, smoking intensity, ideal weight/BMI, diabetes
  control, TNF agent, head/neck cancer, dialysis, silicosis severity, HIV/ART/CD4,
  transplant organ/regimen with tuberculosis cohort adjusted hazard absolute risk;
- ABS life tables2022–2024; UN WPP2024 life tables; WHO GHO life tables licence;
- WHO new recurrent relapse treatment failure incidence imported migration.

No source is selected because of agreement with a model prediction. Seek
Viet Nam and Australian subgroups and at least one additional age-specific
setting; retain unsuitable candidates and reasons in the ledger.

## Eligibility and extraction

Eligible designs: prospective infection/contact cohorts with active disease
exclusion and follow-up; national or defined-population cross-sectional TBI
surveys; informative retrospective registry cohorts with explicit confounding
and ascertainment; systematic/IPD reviews with identifiable included studies.
Extract all relevant ages, distinguishing infants,1–4,5–9,10–14,15+ where
available. Do not transport adult effects to children. Retain sex, calendar
period, geography, residence/birth/arrival and sampling frame.

TBI ascertainment must specify TST dose/threshold or IGRA generation/threshold,
BCG/immunosuppression where known, indeterminate handling and whether positivity
or adjusted latent infection is reported. Disease outcomes must distinguish
prevalent/early screen-detected disease, incident disease, first/recurrent
episodes, treatment failure and imported disease. Follow-up starts, infection
timing intervals, preventive treatment, ART and competing death must be explicit.
Missing details remain missing and restrict use.

Extract citation/DOI/URL/access date/licence; population and exclusions; ages,
sexes and time; numerator/denominator when reported; weighted prevalence and
design-based uncertainty separately from raw sample counts; assay; outcome;
follow-up; reference group; adjusted effect type/estimate/interval/covariates;
absolute rate; transformations; applicability; bias; source-study overlap;
proposed use and review status. Do not reconstruct counts from weighted or
rounded prevalence, or call a disease RR a conditional-TBI hazard ratio.

Exclude numerical adoption when definitions, denominator, provenance/licence,
baseline disease exclusion or infection time are inadequate for the intended
estimand. Such sources may remain contextual or sensitivity evidence. Assess
selection, survivor depletion, treatment selection, reverse causation, residual
confounding, assay misclassification, publication bias and transportability.

## Data and analysis rules

Immutable snapshots preserve source observations separately from transformed
inputs, with source licence, extraction location, byte SHA-256 and scientific
hash. Use exact published values or licensed machine-readable data; no unlogged
figure digitization. Reimport is identical-or-fail, never overwrite. Do not
adopt mortality or assay performance silently.

Predeclare a small set: conservative, central candidate, higher-progression,
plus named one-way assay, tail, backcast, mortality and observation sensitivities.
Published marginal uncertainty bounds are scenarios, not joint curve confidence
bands. Unsupported pediatric/remote components stay external-required. Every
set is proposed for review, never automatically approved.

Direction A: TBI anchors fit only TBI, then independently compare incidence.
Direction B: only a compatible incidence target passing every3C gate may fit
alpha; independent TBI observations then validate, never fit. Otherwise emit a
structured blocked result. Cross-validation splits are by independent cohort
or declared nonoverlapping age records; do not duplicate parent totals, reuse
the same observations, or treat correlated assay results as independent.

Evaluate at the survey year when supported; do not call a 2011 survey a 2024
validation. Preserve open ages, unresolved denominators and subgroup exposure
history. No display-band midpoints or silent within-band interpolation. Any
uniform-within-band approximation is a labelled diagnostic assumption.

Metrics: observed/predicted positivity, absolute/relative error, ratio, age/sex
residuals, population/design-appropriate weighted RMSE, interval inclusion,
incidence discrepancy, unresolved population and source compatibility. Distinguish
source uncertainty, assay/parameter/scenario variation, numerical error and
model discrepancy. Never force discrepancy into alpha or tune assay and alpha
jointly from one observation.

No acceptance threshold is selected centrally. Present candidate absolute-error,
RMSE, coverage and unresolved-population thresholds for reviewer approval;
until approved, report descriptive results and no pass/fail scientific claim.
External expert response is separately recorded as pending. Static historical
pressure cannot substitute for endogenous transmission in Summer/Starsim.

## Prespecified extraction-dependent analysis choices

Recorded after source/design inspection and before computing empirical fits:

- Viet Nam2020 and2021 household-contact QFT-Plus: fit the30–39 band separately;
  validate other disjoint age records. Also fit2020 and compare2021 as a
  district/cohort transport stress test. Keep household clustering explicit.
- US NHANES2011–2012 QFT-GIT: fit US-born15–24; validate other US-born ages.
  Other nativity groups remain descriptive because exposure/migration history
  is missing. No survey counts are inferred from published weighted prevalence.
- Australia premigration child applications: preserve observed ages/tests;
  mixed assay rows, ambiguous published<4 label, selection and exposure history
  prevent latent-profile/progression validation. Do not replace with modelled
  Australian prevalence or extrapolate adult assay performance to children.
- Evaluate finite adult bands using explicitly labelled uniform-within-band
  age distributions; open ages and unreported joint age/sex remain unresolved.
- Metric weights: actual tested denominators for the convenience-contact
  descriptive comparison; equal age-record weights for published NHANES
  estimates. Neither is a national population-weighted performance claim.
- Fixed product-specific USPSTF assay candidates and named source-bound
  one-way sensitivities; no performance fitting. Three adult Menzies-derived
  curve scenarios affect forward outcomes, not TBI-anchored alpha.
- Same-year WHO incidence can be contextual only when episode/population
  compatibility fails. DirectionB is blocked for these empirical targets;
  use a fully specified synthetic fixture to test its implemented mechanics.
