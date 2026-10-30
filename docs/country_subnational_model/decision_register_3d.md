# Milestone 3D decision register

Date2026-10-30. Internal synthesis only; no automatic central approval or external
sign-off. [Reviewer package](milestone3d_review.md) contains response questions.

## Evidence-supported candidates

- Exact licensed age/sex/assay observations from Viet Nam contacts, NHANES and
  Australian premigration applications; source observations and transformations
  are separate, immutable and checksummed.
- Adult source four-point risks and contemporary factor associations are
  reviewable evidence, with effect measures and reference populations retained.
  They do not automatically become country cause-specific hazard parameters.
- ABS2022–2024Australia age/sex period life table is a reproducible licensed
  validation input. It is not adopted as the central competing-mortality model.
- Product-specific adult IGRA performance supports explicitly fixed diagnostic
  scenarios, pending population/BCG/immunosuppression transport review.

## Provisional modelling assumptions

- Three source-knot scenario sets, linear cumulative-hazard interpolation and
  explicit illustrative fixed tails; no pediatric substitution or default set.
- Explicit uniform ages within finite diagnostic survey bands, never midpoints;
  open ages and crossing child/adult bands stay unresolved.
- Contact sample-count weights and equal NHANES age-record weights are
  descriptive metrics, not national population weights or design-correct tests.
- The existing constant pre2000 pressure backcast remains the diagnostic central
  history policy; bounded trend is a named sensitivity, not residual-tuned.
- Legacy RR-as-hazard remains diagnostic. Highest applicable RR and residual
  within-age/sex TBI–any-risk association mathematics are unchanged.

## External review required

- A transportable adult net progression curve;0–4 and5–14 curves; early disease
  exclusion and infection timing; late tail and age dependence.
- Factor-specific exposure definitions, contemporary conditional estimates,
  treatment modifiers and whether any effect may act as a hazard multiplier.
  Broad granuloma definition conflicts with fully calcified-lesion guidance.
- Competing-death cause scope, national/subgroup mortality and open-age policy.
- Assay transport, viable untreated TBI versus persistent positivity, prior-TB
  and recurrence/importation adjustments, migration/exposure histories.
- Study-specific acceptance thresholds, exclusion thresholds, metric weights,
  independent validation splits and external review responses.

No inference configuration is approved by the phrase “central candidate”. The
reported fictional Viet Nam-shaped65.33% TBI remains fictional, not empirical.

## Data required

- Exact licensed Ca Mau community age-specific TBI data; compatible Viet Nam
  risk-factor, life-table, prior-TB/episode and resident exposure inputs.
- Australian-born, recent/long-term arrival and clinical-group TBI records with
  assay, eligibility, age/sex, arrival/residence and exposure details. The mixed
  premigration child table alone cannot supply these.
- Pediatric disease-free incident infection-time data; late remote risks;
  factor-specific TBI-positive cohorts; country/subgroup survival.
- Independent compatible incidence outcomes and defensible recurrence/import
  components. No fractions are fabricated or subtracted.

## Dynamic-model work

High-transmission incidence and pressure are endogenous. This static layer may
provide bounded initialization/diagnostic inputs only. Source/parameter hashes,
subgroups, latent/observed distinctions and missing components must accompany
a later Summer/Starsim handoff. No dynamic integration is undertaken here.

## Remains available / remains blocked

User overall/age/sex TBI anchors, direct total/recent/remote profiles and
population-weighted subgroups remain available without incidence recalibration.
The new validation utilities preserve those inputs and report disagreement.
Empirical country inversion remains blocked, with no production exposure.
The public status remains
`unavailable_pending_progression_and_risk_forward_operator`.

The three3C prospective refinement failures are diagnosed as numerical-resolution
limitations in the checked fixtures. No3B/3C core defect correction or relaxed
tolerance is required. Original failed fits remain unadopted. The independent
reference is a limited numerical audit, not a new clinical forward model.

Next milestone: obtain scientific reviewer responses and targeted matched
datasets, then lock a prospective held-out validation plan. Production UI,
inverse country approval, interventions, economics, DALYs/QALYs and endogenous
transmission remain deferred.
