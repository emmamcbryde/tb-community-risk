# Milestone 3A decision register

These decisions distinguish software policy approval from scientific approval
of effect values.

## Approved allocation policies

- User-entered age/sex factor prevalences; zero, missing, excluded, unreviewed and
  incompatible values remain distinct.
- Explicit allocation and reviewed definition mappings with provenance; no
  automatic exposure transport or interpolation.
- Adult evidence is not automatically applied below 15. Recent infection remains
  derived by the catalytic model, not an ordinary risk factor.
- Conditional independence within age/sex/subgroup cells as provisional overlap;
  complete reviewed joint data can replace it.
- Highest applicable reviewed RR centrally, all ties retained; no addition,
  multiplication or extra multimorbidity multiplier centrally.
- Normalized feasible TBI–any-risk association within age/sex/subgroup cells;
  achieved phi separately, factor/TBI marginals and source recency preserved.
- Population-count combination after subgroup allocation; incomplete groups and
  100+ remain visible, with no arbitrary terminal age.
- Explicit baseline source modes and deliberately blocked incidence inference
  pending an approved progression and risk forward operator.

## Requires later scientific review

- Contemporary progression effects and transport to current populations,
  treatment eras and the eventual progression operator.
- Pediatric effects; no scaled adult substitutes.
- Factor-definition mappings, especially smoking intensity, ideal weight versus
  BMI, CKD versus haemodialysis and cancer site.
- Factor-specific TBI association, joint multimorbidity, state-specific
  association, treatment-modified effects and migration-related profiles.
- Ascertainment, competing mortality, imported disease and selective survival
  assumptions required before active-TB inference is identifiable.

The audited 2008 transcription remains provisional. The default effect policy is
reviewed-only; legacy numeric demonstrations require explicit diagnostic use and
cannot become downstream clinical inputs automatically.

## Deferred

- Active-TB forward operator, incidence fitting and incidence-based TBI inference.
- Individual stochastic generator (explicit RNG/validation contract recorded),
  migration model and endogenous transmission.
- Interventions, screening/treatment, economics, DALYs/QALYs and production UI.

Recommended next milestone: **the age-, time-since-infection- and risk-specific
active-TB progression forward operator**. It needs reviewed natural history and
effects, explicit use of the most-recent-reinfection clock policy, and separate
disease-observation/ascertainment validation. It has not been started here.
