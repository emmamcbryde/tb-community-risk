# Risk-factor schema and effect-measure caveats

Status:
* table and validation: **implemented and validated**;
* engine semantics: **retained for backward compatibility**, awaiting review;
* conversion policy: **designed but not implemented**.

## Fields (per risk factor)

`riskFactorId`, `label`, `enabled`, `prevalence` (a ProfileValue),
`prevalenceBounds` [low, high], `effectEstimate` (a ProfileValue),
`effectBounds` [low, high], `effectMeasure` (`RR`, `HR`, `OR` or null, kept
exactly), `affectedTransition` (`progression_to_disease`, `infection` or
`other`), `evidenceYear`, `populationApplicability`, `evidenceSource`,
`reviewStatus`, `notes`, `engineKey`, `engineApplication`, `userModifiedFields`,
and a derived `userOverride`.

Editing in the interface marks the changed field user-defined and keeps the
other sources. CSV template, import and export use the same columns. A profile
with no risk factors is valid and runs without stratification.

## How the current engine applies effects

| Aspect | Current behaviour |
| --- | --- |
| Affected quantity | Early and late progression hazards (infection to disease) |
| OR | Used as a hazard multiplier without conversion; overstates the effect when outcomes are common |
| RR | Used as a hazard multiplier without conversion; close to a hazard ratio only when risks are small |
| HR | Used directly (consistent with its definition) |
| Unspecified | Used as an unspecified hazard multiplier |
| Joint effects | Multiplied, with no interaction terms; correlated factors may be double counted unless the estimates are mutually adjusted |
| Cap | None |
| Baseline | Progression hazards are recalibrated so that average risk matches the calibration target; effects redistribute risk between people |
| Infection effects | Separate bundled infection ORs (cannabis 3.65, contact 2.53, renal 2.19) multiply the cumulative infection hazard in a complementary log-log model; not editable here |
| User prevalence | Applied uniformly across age groups; demonstration prevalence is age-specific |

## Warnings (never change results)

The application warns when:
* an OR is applied as a hazard multiplier;
* the measure type is missing;
* two or more enabled factors are multiplied;
* the combined multiplier for a person with every enabled factor exceeds 20.

The demonstration profile's combined multiplier is 5 x 2 x 3 x 3.6 x 3 x 3 x 3 =
2,916.

## Conversion policy (inactive)

`engine/profiles/effect_measures.py` defines `ConversionPolicy`. Only `none` is
active. Two other policies are designed but inactive until reviewed:
* OR to RR using a baseline risk: `RR = OR / (1 - p0 + p0 * OR)`;
* RR to HR under a constant-hazard assumption.

Selecting either raises an error.

## Country risk-factor evidence: audit

* **gtbreport2025.** `inc_mort/analysis/attributable_cases.rda` holds country
  prevalence and SD, for 2024 and 199 countries, for diabetes, alcohol use,
  undernutrition and smoking. It also holds population-attributable fractions
  and attributable incidence. `rf.rda` holds attributable incident numbers.
  These are intermediate analysis objects in a repository without a licence.
* **Public WHO dataset.** It publishes attributable cases by risk factor, not
  prevalence or relative risks.
* **Conclusion.** Nothing is imported into the model. Population-attributable
  fractions must not be read as prevalence or relative risks. A suitably
  licensed prevalence source, with year and applicability, is needed before
  country risk-factor profiling.
