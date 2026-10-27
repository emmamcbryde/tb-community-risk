# Staged country TB implementation plan

## Milestone 1 architecture and evidence

Complete the seven proposed contracts, evidence registry, WHO/demographic/code audits, identifiability and association mathematics, isolated validation prototypes and protection checks. Scientific review is the exit gate; software test success alone does not approve the assumptions.

## Milestone 2 frozen data boundaries and scientific calibration design

Implement isolated offline geography/incidence/demography loaders and semantic validators against the proposed contracts. Require an explicit production WHO snapshot ID/hash; reject fixture fallback and unavailable country estimates. Archive source/dictionary files under their terms with verified checksums, inspect the full manifest identity, and retain report-round provenance.

Select and freeze a reviewed demographic source and revision. Resolve duplicate country-age cache keys, open-age support, national/subnational definitions, historical demography, mortality and migration. Do not build production UI or connect either engine.

In parallel within the scientific review process, agree on child/adult progression kernels, HIV/ART handling, infection-versus-test positivity, reinfection/clearance, recent window, historical factor profiles and at least two defensible pre-WHO scenarios. Audit modern evidence for all eleven factors and define how legacy RRs may enter an approved risk representation.

Deliver executable cross-artifact checks and small invented calibration fixtures with analytical expected behavior. If evidence decisions remain open, keep outputs explicitly exploratory and readiness false. Exit only when the data boundaries and calibration parameter-identification argument are reviewable.

## Milestone 3 infection history and calibration

Implement the cohort history and observation model independently of interventions. Use incidence-informed shape and explicit natural history; fit the minimal scale model, then optional observed TBI anchors. Test recovery on synthetic data, incompatible anchors, age/sex aggregation, historical coverage, boundary hits and structural sensitivities. Display profiles and diagnostics through a research review artifact, not production pages.

Require a scientific review of plausible fitted histories and uncertainty. Do not fix an inconsistent age profile by freely fitting a parameter for every band.

## Milestone 4 common synthetic baseline

Generate one reproducible population per calibration draw using age-specific factor marginals, declared union allocation and feasible residual association. Compare expected and realized marginals and multimorbidity. Finalize mortality, prior/active TB, infectiousness and eligibility annotations; hash all parents and seed streams.

Demonstrate two independent adapters read the identical baseline bytes/identity without regenerating it. A dynamic adapter cannot claim readiness until its initial disease-state requirements are met.

## Milestone 5 projection interfaces

Implement the static/direct-effects engine boundary first, with explicit future exposure mode. Validate person-time, competing events and direct outcome definitions. Separately implement a dynamic transmission adapter only with reviewed mixing, reinfection, detection/treatment and calibration; a Starsim prototype is evidence of software feasibility, not a calibrated engine.

Compare the two engines under zero-transmission/direct-only conditions where equivalence is expected. Differences after introducing transmission must be attributable to declared processes, not different baseline populations.

## Later milestones

Add intervention pathways after baseline and projection review. Build the proposed country UI as an isolated application only after input semantics are stable. Screening, preventive treatment, diagnostics, health economics, DALYs, ICERs and deployment each need their own explicit scope and review; none are part of Milestones 1–4.

## Decisions requiring review

| Decision | Required resolution |
| --- | --- |
| Pre-2000 exposure | Country-specific backcast choice and alternatives; no silent default |
| Trend link | Elasticity, smoothing, COVID handling, discrepancy and local transport assumptions |
| Natural history | Age at/current age/time since infection, survivor selection, clearance and reinfection |
| Pediatric progression | Infant and 5–14 evidence, adolescent boundary and factor-specific permissions |
| Medical effects | Contemporary RR/HR interpretation, HIV/ART, factor overlap and reference-risk standardization |
| Multimorbidity | Provisional highest policy versus supported interactions; validate union patterns |
| TBI anchors | Survey definition/representativeness, assay mapping, covariance and incompatibility rules |
| Demography | Source revision, duplicates, oldest ages, sex, migration and mortality |
| Disease initialization | Active/prior TB and infectiousness; do not infer prevalent disease from incidence alone |
| Source rights | WHO subset redistribution clarification and rights to any third-party demographic inputs |

Recommended next step is Milestone 2: production snapshot guards, frozen demography and semantic validation, alongside the scientific calibration design. It is not a request to wire screening or launch a production app.
