# Milestone 3B evidence review, 29 October 2026

This is a focused source audit for a forward-operator prototype, not a clinical
guideline or approval of country estimates. Source checking, numerical assumptions
and approval for parameter transport are distinct. No evidence value is adopted
merely because it existed in an earlier branch.

## Citations

1. Menzies NA, Swartwood N, Testa C, Malyuta Y, Hill AN, Marks SM, Cohen T,
   Salomon JA. Time Since Infection and Risks of Future Disease for Individuals
   with Mycobacterium tuberculosis Infection in the United States. *Epidemiology*
   2021;32:70–78. [doi:10.1097/EDE.0000000000001271](https://doi.org/10.1097/EDE.0000000000001271),
   [primary manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC7707158/).
2. Sloot R, Schim van der Loeff MF, Kouw PM, Borgdorff MW. Risk of tuberculosis
   after recent exposure: a 10-year follow-up study of contacts in Amsterdam.
   *Am J Respir Crit Care Med* 2014;190:1044–1052.
   [doi:10.1164/rccm.201406-1159OC](https://doi.org/10.1164/rccm.201406-1159OC).
3. Martinez L, Cords O, Horsburgh CR, Andrews JR; Pediatric TB Contact Studies
   Consortium. The risk of tuberculosis in children after close exposure:
   a systematic review and individual-participant meta-analysis. *Lancet*
   2020;395:973–984. [doi:10.1016/S0140-6736(20)30166-5](https://doi.org/10.1016/S0140-6736(20)30166-5),
   [primary IPD report](https://pmc.ncbi.nlm.nih.gov/articles/PMC7289654/).
4. Vynnycky E, Fine PEM. The natural history of tuberculosis: the implications
   of age-dependent risks of disease and the role of reinfection.
   *Epidemiol Infect* 1997;119:183–201.
   [doi:10.1017/S0950268897007917](https://doi.org/10.1017/S0950268897007917).
5. Andrews JR, Noubary F, Walensky RP, Cerda R, Losina E, Horsburgh CR. Risk
   of progression to active tuberculosis following reinfection with Mycobacterium
   tuberculosis. *Clin Infect Dis* 2012;54:784–791.
   [doi:10.1093/cid/cir951](https://doi.org/10.1093/cid/cir951).
6. Ekramnia M, Li Y, Haddad MB, et al. Estimated rates of progression to
   tuberculosis disease for persons infected with Mycobacterium tuberculosis
   in the United States. *Epidemiology* 2024;35:164–173.
   [doi:10.1097/EDE.0000000000001707](https://doi.org/10.1097/EDE.0000000000001707),
   [CDC primary manuscript](https://stacks.cdc.gov/view/cdc/147670).
7. Vega V, Rodríguez S, Van der Stuyft P, Seas C, Otero L. Recurrent TB: a
   systematic review and meta-analysis of the incidence rates and the proportions
   of relapses and reinfections. *Thorax* 2021;76:494–502.
   [doi:10.1136/thoraxjnl-2020-215449](https://doi.org/10.1136/thoraxjnl-2020-215449).
8. Menzies D, Adjobimey M, Ruslami R, et al. Four months of rifampin or nine
   months of isoniazid for latent tuberculosis in adults. *N Engl J Med*
   2018;379:440–453. [doi:10.1056/NEJMoa1714283](https://doi.org/10.1056/NEJMoa1714283).
9. Austin PC, Lee DS, Fine JP. Introduction to the analysis of survival data
   in the presence of competing risks. *Circulation* 2016;133:601–609.
   [doi:10.1161/CIRCULATIONAHA.115.017719](https://doi.org/10.1161/CIRCULATIONAHA.115.017719).
10. WHO. *Global tuberculosis report 2025*, technical appendix, definitions.
    [Primary appendix](https://cdn.who.int/media/docs/default-source/global-tuberculosis-report-2025/technical-appendix-gtb-report-2025.pdf).
11. ABS. *Life expectancy methodology, 2022–2024*.
    [Official methodology](https://www.abs.gov.au/methodologies/life-expectancy-methodology/2022-2024).

## Evidence-to-model decisions

| Citation | Population, age, setting and follow-up | Outcome, effect measure and uncertainty | Limitation, proposed use and review status |
| --- | --- | --- | --- |
| 1 | US evidence-synthesis model; focal infection at age 20 without added progression factors; 25-year trajectory | Cumulative risks at 1/5/10/25 years: 3.8/6.6/7.2/7.9%; reported 95% UIs: 3.2–4.5/5.8–7.7/6.4–8.3/7.0–8.9% | Model outputs, not direct contact-cohort risks. No reinfection or preventive treatment in the focal curve. US mortality explicitly enters remaining-lifetime results; transportable cause-specific hazards are not established here. Source checked; four-knot **diagnostic surrogate only**. The 5.0% two-year knot is unverified and omitted. |
| 2 | Amsterdam pulmonary-TB contact records, mixed ages, 2002–2011, up to ten years | Reported five-year incident disease risk among LTBI contacts not starting prevention: 2.4% (95% CI 1.2–4.7%) | Exposure/enrolment time is not known infection time; incident and prevalent disease differ. Use to review ascertainment and transport, not to splice a new universal curve. Source checked; no numerical default adopted. |
| 3 | 137,647 exposed children under 19, 46 cohorts in 34 countries; two-year risks and longer incident follow-up | Infected untreated under-5 group: two-year risk 19.0% (95% CI 8.4–37.4%). This cumulative outcome includes prevalent and incident disease; 83% of all under-5 cases were diagnosed within 90 days | Does not provide an unqualified post-baseline disease-free, known-infection-time hazard for 0–4 or 5–14. Preserve age groups separately, require reviewed external curves. Source checked; pediatric numerical defaults **not adopted**. |
| 4 | Age-structured reconstruction of England/Wales since 1900; historical population and reinfection processes | Model-estimated age-dependent disease mechanisms; not one contemporary prospective survival curve | Supports explicit age and reinfection structure. Historical context, identifying assumptions and population transport limit parameter reuse. Source checked for structure; no parameters imported. |
| 5 | 23 paired historical cohorts, 19,886 participants, largely occupational/adult exposure; variable cohort follow-up | Adjusted disease IRR after reinfection versus initial infection 0.21 (95% CI 0.14–0.30) | Acquisition protection, progression protection and host selection are not separated by this IRR. It does not identify clock reset or justify multiplying infection pressure by 0.21. Source checked; sensitivity rationale only. |
| 6 | US surveillance, IGRA prevalence and population estimates, principally 2011–2012; population person-years rather than known individual infection times | Stratum-specific inferred reactivation incidence rates with uncertainty from infection prevalence and ascertainment | Contemporary publication does not make old exposure definitions or RR-to-hazard transport automatic. Supports later external reconciliation; no effect estimate or fixed remote hazard adopted. |
| 7 | Successfully treated TB cohorts, 145 studies; mean follow-up about 2.3 years | Recurrence 2.26/100 person-years (95% CI 1.87–2.73), with substantial heterogeneity | Recurrence mixes relapse/reinfection and differs by setting. Excluded from this absorbing first-post-baseline-event operator; requires its own later model. Source checked, no recurrence parameter imported. |
| 8 | Adult LTBI trial in nine countries; 28-month confirmed-TB endpoint | Noninferiority comparison of preventive regimens, with confidence intervals; not an untreated natural-history curve | Treatment completion and residual disease risks are intervention-specific. Future evidence register only; no preventive-treatment effect, screening or treatment simulation added. |
| 9 | Statistical methods with clinical examples; follow-up and risks depend on the example | Cause-specific hazards versus cumulative incidence under competing events | Use cause-specific progression and death hazards with absorbing risk-set depletion. Do not transform a competing-risk CIF into a cause-specific hazard without further evidence. Methods checked; equations independently tested. |
| 10 | WHO annual population surveillance/estimation, all ages and sexes | Incidence includes new and recurrent episodes; source uncertainty retained | First-episode cohort predictions exclude additional WHO processes. Comparison is validation-only, never an inverse fitting target in 3B. Official definitions checked. |
| 11 | Australian period age/sex life-table methodology | Period survival; not a cohort forecast or a crude life-expectancy replacement for age-specific survival | Candidate source for a future immutable import only. Methodology checked; no life-table values downloaded or adopted. Not applicable automatically to Viet Nam. |

The baseline curve is therefore not based on a secondary website. Its primary
source XML was retrieved through Europe PMC after browser access was unavailable;
SHA-256 is `167a755a17de17b641222adb9513814a1af72bbf9e0cbeefa30279bd9c558f3d`.
The unmodified research copy remains under the ignored audit output, under the
article's CC BY-NC-ND terms. This package records derived assumptions and citations,
not a redistributed or modified research manuscript. Historical trial raw data and
the authors' full posterior model were not independently reproduced.

First-episode-only counting and baseline exclusion in every historical component
cohort were not independently verified. The focal model trajectory is not taken
as proof of compatibility with this first-post-baseline endpoint. This is another
reason the unqualified candidate is rejected.

## Adult and pediatric decisions

The requested five-knot adult candidate is **rejected for unqualified adoption**:
the two-year point was not verified, and cohort/competing-risk transport cannot be
assumed. `audited_adult_candidate()` exposes the verified four-point *diagnostic
surrogate*, never a clinical/country default. External reviewed curves must declare
net risk without competing mortality before conversion to cause-specific hazard.

The provisional central tail policy is **no extrapolation**. Constant final-slope,
declining-to-plateau and fixed-remote tails are explicit sensitivities. Source
uncertainty, curve interpolation and extrapolation uncertainty are separate.

Both pediatric groups remain `external_required`. No adult curve is assigned to
them automatically. Synthetic child curves occur only in mathematical tests.

## Mortality inventory

The installed frozen ABS/WPP files are demographic stocks; the WHO snapshot is
incidence data. No suitable reviewed immutable national life-table snapshot was
found. Existing `mortality.csv` and `data/mortality.csv` have columns
`calendar_year,Age,Prob`, without the required reviewed country/sex/life-table
provenance and licensing contract. `tb_mortality.csv` and `sae_mortality.csv` concern
different outcomes. They are not adopted or modified. Their headers and checksums
are recorded in `abm/output/country_tb_m3b_audit_20261029/mortality_inventory.json`.

External age/sex/subgroup one-year survival is supported. A national life-table
contract additionally requires reviewed licensing terms/URL and an immutable
snapshot SHA-256. No new mortality source was downloaded. No-mortality outputs are
mathematical diagnostics, labelled `competing_mortality_not_modelled`, and generally
overstate long-horizon TB risk and disease-free survival relative to explicit
competing death under the same progression hazard.
