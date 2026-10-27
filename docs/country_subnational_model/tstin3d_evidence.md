# TSTin3D candidate progression evidence

## Scope and source interpretation

The eleven candidate values below are the requested legacy inputs. They are not validated contemporary hazard ratios. Every numerical effect remains sensitivity-only until its definition, target population, reference risk, confounding, treatment era and age applicability are reviewed. No medical-factor engine is connected in this milestone.

[TSTin3D calculator version 3](https://www.tstin3d.com/en/calc.html) describes an adult tool. Its form includes broader underweight criteria (BMI at most 20 as an alternative), diabetes of any type, and TNF inhibitors beyond infliximab. These broadenings need definition matching before using the corresponding legacy estimates.

[Menzies et al. 2008](https://pubmed.ncbi.nlm.nih.gov/18419884/) developed an adult positive-TST interpretation aid from published summaries. M08 below means its [Table 2 and numbered reference chain](https://www.tstin3d.com/thinking3D_2008.pdf). Original studies listed through that chain have not been independently reappraised here. The table is an evidence inventory, not a modern systematic review.

## Candidate factor inventory

P means a candidate modifier of progression after infection. B means acquisition/exposure and progression associations may both contribute to observed disease risk; the reported RR does not separate them. These classifications are proposed causal interpretations for review, not additional numerical acquisition multipliers. A radiographic marker can indicate prior disease or infection history rather than a causal exposure.

| Stable factor ID | Definition to record | Candidate RR | Evidence citation in M08 | Population and age applicability | Acquisition or progression interpretation | Analysis status and required modern review |
| --- | --- | ---: | --- | --- | --- | --- |
| cxr_granuloma | Granulomatous abnormality on chest radiography | 2.0 | refs 24–25; Horwitz 1969, Grzybowski 1971 | Adult radiographic cohorts; pediatric morphology needs separate review | P marker of past disease/exposure; not an acquisition coefficient | Sensitivity; central pending imaging definition and treated/untreated history |
| heavy_smoking | More than one pack of cigarettes daily | 3.0 | ref 26; Maurya 2002 | Adult smokers; duration/intensity and former smoking matter; no pediatric default | B; exposure confounding and progression cannot be disentangled by this value | Sensitivity; central pending contemporary intensity-specific review |
| underweight | Weight below 90% of ideal; retain criterion separately from BMI | 3.0 | ref 27; Comstock 1975 | Adult anthropometry; do not translate to child BMI or growth faltering | B; reverse causation from undiagnosed TB and deprivation must be assessed | Sensitivity; central pending nutritional definition, temporality and pediatric growth evidence |
| diabetes | Clinically defined diabetes; specify type, control and duration | 3.6 | refs 28–31; Kim 1995, Silwer 1958, Pablos-Mendez 1997, Boucot 1957 | Predominantly adult disease cohorts; child type 1 diabetes needs separate evidence | P candidate; disease incidence associations can include acquisition differences | Sensitivity; central pending modern adjusted estimates and glycaemic control |
| tnf_inhibitor | TNF inhibitor therapy with drug, dose and duration | 4.0 | ref 32; Keane 2001 | Original infliximab context; other drugs and children cannot inherit it automatically | P; new acquisition possible during therapy but no acquisition RR specified | Sensitivity; central pending agent-specific effects, screening and prophylaxis era |
| cxr_fibronodular | Fibrotic/nodular radiographic disease consistent with past TB | 6.0 | refs 25, 33–34; Grzybowski 1971/1966, Nolan 1988 | Adults with inactive disease/refugee cohorts; exclude active disease clinically | P marker; may encode prior disease and recency | Sensitivity; central pending imaging criteria, past treatment and overlap with granuloma |
| head_neck_carcinoma | Malignancy of head or neck, with treatment status | 16 | ref 37; Rieder 1989 | Adult malignancy populations; no ordinary childhood prevalence default | P candidate; smoking, nutrition and treatment confounding | Sensitivity; central pending cancer site, stage and modern therapy review |
| haemodialysis | Renal failure requiring maintenance haemodialysis | 25.3 | refs 38–41; Malhotra 1981, Lundin 1979, Andrew 1980, Pradhan 1974 | Dialysis cohorts; not all chronic kidney disease; child evidence separate | B; immune impairment and healthcare exposure may coexist | Sensitivity; central pending dialysis setting, comparator and modern infection control |
| silicosis | Diagnosed pulmonary silicosis | 30 | refs 42–43; Hong Kong trial 1992, Cowie 1994 | Occupational/mining cohorts; ordinarily unavailable as a child default | B; occupational exposure and progression both plausible | Sensitivity; central pending severity, dust exposure and HIV adjustment |
| hiv | HIV infection with stage, CD4 and ART status recorded when known | 50 | refs 44–47; Antonucci 1995, Guelar 1993, Selwyn 1989, Wood 2000 | Adult HIV cohorts across older treatment eras; pediatric HIV separate | P dominant; acquisition differences cannot be inferred from this RR | Sensitivity; central requires ART/CD4-specific contemporary review and consistency with TB/HIV calibration |
| transplant_immunosuppression | Transplantation requiring immunosuppressive treatment | 74 | refs 48–51; Sakhuja 1996, Aguado 1997, Miller 1995, Meyers 1994 | Adult organ-specific transplant cohorts; children and stem-cell transplants require review | B; donor-derived/acquired infection and progression can coexist | Sensitivity; central pending organ, regimen, time after transplant and prophylaxis |

All eleven require modern evidence review. A numerical value for a historical adult disease association must not become an acquisition hazard, a pediatric progression coefficient or an adjusted causal HR merely because it appears in this registry.

## Multiple simultaneous factors

The paper's two-factor example is consistent with adding the two RRs (2 + 3 = 5), not multiplying them. More directly, read-only inspection of the current [calculator script](https://www.tstin3d.com/en/script.js), accessed 2026-10-27, found the comment that multiple selected factor values are summed and the calculation r += selected_value. Script SHA-256: 35fd9e21360039ed09cef8db4e1de1f0135c75c4268f24a052b2147bc52af05b. No script code is imported or redistributed.

That website behavior does not validate summation scientifically. Summing raw RRs also differs from adding excess risks: for two factors with RRs 2 and 3, raw summation gives 5, while 1 + (2−1) + (3−1) gives 4.

| Policy | Definition | Interpretation, bias and identification |
| --- | --- | --- |
| Highest applicable | max(1, RR_j) | Parsimonious provisional central policy after effect review; avoids unsupported compounding but can miss synergy and underestimate disease |
| Multiplicative | product(RR_j) | Independent multiplicative effects on a specified scale; confounded univariate RRs do not identify it; can inflate risk severely |
| Additive excess | 1 + sum(RR_j−1) | Additive excess on a declared risk/rate scale; interactions and reference population remain unresolved |
| Capped multiplicative | min(cap, product(RR_j)) | Limits extremes but cap needs evidence and creates arbitrary saturation; not a probability repair |
| User combined risk | Directly supplied joint-pattern risk with horizon and reference | Useful with joint cohort evidence; otherwise adds unidentifiable parameters |
| Fitted multivariable model | Validated conditional prediction model | Can address covariance/interactions; requires external validation, transport and calibration and must not also apply the same individual RRs |
| Website raw sum | sum(RR_j), with baseline one if no factors | Audited legacy behavior only; counts the baseline more than once when several neutral factors are present |

Recommended sensitivities are additive excess, a justified cap, and observed combined-risk scenarios. Unrestricted multiplication is not a default. The old OR product of 2,916 must not be reactivated. Baseline progression, HIV strata and past-TB/radiographic history need reconciliation before any combination policy is numerically applied.

## Exclusions

Exclude recent infection RR 15 and recent close contact RR 15 from medical-factor prevalence and combination. Recency comes from infection history; contact may inform acquisition or exposure timing. Applying either again to a recency-dependent progression curve risks double-counting.

Do not import TSTin3D's fixed annual baseline risk of 0.1%, linear lifetime calculation to age 80, historical Styblo conversion, false-positive TST rules, or BCG/NTM interpretation. The site also offers AIDS, glucocorticoids and young infection age; these are outside the eleven-factor candidate registry. AIDS overlaps HIV and must not be added as an independent duplicate factor.

## Children and adults

Use separate reviewed baseline progression kernels for ages 0–4, 5–14 and 15+, preserving finer ages internally. Review [Martinez et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32199484/) as child-specific evidence: its contact-cohort population and age/time definitions must be reconciled before applying it to a general community. It does not validate adult factor RRs in children.

Pediatric factor prevalence may be stored while its progression effect remains unresolved. A reviewed numerical analysis must fail readiness when it needs an unreviewed child effect; it must not silently assign the adult value or treat the factor as absent.
