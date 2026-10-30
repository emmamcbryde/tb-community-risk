# Milestone 3D empirical comparisons

**Internal diagnostic validation only. No country TBI estimate, central parameter
approval or external validation is claimed.** The existing user-TBI pathways
remain available; production incidence inference remains unavailable. Read the
[protocol](milestone3d_protocol.md) before interpreting these results.

## Immutable observations

Snapshot: `data/country_tb_validation/snapshots/tbi-validation-fe05d6d1ab8c/`.
It contains96 selected source observations across three studies, including
overlapping descriptive age/sex/assay totals. **96 is not a count of independent
surveys.** Splits prohibit using paired assays, parent totals or overlapping
nativity groups as independent validation. Household clustering remains explicit.

Each source record retains citation, DOI/URL, access date, population, survey
period/year, age bounds, sex, raw numerator/denominator when actually reported,
rounded or weighted percentage, source95%CI if present, assay generation,
threshold, sampling frame/design, migration fields, exact table cells, licence
and extraction notes. Counts, weighted estimates and missing values remain
distinct. Transformed inputs carry a SHA-256 of the original record and a
scientific hash; source JSON and derived JSON have separate byte checksums.
Original full texts/workbooks are not added to the repository.

| Study | Selected observations and ascertainment | Appropriate interpretation |
| --- | --- | --- |
|[Innes AL, Nguyen ST, Lebrun V, et al. Public Health Action.2023. doi:10.5588/pha.23.0020](https://pmc.ncbi.nlm.nih.gov/articles/PMC10446657/)|Viet Nam household-contact campaigns2020/2021; Table3 exact paired-valid QFT-Plus and TST5/10mm counts, age/sex/overall. QFT-Plus manufacturer .35IU/mL algorithm;5TU TST. CC BY4|Convenience contact subsets in different districts, not a national sample or longitudinal population. Baseline CXR/Xpert screening; TST eligibility excludes under5/pregnancy; indeterminates not positivity |
|[Miramontes R, Hill AN, Yelk Woodruff RS, et al. PLoS One.2015;10:e0140881](https://doi.org/10.1371/journal.pone.0140881)|US NHANES2011–2012 Table4 QFT-GIT weighted prevalence/CI by age,sex,nativity. Public-domain article|Probability survey provides an informative contrasting low-incidence setting. Published weighted percentages do not imply binomial counts. US-born analysed separately; foreign-born exposure/arrival histories missing |
|[Laemmle-Ruff I, Graham SM, Williams B, et al. Emerg Infect Dis.2022;28. doi:10.3201/eid2809.212426](https://doi.org/10.3201/eid2809.212426)|Children's Australian visa applications, Table2 valid application-level results; age/sex and separate test-level assay totals. Public-domain/attributed source|Premigration screened applicants, not Australian residents or unique people. Mixed assays; repeated applications possible. No adult natural-history or assay performance transported to children |

The Australia table's youngest label is literally `<4`, followed by5–9. Its
coverage ambiguity is preserved in `transformation_notes`; no inference about
age4 is adopted. The age bounds are a provisional indexing convention only and
the records are unavailable for model anchoring. Test-level counts are not
substituted for application-level denominators.

Other candidates remain documented rather than invented: Marks et al.'s2018
Ca Mau community sample requires exact age-table access/rights; Vo et al.'s
two-province cohort is mixed/selected and its noncommercial licence is not
assumed to permit repository redistribution; Australian refugee clinics need
age/assay extraction and eligibility review. The Dale Australian TBI model is
excluded as independent validation because it already infers TBI from incidence
and migration. Korea is an intermediate-setting data-required candidate.
NHANES reanalyses of the same cohort are not additional independent datasets.

## Assay model and age allocation

Predicted positivity = Se×latent TBI +(1-Sp)×(1-latent TBI).
The [USPSTF2023 evidence review](https://doi.org/10.1001/jama.2023.3954)
supports fixed adult candidate operators QFT-Plus Se.89/Sp.98 and QFT-GIT
Se.81/Sp.99. Their sensitivity largely uses active-disease reference samples;
there is no gold standard for viable latent infection. Country/BCG,
immunosuppression and pediatric transport remain unapproved.
Named one-way source-bound scenarios retain their source review, including
[Oh et al.2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8423471/).
These are not joint assay confidence intervals. Se/Sp are never fitted withalpha.
TST, QFT-GIT, QFT-Plus and mixed assays are not interchangeable.

Analysis is at the survey year, not labelled2024. Within finite supplied age
bands, uniform continuous age weights are an **explicit diagnostic assumption**;
there are no band midpoints. No sex counts are invented for both-sex records.
Bands crossing childhood and open-ended bands remain unresolved. Sex marginals
are retained, but actual joint age/sex residuals require joint data/weights;
they are unavailable in these age-only fits, not assumed zero.

Pressure uses the retained deterministic WHO log-incidence smoother and normalized
shape, not an estimated force of infection. The earliest observed normalized
shape is held constant before2000; named bounded-backward-trend sensitivity has
slope−.01/year, bounds.01–100, reconstructed from1900. It is not fitted to
validation residuals. No WHO shape is automatically extrapolated after2024.

## Direction A: TBI fit, held-out TBI and forward comparison

Prespecified splits fit one age band, then compare disjoint remaining ages.
No observed incidence enters the TBI scale fit. Contact metrics use tested
sample counts as descriptive weights; NHANES uses equal age-record weights
because appropriate population weights are not published in the selected table.
Neither metric is a national population-weighted performance estimate.

| Data/split | Fittedalpha/year | Held-out weighted RMSE | Unresolved held-out tested count |
| --- | ---: | ---: | --- |
| Viet Nam contacts2020: fit30–39 |.01121195425|8.2549 percentage points|302:0–19 banded and70+ |
| Viet Nam contacts2021: fit30–39 |.00487934296|2.7085 points|354:0–19 banded and70+ |
| US-born NHANES: fit15–24 |.000374332379|1.1842 points|Unknown: weighted counts not reconstructed;6–14 and65+ unresolved |
| Fit2020 contacts, compare separate2021 sample |.01121195425|22.1399 points|Children/mixed-child bands and70+ unresolved |

Counts excluded from held-out metrics are not excluded from snapshot population
reconciliation. Unreported ages/sexes in the source are separately noted; the
age/sex marginals need not sum to the overall sample. For example, the2020
paired sample has1330 records but1327 with reported age. Do not “repair” these
source missingness differences. No full-population prevalence is reported.

| Sample | Age band | Observed assay positivity | Predicted assay positivity | Residual(predicted−observed) |
| --- | --- | ---: | ---: | ---: |
| Viet Nam2020 |20–29|30.400%|29.302%|−1.098points |
| |40–49|42.661%|47.408%|+4.747 |
| |50–59|44.656%|54.284%|+9.627 |
| |60–69|48.768%|60.023%|+11.254 |
| Viet Nam2021 |20–29|17.708%|14.940%|−2.769 |
| |40–49|23.973%|25.718%|+1.745 |
| |50–59|27.027%|30.504%|+3.477 |
| |60–69|37.013%|34.928%|−2.085 |
| US-born NHANES |25–44|1.900%|3.240%|+1.340 |
| |45–64|3.600%|4.605%|+1.005 |

US predictions fall within both published held-out95%CIs: descriptive interval
inclusion2/2, not an independent coverage probability or acceptance claim.
The contact table does not report these intervals; none are invented. Age
decreases in source/prediction are exported, including descriptive source open
ages. Increasing age is not itself a validation criterion: migration, cohort
effects, selection and assay variation can produce non-monotone observations.
The growing2020 older-age residual and cross-campaign mismatch challenge
transport of a shared national historical shape to selected contact populations.
Children, old ages and sex-specific fit cannot be claimed validated.

| One-way scenario |2020 contactalpha / RMSE points|2021 contactalpha / RMSE points|US-bornalpha / RMSE points|
| --- | --- | --- | --- |
| Lower specificity |.01062401 /7.6469|.00428519 /2.5066|.0000627716 /.7774 |
| Higher specificity |.01141952 /8.4675|.00508911 /2.8803|.0005286969 /1.8724 |
| Bounded backcast |.01076152 /11.0793|Not prespecified|.0003654190 /1.5831 |

Low US positivity is especially sensitive to specificity: a better residual
under one operator does not select or validate that operator. The source assay
parameters remain fixed. Backcast changes older-age profiles and must not be
chosen solely to improve fit.

Forward contextual rates use the unchanged TBI fit, adult RR1 reference curve,
finite adult bands and explicit sample/equal-record weights. No risk-factor
prevalences are invented. The independent last-event quadrature integrates
instantaneous baseline hazard, so mortality does not enter this rate.

| Sample |Conservative|Central candidate|Higher progression|Same-year WHO contextual rate|
| --- | ---: | ---: | ---: | ---: |
| Viet Nam contacts2020 |88.629|100.761|114.000|171 |
| Viet Nam contacts2021 |40.144|45.713|51.625|171 |
| US-born NHANES2012 |4.344|4.920|5.518|3.8 |

Units are expected instantaneous first-episode cases/100000 disease-freePY
for the covered diagnostic adult reference sample. WHO rates are annual
country new/recurrent episodes and cover different populations. Differences
and ratios are exported but are **not compatible incidence validation** and
not evidence of country under-/overprediction. Mortality, children, migration,
prior disease, risk factors and joint demography remain missing for a matched
country projection. No cumulative or long-horizon country estimate is made.

## Direction B: incidence fit, independent observed TBI

The utility uses the existing3C gate, not a new exception. It requires explicit
provenance that validation TBI was not used to construct the incidence target;
country mismatch and participant leakage are refused. It produces achieved
versus observed positivity only if a compatible, fully reviewed diagnostic
configuration passes. A fictional closed-population test recoversalpha.017
within1e-10; independent fictional TBI remains unchanged despite disagreement.

**Actual Viet Nam and Australian WHO targets are blocked_scientific_gate.**
No empirical direction-B scale or TBI profile is adopted. Observation episodes,
acquisition, pediatric/remote progression, effects, mortality, migration,
historical depletion and reviewer thresholds are not resolved. Country gate
demonstrations are clearly labelled incomplete placeholders, not national models.
The old fictional Viet Nam-shapedalpha≈.0181993/TBI≈65.33% is never used as an
empirical anchor or presented as a national estimate.

## Interpretation and discrepancy

These analyses provide empirical **model checks**, not external sign-off.
The code reports absolute/relative error, ratio, age/sex residuals when defined,
explicit-weight RMSE, interval inclusion, incidence mismatch, unresolved tested
counts and compatibility. Missing observations are not zero. Source-data,
assay/parameter, scenario and numerical uncertainty remain separate.

The discrepancy contract explicitly records WHO estimation, observation
incompatibility, migration, importation, recurrence, prior TB, omitted factors,
assay error, incidence revisions, shape misspecification, endogenous transmission,
aggregation and mortality. Magnitudes are unknown without component-specific
data; their total is not fitted intoalpha. No causal explanation is inferred
from a residual alone.

Candidate review thresholds are absolute-error2/5/10points, RMSE2.5/5/10points,
interval inclusion80/90/95%, and unresolved-population0/.1/1%. None is selected
or approved. Study design and sampling error must inform any eventual choice.
The [reviewer table](milestone3d_review.md#reviewer-response) records pending
judgements separately from software tests.
