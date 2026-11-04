# Menzies 2008 publication reference

This is historical model reproduction and evidence review, **not clinical
validation or a contemporary recommendation**. It is not an implementation of
the current TSTin3D website. No website source code is copied.

Both formulations concern infection present at baseline. Neither generates future
acquisition, reinfection, transmission, migration, changing exposure, recurrence
after prior TB, or future changes in factors. Those mechanisms are not approximated
implicitly. Current active TB, prior disease and preventive treatment need a
separate model and remain unavailable in the ordinary first-event risk set.

## Published equations

For the historical annual infection probability ARI and contact term c:

```
pLTBI = 1 - (1-c)*(1-ARI)^age
PPV_TST = pLTBI / (pLTBI + FP_BCG + FP_NTM)
ARD = PPV_TST * 0.001/year * RR
CRD_to_age80 = ARD * (80 - age)
```

Units are explicit proportions. If the caller supplies historical smear-positive
incidence per 100000/year, the publication's coefficient49 yields percent ARI;
divide by100 again for a probability. [WHO TB impact measurement (2009), Annex 3](https://iris.who.int/bitstream/handle/10665/44231/9789241598828_eng.pdf),
clarifies the percent unit and the historic nature of this relationship. This
conversion exists only in the explicitly selected historical formulation. There
are no country downloads or lookup defaults.

The package records contact c=.40 for close and .08 for casual contact; explicit
none gives c=0. The caller must supply a source-bearing ARI or incidence and an
explicit full-age constant-exposure policy. Immigration mapping is unresolved;
provided immigration metadata blocks full PPV rather than truncating history.
Known country of birth/residence does not create an exposure history.

For non-BCG 10-14mm reactions, `FP_NTM=.027*supplied_NTM_sensitization`. For >=15mm,
the paper's NTM crossreaction rule gives zero NTM false-positive component.
Neither missing NTM prevalence nor unknown BCG becomes zero.

Table 1 supplies aggregate BCG false-positive values .085/.418, while the methods
describe separate size categories. That exact size-specific binding is not
recoverable from the publication. These values are retained in the package but
**not substituted automatically** for an individual 10-14 or >=15mm result.
BCG fixtures therefore have unresolved PPV. This deliberately avoids importing
the different numbers on today's About page. The paper implements positive TST
>=10mm; 5-9mm and IGRA are discussed as future work. IGRA-only profiles remain
unavailable. A supplied IGRA alongside TST is retained with an explicit TST-only
warning, never interpreted using current-site rules.

## Progression and limitations

The baseline is .001/year. Eleven medical/radiographic RRs are documented in the
[parameter table](parameter_evidence_m2a.csv), each with source, unit, definition,
population, evidence scope and effect-measure qualification. Smoking requires
>1 pack/day (the intersection of Table 2 >=1 and Figure 1 >1; exactly1 is
unresolved); underweight requires <90% ideal body weight; dialysis is not all
renal disease; infliximab is not every TNF inhibitor; transplant RR 74 requires
immunosuppression. These are historical associations, not verified causal HRs.

A single factor can be reproduced. The paper does not state a general rule for
concurrent factors; Table 3's granuloma/underweight pair is consistent with an
RR near5, but that is not evidence for a general sum rule. Compatibility mode
blocks multiple factors, including medical plus historical recency, instead of
choosing a modern policy. The modern default remains highest applicable RR;
additive-excess and multiplicative are explicit sensitivities only.

Table 2 assigns RR15 to recent infection/close contact for the first two years,
then the baseline rate. Confirmed close contact is treated as recent only under
this explicit historical policy. A supplied most-recent-infection duration is
needed to locate the remaining early period. Exactly two years is blocked because
the table and footnote do not define that annual-rate boundary consistently. Annual and wholly early-period
arithmetic can be reported. The long-horizon equation, footnote and contact
Table 3 outputs do not reconcile: crossing-period results and recent lifetime
risk stay unavailable. The separate raw footnote schedule is diagnostic, not an
adopted reconstruction. For example, hypothetical age35/tau0 gives .073 lifetime
risk at PPV1, whereas the paper prints .090 for Canadian contacts. No 5% rule from
the current website is imported to force agreement.

For ordinary non-recent profiles, requested risks use literal linear multiplication.
Raw results beyond probability bounds are retained in `historical_arithmetic`;
probability fields become unavailable, without clipping. Age80 gives zero remaining
historical lifetime. Above80 yields a negative raw expression and an unavailable
risk, never a silent zero. Requested horizons extending beyond80 are unavailable.
The assumption of survival to80 is historical bookkeeping, not a life table.
Treatment-benefit calculations are unavailable because complete equations are
not documented; no preventive-treatment recommendation is added.

## API and packages

```python
from tb_individual_risk import IndividualTBRiskModel, load_parameter_package
reference = IndividualTBRiskModel(
    load_parameter_package("menzies_2008_published_reference_v1"),
    formulation="menzies_2008_compatibility_v1",
)
hazard = IndividualTBRiskModel(
    load_parameter_package("synthetic_demonstration_v1"),
    formulation="hazard_based_progression_v1",
)
```

The default remains hazard-based with external-required clinical parameters,
exactly as in M1; selecting a package never silently selects a formulation.
Packages are frozen canonical records with file checksums pinned by the loader.
Returned records are independent copies. Editing an asset or mixing the historical
package into hazard mode fails. A change requires a new package version.

The hazard formulation remains `delta_H=H0(tau+horizon)-H0(tau)` and
`risk=P(TBI)*(1-exp(-m*delta_H))`, with explicit optional competing death. It is
not tuned to match historical linear outputs. Actual synthetic-package projections
and controlled same-number exponential calculations are distinguished in the reports.

Every result records formulation, package, source version, horizon, combination
policy, complete parameter record, scientific hash, warnings and unresolved
assumptions. All six input states remain distinct. There are no hidden historical
blank-field defaults or imputation priors. The old screenshot's combined
never-BCG/unknown option is documented but does not collapse unknown in this API.
