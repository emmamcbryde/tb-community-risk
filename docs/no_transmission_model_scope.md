# Scope and intended use of the no-transmission screening model

Status: scope statement adopted for branch `feature/generic-no-transmission-model`.
The applicability framework in section 5 is **provisional and requires scientific
review**. This document changes no model calculations, epidemiological or economic
logic, or frozen artifacts.

## 1. Scope statement

This model evaluates TB screening, diagnostic, preventive-treatment and
health-economic strategies for a defined population. It rests on one structural
assumption:

> Onward local TB transmission generated during the analysis horizon does not
> materially affect the comparison between the strategies evaluated.

The model therefore estimates **direct outcomes among the modelled population**:
people screened, test results (true and false positives), preventive treatments
started and completed, adverse events, infections effectively treated, active TB
cases in the modelled population under each strategy, and the resulting costs,
DALYs and cost-effectiveness.

The model does **not** estimate:

* indirect (transmission-mediated) effects of any strategy;
* herd or population-level protective effects;
* secondary infections or secondary cases prevented;
* changes in force of infection, annual risk of infection or community incidence;
* the effect of interventions on outbreaks or transmission clusters.

A difference in active TB between strategies is a difference in disease
**within the modelled population**. It is not an estimate of the TB burden
removed from the wider community.

The model is intended primarily for **low-incidence communities, or other clearly
defined populations, where the no-transmission assumption is defensible** (for
example, a screening programme for recent migrants in a setting where most disease
arises from reactivation of infection acquired elsewhere). It is for planning and
sequencing, not for denying care.

## 2. What the no-transmission assumption means in this engine

As implemented in the individual-based engine (`engine/apy/simulation.py`):

* infection status is assigned at baseline from the prevalence of TB infection
  (with an optional recent-versus-remote history, see
  `recent_remote_infection_history.md`);
* no new infection and no reinfection are simulated during follow-up;
* active TB in one modelled person has no effect on the risk of anyone else,
  inside or outside the modelled population;
* estimated TB disease incidence attached to a population profile (WHO or local)
  is descriptive only (`incidenceUsedByEngine: false`); it does not drive
  infection pressure.

Two distinct assumptions follow, and both must hold for results to be
interpretable:

1. **No outward transmission effect.** Cases prevented in the modelled
   population would not have generated enough onward transmission, over the
   horizon, to change the ranking or cost-effectiveness of the strategies.
2. **No material ongoing exposure.** Members of the modelled population face
   little new infection or reinfection risk during the horizon. Where local
   exposure is substantial, the benefit of preventive treatment can be
   overstated (treated people may be reinfected) and baseline risk understated.

An optional *background exposure to TB infection* (an exogenous, catalytic force of
infection with no feedback) could relax assumption 2 without adding transmission dynamics. It
is specified, not implemented, in `catalytic_infection_pressure_spec.md`.

Assumption 1 usually makes results **conservative** for interventions that
reduce disease (indirect benefits are omitted). Assumption 2 can bias results in
**either direction**. Neither bias is quantified by this model.

## 3. When the model may be inappropriate

The assumption is about the population and setting, not about a country's
headline incidence. The model may be inappropriate even in a nominally
low-incidence country when, for example:

* there is a current outbreak, or genomic or epidemiological evidence of
  sustained local transmission clusters;
* the modelled population is a high-incidence sub-population (for example some
  First Nations, remote, homeless, incarcerated or institutional populations)
  inside a low-incidence country;
* the population mixes intensively with a high-transmission group (household,
  institutional or occupational contact);
* the intervention is large enough relative to the local epidemic that its
  indirect effects could plausibly matter;
* the analysis horizon is long enough for onward transmission chains to
  accumulate.

**The model should not be selected on incidence alone.** Incidence is one input
to the applicability assessment in section 5.

## 4. Relationship to the dynamic (Starsim) workstream

A separate dynamic transmission model, developed with Starsim in a separate
workstream, is the appropriate pathway for high-incidence settings and for any
population with meaningful ongoing transmission, where interventions can change
the force of infection and produce indirect population effects.

* Starsim is **not** a dependency of this repository branch, and no Starsim
  files are included or modified here.
* The Python dynamic model in `engine/dynamic/` (used by the legacy
  `streamlit_app.py` research pages) and its readiness plan
  (`dynamic_model_readiness_spec.md`) predate this scope decision. They are
  retained unchanged for backward compatibility. Transmission modelling for new
  analyses is expected to go through the Starsim workstream; whether the
  in-repository dynamic model is retired, kept for comparison only, or aligned
  with Starsim is an open decision (section 8).
* The APY results field `technical.dynamicComparison` is a re-expression of
  individual-based results in the dynamic model's output vocabulary. It is not
  produced by a transmission model and contains no transmission effects.

## 5. Provisional applicability framework (requires scientific review)

No incidence value is presented here as a validated cut-off. Any numerical
threshold is **provisional** until the validation work in section 7 is done.

### 5.1 Terminology: WHO categories are not an applicability boundary

* **WHO low-incidence terminology.** The WHO action framework for low-incidence
  countries uses fewer than 100 cases per million (10 per 100,000) per year for
  "low incidence", with lower thresholds for pre-elimination and elimination.
  Separately, WHO guidance on programmatic TB-infection testing and treatment
  has used an estimated incidence below 100 per 100,000 per year as one criterion
  for systematic testing in some risk groups. National programmes use other
  values (for example, about 40 per 100,000 per year for some migrant screening
  policies). *Exact values and citations are to be confirmed against current
  WHO and national sources during review.*
* **Model-specific applicability boundary.** The question for this model is
  different: *could transmission effects over this horizon change the
  comparison between these strategies in this population?* That depends on
  transmission, not on a programme-eligibility label. A WHO category, or a
  value such as 40 per 100,000, is at most a screening flag for further
  assessment, not a test of model validity.

### 5.2 Assessment domains

An applicability assessment should record, with sources and uncertainty:

| Domain | Evidence to consider | Signals against the no-transmission model |
| --- | --- | --- |
| Incidence | Estimated incidence in the modelled population and in the surrounding community; trend | High or rising incidence in either |
| Local transmission | Notifications in children (a proxy for recent transmission), contact-investigation yield, recent-infection indicators | Sustained paediatric cases; high contact yield |
| Clustering | Genotyping or whole-genome sequencing cluster rates and sizes | Large or growing clusters; recent cluster links |
| Imported versus locally acquired disease | Share of cases in people born in high-incidence countries; time since arrival; genomic linkage to local cases | Majority locally acquired |
| Population mixing | Household, institutional and social contact between the modelled population and higher-transmission groups | Intense mixing with high-transmission groups |
| Intervention scale | Share of the local infectious pool the intervention could reach | Intervention large relative to the local epidemic |
| Analysis horizon | Length of follow-up relative to generation time and cluster growth | Long horizons with non-negligible transmission |

### 5.3 Provisional decision framework

| Setting | Provisional guidance |
| --- | --- |
| **Very low incidence with negligible local transmission** (disease mainly from reactivation of infection acquired elsewhere; little clustering; no outbreaks) | The no-transmission model may be appropriate. Record the assessment and state the assumption in every output. |
| **Intermediate or uncertain** (moderate incidence, incomplete evidence on transmission, clustered sub-populations, or long horizons) | Requires an explicit, documented assessment. Preferably compare with the dynamic model; if the strategy ranking or cost-effectiveness conclusion changes, report the dynamic result as primary or present both. |
| **High incidence, or meaningful sustained local transmission** at any incidence | The dynamic model is preferred. Results from this model, if shown, must be labelled as direct effects only and not used as the primary estimate. |

The boundaries between these rows are deliberately qualitative. They are not
implemented as a gate in the software.

## 6. Outputs this model must never describe as estimated

Unless a transmission model is used, outputs, reports, figures and user-facing
text must not present the following as estimated by this model:

* secondary infections prevented (averted);
* secondary active TB cases prevented, or cases prevented in contacts;
* population-level or community-level transmission reduction;
* change in force of infection or annual risk of infection;
* change in community (population) TB incidence, or incidence trajectories
  under intervention;
* herd, indirect or spill-over effects;
* progress towards elimination targets attributable to the intervention;
* outbreak prevention or cluster reduction;
* total cases averted "in the community" (as opposed to in the modelled
  population).

Acceptable wording names the population: for example, "active TB cases averted
**in the screened cohort** (direct effect only)". Where it helps, outputs may
say that indirect benefits, if any, are not included.

## 7. Evidence and validation needed before a formal applicability gate

Before any threshold or rule is implemented as a software gate:

1. **Literature and guideline review.** Confirm WHO and national thresholds and
   their intended purpose; review evidence on the share of disease from recent
   local transmission in low-incidence settings (genomic epidemiology,
   paediatric notifications, migrant cohorts).
2. **Structural comparison with the dynamic model.** Run matched scenarios in
   this model and the Starsim model (same population, natural history, test
   accuracy, cascade, costs and horizon) across a grid of incidence, local
   transmission share, mixing, intervention scale and horizon.
3. **Decision-relevant error metric.** Define in advance what counts as a
   material difference (for example a change in strategy ranking, an ICER
   crossing a willingness-to-pay value, or a relative difference in direct plus
   indirect cases averted above an agreed tolerance).
4. **No-transmission limit check.** Confirm that the dynamic model with
   transmission switched off reproduces this model's direct outcomes within
   Monte Carlo error, so differences can be attributed to transmission.
5. **Boundary estimation with uncertainty.** From steps 2-4, estimate the region
   where the no-transmission approximation is adequate, with uncertainty, and
   state which domains (section 5.2) matter most.
6. **Review and sign-off.** Scientific review of the boundary, its wording and
   its evidence, before it is shown to users as more than advisory.
7. **Implementation.** Only then add a recorded applicability assessment to the
   population profile and an advisory (not blocking) check in the interface,
   with tests.

## 8. Open decisions

1. Adoption and wording of the provisional framework in section 5.
2. Whether any numerical flag (for example 10 or 40 per 100,000 per year) is
   shown to users, and how it is labelled.
3. The future of `engine/dynamic/` and `dynamic_model_readiness_spec.md` given
   the Starsim workstream.
4. The scientific decisions for optional background exposure
   (`catalytic_infection_pressure_spec.md`, section 4).
5. Standard direct-effects wording now lives in `engine/model_scope.py` and
   `app/general/terminology.py`; reports built outside the general application
   should adopt the same wording.
