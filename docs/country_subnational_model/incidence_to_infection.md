# Incidence to infection identifiability

## Distinct quantities

Disease incidence I(t) is new TB disease per population-time. Force of infection lambda(a,s,t) is the instantaneous acquisition hazard per susceptible person-year; annual risk of infection is 1 − exp(−integral(lambda)). Dividing an incidence rate by 100,000 changes its units, not its biological meaning. Neither that division nor an unreviewed Styblo conversion identifies infection pressure.

The reconstruction combines recent disease following infection, reactivation after remote infection, reinfection, age at exposure, age-specific progression, changing medical factors, HIV, detection/treatment, migration and survival. An apparent incidence trend can arise from several of these without a corresponding proportional change in infection pressure.

## Proposed infection trend model

Fit a descriptive positive trend g(t) to an explicitly chosen local or WHO incidence series; normalize it at a reference year or over a stated period. Candidate:

    log lambda(a,s,t) = log k + beta log g(t) + log w(a,s,t) + delta(t)

Here k is an absolute annual infection-hazard scale, beta is the incidence-trend elasticity, w is an age/sex exposure function with population-weighted normalization, and delta is a constrained discrepancy term. A parsimonious first fit fixes beta and w under reviewed scenarios and sets delta to zero; beta = 1 is an explicit hypothesis, not an identity. With zero incidence, do not invent a small epsilon: use a nonnegative smoother with a documented observation error model or reject unsupported log fitting.

Retain WHO source estimates and bounds unaltered. Smoothing, COVID-period handling, missing years and territorial changes are analytical choices with their own provenance. Avoid over-weighting 25 highly correlated modelled annual estimates as if they were independent counts. A covariance-aware likelihood is preferable; otherwise use explicitly weighted calibration targets and label the objective as a discrepancy score rather than a true sampling likelihood.

## Cohort model and observation operator

Track age, calendar year, sex and time since infection. The susceptible compartment receives first infection at lambda. Infected people age through infection-duration states, may progress with hazard h(a_at_infection, a_current, tau, sex, medical history), may die, and may experience reinfection at susceptibility sigma times lambda. Reinfection updates the most-recent exposure clock under a named policy; partial immunity sigma, clearance and treatment states are reviewed parameters.

Predicted disease incidence is a forward sum over surviving infected cohorts and their progression hazards:

    predicted_cases(t) = sum_over_strata integral n_infected(a,s,t,tau,z)
                        * h(a−tau, a, tau, s, z, t) d(tau)
                        + explicitly modelled disease pathways outside that kernel
    predicted_incidence(t) = 100000 * predicted_cases(t) / population(t)

Do not count recent progression both in the infection flow and again in the duration kernel. WHO estimated incidence is compared with predicted disease onset using matched case definitions. A notification observation model, if used, additionally requires detection probabilities and delays and must not substitute notifications for WHO incidence.

For a transparent infection-only diagnostic, with no clearance, selective mortality or migration and a Poisson exposure process unaffected by prior infection, define H(a) as lifetime cumulative hazard and H_recent(a) as hazard within the chosen recent window r (bounded by birth). Then:

    P(uninfected) = exp(−H)
    P(recent) = 1 − exp(−H_recent)
    P(remote_only) = exp(−H_recent) − exp(−H)

These classify by most recent event; remote_only excludes a recent reinfection. For a first-infection-only model the corresponding recent fraction is exp(−H_before_recent) * (1 − exp(−H_recent)), which is different. Neither diagnostic is automatically the clinical TBI prevalence among survivors once progression, clearance or differential mortality is included. The forward cohort model must implement those exits consistently. Record first and last event semantics and do not silently exchange them.

Recent means time since the most recent modelled acquisition/reinfection is less than r; remote means at least r. The cutoff is a reporting/configuration choice with sensitivity at 2 and 5 years, not a fixed independent RR. Continuous time since infection drives the hazard kernel, avoiding an artificial homogeneous recent group.

## Calibration objective and anchor

Use the disease series to constrain the trajectory and an explicit natural-history kernel to connect scale to disease. With observed TBI:

    objective(theta) = incidence_discrepancy(theta)
                       + TBI_observation_discrepancy(theta)
                       + prior_or_regularization_penalties(theta)

TBI likelihoods may use binomial counts only when sampling assumptions justify them; weighted surveys need effective sample sizes/design effects or supplied covariance. Map modelled infection to the assay-defined observation through an approved measurement model. Aggregate predicted single-age prevalence using the survey's age/sex weights, not the model's convenience bands.

First estimate k with other assumptions fixed. Explore profile likelihood or posterior distributions, boundary hits, parameter correlations and prediction intervals. A free k plus free overall progression multiplier can trade off and must not both be claimed identified by incidence alone. Additional age exposure gradients, beta, backcast and clearance parameters require sufficient independent anchors or explicit priors.

Without an observed TBI anchor, fit only a parsimonious conditional model across declared natural-history scenarios and report the resulting prevalence as model-estimated and uncertain. Do not call the result an empirically identified TBI estimate.

## What is and is not identified

| Evidence | Information provided | Remaining ambiguity |
| --- | --- | --- |
| Annual disease incidence alone | Disease burden trajectory | Infection scale versus progression, historical stock, detection, HIV and other composition changes |
| Incidence plus fixed natural history | Conditional infection-pressure reconstruction | Dependence on fixed kernel, pre-WHO history, mixing, migration and survival |
| TBI by age at one survey | Cumulative exposure patterns | Cohort versus age effects; timing within life; viable infection versus persistent immune response |
| Serial TBI surveys or conversion cohorts | More temporal information | Measurement error, exposure selection and reinfection still matter |
| Age/sex disease plus HIV and treatment history | Additional natural-history constraints | Does not fully identify historical infection pressure or recent/remote stocks |

Disease incidence cannot uniquely separate recent and remote TBI. Even an age-specific TBI anchor constrains cumulative exposure more strongly than precise recency. Report uncertainty in the recent/remote decomposition separately from total prevalence. Local adult infection may have occurred elsewhere; migration histories can dominate a national backcast.

## Required pre-WHO scenarios

Let t0 be the first usable WHO year for the selected geography (often 2000, not universally).

| Option | Explicit rule | Required sensitivities |
| --- | --- | --- |
| Constant earliest level | Extend inferred infection trend at g(t0) back to cohort births | Compare earliest-year averaging window and all other defensible histories |
| Earliest fitted trend | Fit slope on an explicit earliest window and extrapolate before t0 | Window length, slope uncertainty and declared plausibility constraints; no silent caps |
| User historical series | Ingest dated incidence or infection-pressure observations with quantity and units | Coverage gaps, joins, ascertainment and incidence-to-infection mapping |
| Reviewed country scenario | Named country/history reconstruction approved by a reviewer | Alternative trajectory, migration and historical demographic assumptions |

A supplied historical incidence series still requires the calibration link; it is not automatically an infection-hazard series. Require complete history through each cohort's birth, territorial consistency and declared joining rules. Histories that imply implausible exposure or fail anchors return a failed scenario, not a clipped fit.

## Evidence and provenance

The proposed equations are a new architecture, not copied code or a validated fitted model. Historical exposure reconstruction combined with demography has precedent in [Houben and Dodd 2016](https://pubmed.ncbi.nlm.nih.gov/27780211/); its historical ARI estimates are not adopted as present-day defaults.

Record source snapshot IDs, hashes, population denominators, observation definition, backcast, earliest trend window, recent window, progression evidence, HIV/medical history scenarios, migration, mortality, clearance, reinfection susceptibility/reset policy, fitting objective, priors, optimizer version, seed, parameter bounds and compatibility diagnostics in the calibration artifact. Incidence and TBI uncertainty must propagate into the common population draw manifest.
