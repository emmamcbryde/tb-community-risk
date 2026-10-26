# Explicit recent/remote TBI v2 — Milestone 3G status

Branch: **feature/explicit-recent-remote-tbi-v2**

Starting HEAD: **136de0556d0ce6ce0ce9cb76b04e22efa5e2cffc**

Starting working tree: clean; tracking origin/feature/explicit-recent-remote-tbi-v2,
0 ahead and 0 behind. Directory:
C:\Users\emmas\Documents\GITHUB\tb-community-risk-tbi-calibration.
Python 3.12.4; installed Streamlit 1.32.0. The repository pins Streamlit 1.39.0,
also used for supported-version rendered verification in an isolated validation
install without changing the installed environment.

Frozen SA Health release branch **release/sa-health-apy-he-v1.0.0** and peeled tag
**sa-health-apy-he-v1.0.0** remain at
**03cc16e4a52a55e10019dab16bd4f599571c4b90** locally and remotely.
The annotated tag object is **e99d36725dfc3ff40edc0ed68127e0b6cd927da5**.

This document replaces the previous milestone status with the current summary;
the earlier milestone summaries remain in Git history.

## Purpose and scientific boundary

Milestone 3G introduces an explicit population-level progression-risk scenario
among people who already have TBI at baseline. The motivating risk-enriched
population might have roughly 20 future comparator cases, whereas the unchanged
standard scenario projects about 2.22. This motivation does not supply an
observed target or identify a fitted multiplier.

The isolated app now defaults to a provisional 8× hazard scenario. The 8× value
is a working scenario assumption, not a calibrated, observed, estimated,
independently estimated or validated relative risk. The model does not solve
for 20 cases. Actual output is reported below.

Risk-factor-specific progression policy remains **none**. There is no inherited
prevalence/disease OR multiplication, product of multiple risk effects, use of
the old maximum 2,916 multiplier, retired continuous Markov pathway, or
Rising/Steady/Falling control in this explicit workflow.

## Versioned contract and computation

Contract identifier:

    population_tbi_progression_risk_adjustment_v1

Canonical configuration location:

    explicitRecentRemoteTBI.populationTBIProgressionRiskAdjustment

Required fields are contractVersion, version, mode, multiplier, source,
reviewStatus and notes.

| Mode | Resolved multiplier |
| --- | ---: |
| standard | Exactly 1.0 |
| risk_enriched | Exactly 8.0 |
| custom | Explicit positive finite value |

Zero, negative, missing, non-finite, boolean and non-numeric custom values are
rejected without clipping. A fixed mode with a conflicting supplied multiplier
is also rejected. The app explicitly selects risk_enriched, source
**Working scenario assumption**, reviewStatus
**provisional_requires_scientific_review**. General explicit runner builders
and pre-3G runner configurations retain the standard 1× reference behavior.

With unchanged central cumulative hazard H0 and baseline progression clock u:

    P(TB within t | baseline clock u) = 1 - exp(-m * [H0(u+t) - H0(u)])

The adjustment multiplies hazard in every future time-since-infection segment,
including the post-final-anchor extrapolation. It never computes m times the
final probability or clips such a product.

Stochastic inversion uses the same uniform draw and an unscaled future
cumulative-hazard threshold of -log(1-U)/m. Comparator and intervention continue
to share untreated progression potential. Effective TPT prevents only future
progression after its effective time.

Competing-mortality calculations multiply the TB hazard and TB cumulative
hazard by m while leaving the death hazard/survival unchanged. Constant,
externally supplied survival, and callable age/time mortality routes use that
same cause-specific semantics. No mortality data or new mortality control is
enabled in the isolated app.

Reset/no-reset clock selection precedes the adjustment. The deterministic
intervention helper now correctly retains the remote clock for recently
reinfected people with prior remote history under no-reset; they remain in the
recent prevalence group. This matches the existing stochastic clock contract.

The adjustment does not change TBI targets, calibration equations, fitted
exposure hazards, age weights, state allocation or infection-time distributions.
Uninfected people have zero prospective progression in this no-new-infection
workflow. Baseline active TB is separate and excluded from incident projection.
Test characteristics, coverage, initiation/completion probabilities, efficacy
and mortality inputs are unchanged. Progression before screening can change
how many people remain latent and eligible for treatment.

## Actual working-default projections

Configuration: 1,500 people; recent TBI 2%; remote-only TBI 8%; baseline active
TB 0%; 20 years; IGRA/3HP and the unchanged intervention/cascade assumptions;
reset clock; no competing mortality; no future infections or endogenous
transmission. Deterministic quadrature remains 16 infection-time by 12
screening-time points.

| Progression scenario | Recent comparator cases | Remote-only comparator cases | Total comparator cases | Intervention cases | Direct cases averted |
| --- | ---: | ---: | ---: | ---: | ---: |
| Standard 1× | 0.9650788868031669 | 1.250218587353016 | 2.2152974741561837 | 1.937301993328546 | 0.27799548082763775 |
| Risk-enriched 8× | 6.737556424892722 | 9.641908129825808 | 16.379464554718503 | 14.336705091021601 | 2.042759463696902 |

At ordinary display precision this is **2.22 versus 16.38 comparator cases**.
The 8× result is below eight times the standard count because of nonlinear
hazard-to-probability conversion. It is not forced to 20.

The enriched decomposition is:

| Baseline group | People | Average 20-year progression probability | Expected comparator cases |
| --- | ---: | ---: | ---: |
| Recently infected within 5 years | 30 | 22.5% | 6.74 |
| Remote-only infection | 120 | 8.0% | 9.64 |
| Uninfected | 1,350 | 0.0% | 0 |
| Baseline active TB | 0 | Not applicable | Excluded from incident projection |
| Total | 1,500 | — | 16.38 |

Counts and probabilities retain full precision internally. Hard accounting
checks verify baseline population totals, entered recent/remote TBI targets,
zero uninfected progression, baseline-active exclusion, group probability/count
identities, recent plus remote cases equalling comparator incidence, and
intervention plus averted cases equalling comparator cases. Checks run before
returning and before displaying results. Display-only component rounding
preserves the displayed case total and never feeds back into the model.

## Interface, sensitivity and state

Run from the repository root:

    python -m streamlit run explicit_tbi_ui\app.py --server.port 8502

The removed root explicit_tbi_app.py has not been restored. The isolated
directory has no pages subdirectory. Its four workflow entries remain:

1. Set up
2. Run analysis
3. Results
4. Evidence and technical information

Set up contains **Progression risk among people with TBI** beside known TBI
prevalence, with Standard progression risk — 1×, Risk-enriched population — 8×,
and Custom multiplier. Non-standard modes show the requested assumption
warning. One existing Restore working defaults action restores 8×.

Navigation and risk controls use stable widget identities and canonical
application state. Mode/custom selection survives reruns/navigation and marks
a completed result stale. The version-2 session envelope rejects older,
unmarked, frozen and legacy contamination rather than relabelling results.
State still uses the existing explicit-only namespace, preserving the legacy
economics guard.

Results show the expected-future-TB headline, selected mode/multiplier, group
decomposition and an analytical standard-versus-selected comparison. Every
other completed population, age, prevalence and intervention assumption is
identical in that comparison. It is explicitly labelled analytical even when
the main analysis is stochastic.

One click dispatches one workflow analysis. A selected multiplier other than
1× computes two intentional analytical projections: selected risk and its
standard-risk sensitivity. A resolved 1× run reuses its selected projection. Navigation and
ordinary reruns do not dispatch analysis. This distinction is tested rather
than misreporting the sensitivity as a second user run.

Mode, resolved multiplier and source/review/notes appear in canonical
configuration, runner hashes, the actual shared-runner calibration-cache key,
health signatures, staleness checks, stochastic replicate/aggregate metadata,
technical JSON and provenance text. Internal identifiers are confined to
technical evidence/downloads. Ordinary percentages/counts retain the 3F
formatting contract; raw JSON/CSV retain full precision.

Active-TB observations remain validation-only. Changing their numerator,
denominator or observation period does not change the multiplier, curve,
configuration health signature, or projection. Any observed-versus-expected
comparison is validation, not calibration.

## Numerical and rendered validation

The new suite covers the contract modes and invalid custom values; exact
hazard transformation and inverse sampling; preserved prevalence, age and
infection histories; both clock policies; mortality cause-specific math;
baseline/uninfected exclusions; zero coverage/efficacy; accounting failures;
staleness, reset and contamination; cached configuration identity; precision;
four-only navigation; and economics/ledger rejection.

For Monte Carlo comparison, the acceptance rule was declared before running:
20 seeds, N=1,500 per seed, and absolute mean-count difference at most
max(0.5 cases, 4 × MCSE). Deterministic references use the workflow's 16×12
quadrature.

| Mode | Metric | Deterministic | Stochastic mean | MCSE | Acceptance tolerance |
| --- | --- | ---: | ---: | ---: | ---: |
| 1× | Comparator cases | 2.2153 | 1.70 | 0.2417 | 0.9668 |
| 1× | Intervention cases | 1.9373 | 1.45 | 0.2458 | 0.9830 |
| 1× | Direct cases averted | 0.2780 | 0.25 | 0.0993 | 0.5000 |
| 8× | Comparator cases | 16.3795 | 17.05 | 1.0525 | 4.2100 |
| 8× | Intervention cases | 14.3367 | 14.55 | 0.7996 | 3.1984 |
| 8× | Direct cases averted | 2.0428 | 2.50 | 0.4947 | 1.9788 |

All declared reconciliation metrics passed. Additional no-reset/constant
mortality tests check the three epidemiological case outcomes using the same
MCSE rule, and verify identical death draws at 1×/8×. Existing ancillary
cascade approximations under externally supplied mortality are not relabelled
as validated by those three case checks; the isolated UI still uses no
mortality data.

The central curve/anchors are compared directly with the starting commit, and
the standard reset projection is compared with the starting runner. Both are
unchanged. Shared runner/simulation, frozen reference integrity, and bounded
legacy event-ledger/economics regressions are also run.

No interactive browser surface was available. AppTest provides rendered UI
checks and the live localhost protocol audit provides real rendered elements,
tables, page-registry and dispatch-count checks, not HTTP 200 alone.

Generated evidence remains under **abm/output/milestone3g/**, including the
full working-default technical JSON and MCSE reconciliation records.

### Completed verification

| Test group | Result | Runtime |
| --- | --- | ---: |
| New risk mathematics plus installed-version UI/isolation regression | 47 passed, 1 version-dependent frozen-router skip | 264.73 s |
| Existing explicit mathematics, shared runner/simulation, frozen reference | 171 passed | 149.46 s |
| Existing bounded event-ledger/economics regressions | 71 passed | 597.40 s |
| Complete supported Streamlit 1.39 isolated/risk/reference suite | 66 passed | 469.12 s |
| Added saturated-custom, accounting-failure and rendered regression | 8 passed | 69.30 s |
| Actual Streamlit 1.54 progression-risk rendered checks, final UI | 6 passed | 55.33 s |

The extra saturated-custom test accepts an explicitly supplied 1,000,000×
multiplier without clipping. Its derived group average can exceed one by a
few floating-point ulps after weighted quadrature/count division. The hard
check allows only a 1e-12 upper-bound numerical tolerance and retains the raw
value. This changes neither the projection nor the multiplier. Small positive
custom multipliers display in scientific notation when needed instead of
appearing as 0×. Full-precision configuration is retained.

The originally running port-8502 app was in a different environment: Python
3.12.12, Streamlit 1.54.0, NumPy 1.26.4 and pandas 2.2.2. That existing isolated
development server was restarted in the same environment to remove cached 3F
imports. No package installation or repository-wide Streamlit change was made.
The live audit supports explicitly attaching to that server and the 1.54
string-valued radio protocol; it checks the rendered 8× selection, not merely
the radio's default-index field. Streamlit 1.54 identifies a single root page
with an empty name, whereas older versions call it app.

The frozen live server again rendered its original nine navigation entries
and populated Results tables using Streamlit 1.39. In 1.54 a direct request
for an unregistered legacy page can echo the requested alias/hash in navigation
metadata without emitting the older page_not_found message. The audit verifies
that mainScriptPath still resolves to explicit_tbi_ui/app.py, the four-option
workflow remains, and the actual body is the isolated pre-run summary with no
Economics page. The echoed alias does not execute a legacy page.
The separate direct legacy-economics AppTest rejects explicit state before any
economic display or ledger consumption.

The final complete live Streamlit 1.54 audit passed and is recorded in
abm/output/milestone3g/live-complete/isolated-1.54.0-audit.json. Initial render
took 1.17 s. The first and second unchanged deterministic clicks took 34.34 s
and 33.34 s; dispatch counts were exactly 1 and 2 respectively. Each dispatch
deliberately includes the selected 8× projection and a separate 1× analytical
sensitivity, not a duplicated user run. Selected quadrature took 16.54/15.13 s,
the analytical sensitivity 16.39/16.52 s, calibration 0.028/0.014 s and result
construction 1.33/1.60 s. Results rendering took 0.09/0.12 s. The rendered
recent and remote contributions were 6.74 and 9.64 cases, summing to 16.38.
The app remains available on port 8502. No interactive browser surface was
available; AppTest and live Streamlit protocol/Arrow rendering were verified.

Python compilation and Git whitespace checks pass. No timed-out or unrun test
is counted as passed. The installed-version frozen-router skip is covered by
the supported-version live frozen-router check.

Changed files: app/explicit_tbi_workflow.py, app/explicit_tbi_session.py,
app/explicit_tbi_formatting.py, engine/apy/population_tbi_progression_risk.py,
engine/apy/explicit_recent_remote_progression.py,
engine/apy/explicit_recent_remote_runner.py, explicit_tbi_ui/app.py,
explicit_tbi_ui/README.md, scripts/validate_explicit_tbi_ui.py,
tests/test_population_tbi_progression_risk.py,
tests/test_explicit_recent_remote_tbi.py,
tests/test_explicit_tbi_streamlit_app.py,
tests/test_explicit_tbi_ui_isolation.py, this status document and the explicit
TBI foundation specification. Historical scope tests now allow the existing
legacy economics rejection guard while continuing to prohibit legacy pages
from importing or exposing the isolated workflow.

## Future risk-factor modelling and guardrails

Future risk-factor-specific progression modelling needs reviewed RRs/HRs,
correlation handling and prevalence among TBI-positive people. It should
replace or explicitly re-normalise/contextualise the aggregate scenario and
avoid double counting population enrichment. The inherited OR product is not
a shortcut to that work.

No economics, ICER, DALY, QALY, cost change, costable ledger, new infection,
endogenous transmission, MATLAB change, dynamic-model change, deployment,
merge or tag is introduced. The deliberately unavailable/non-costable ledger
marker remains. Frozen SA Health files, release pointers and artifacts are
preserved. This model supports planning and sequencing, not denial of care.
