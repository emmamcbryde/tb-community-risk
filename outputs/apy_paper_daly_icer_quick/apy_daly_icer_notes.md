# APY Paper DALY/QALY/ICER Outputs

Purpose: offline health-economic outputs for the APY paper/report.

Model basis: Python APY closed-cohort, person-level prevention model using direct prevented TB cases. Downstream transmission effects are not included.

Costing: Australian health-care system perspective, 2019 AUD, KWAB150/Dale-informed unit costs, 3% discounting assumption recorded for interpretation.

Program setup, program running, and false-positive incremental costs are set to zero in this offline paper script because the KWAB150 preset does not specify them. This should be sensitivity-tested if programme implementation costs are added.

DALY assumptions: active TB disability weight 0.333, active TB duration 0.5 years, configurable TB case fatality risk and YLL per TB death. YLL assumptions require explicit review.

QALY assumptions: healthy utility 0.8733, active TB utility 0.8182, SAE utility 0.8685, and LTBI treatment decrement sensitivity based on Dale-informed assumptions.

ICER formulas: incremental cost divided by DALYs averted, QALYs gained, or active TB cases prevented. Net monetary benefit equals health gain times WTP threshold minus incremental cost.

Limitations: not an online tool output; not a full transmission model; not a final policy gate; mortality/YLL assumptions are sensitivity inputs.

Lowest cost per DALY averted among calculable base-case strategies: prevent.

Lowest cost per QALY gained among calculable base-case strategies: prevent.

The APY model supports sequencing and prioritisation, not exclusion from care.

Run settings: quick=True, N=100, nReps=5, seed=1, screenCoverage=0.3.

Sensitivity rows generated: 28.
