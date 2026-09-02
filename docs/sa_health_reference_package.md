# SA Health APY Working-Reference Package

This package freezes a non-dynamic APY LTBI screening health-economic
working-reference run for planning. The primary epidemiological anchor is a
Python stochastic event-ledger run using MATLAB-v9-compatible implicit
early/late natural-history semantics, validated distributionally against the
frozen MATLAB APY v9 reference scenario underlying the earlier plain SA Health
report.

Generate it from the repository root:

```powershell
python scripts/build_sa_health_reference_package.py
```

The generator uses the existing authoritative path:

```text
APY configuration
  -> MATLAB-v9-compatible Python stochastic APY model
  -> ltbi_screening_event_ledger_v3
  -> ltbi_health_economics_results_v3
  -> report-ready tables, figures and workbook
```

It does not use the dynamic community transmission model and does not modify
the MATLAB reference file.

## Reference Comparison

- Comparator: current practice, including ordinary background TB care, without
  additional systematic LTBI screening.
- Intervention: IGRA + 3HP targeted systematic LTBI screening.
- Targeting: prioritise people most likely to avoid active TB.
- Population: APY demonstration population of 1,500 residents.
- Eligibility: all 1,500 people are eligible; targeting determines who is
  screened.
- Coverage: 30% over a two-year screening window.
- Follow-up: 20 years.
- Primary model: stochastic APY v9-compatible simulation, 2,000 simulations,
  seed 1. Expected costs and effects use means across paired replicates.

## Economic Approach

- Perspective: Australian health-care system.
- Currency and price year: AUD 2019.
- Primary discounting: 3% for costs and health outcomes.
- Comparison discounting: 0%.
- Primary outputs: cost-consequence results, gross intervention delivery
  expenditure ratios and net health-system ratios after active-TB care offsets.
- Secondary outputs: provisional DALYs and ICER classification.
- No willingness-to-pay threshold is supplied; NMB and probability
  cost-effective are unavailable.

Programme setup, programme running, travel, outreach and staff support are
not yet locally costed. The primary working run uses the existing zero
compatibility values only so the calculation can proceed; this is not evidence
that implementation is costless. An illustrative AUD 500,000 setup scenario
is exported separately.

## Provisional Epidemiological Anchor

The inherited `10/770` active-TB target is treated only as a MATLAB software
calibration target for active TB observed over the historical two-year
screening window. Its scientific provenance and interpretation remain
unresolved. The implicit early phase is not measured APY recent-LTBI
composition. The older deterministic explicit recent/remote scenario is
preserved as a technical scenario and is not the primary SA Health
epidemiological estimate.

## Outputs

The default output directory is:

```text
outputs/sa_health_reference/sa_health_apy_matlab_v9_compatible_working_reference_igra_3hp_prevent_30pct/
```

The generated directory contains:

- `scenario_config.json`
- `economics_config.json`
- `manifest.json`
- `evidence_registry_snapshot.csv`
- `event_ledger_totals.csv`
- `event_ledger_annual.csv`
- `economic_replicates.csv`
- `economic_annual_by_arm.csv`
- `economic_summary.csv`
- `technical_recent_remote_scenario.json`
- `sa_health_apy_working_reference_outputs.xlsx`
- report-ready CSV tables in `tables/`
- portable SVG figures in `figures/`

The manifest records the code commit, configuration hashes, evidence-registry
hash, model and contract versions, unresolved inputs and interpretation
guardrails.
