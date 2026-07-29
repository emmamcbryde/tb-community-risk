# APY Epidemiology Report Release v1

This directory freezes the APY targeted LTBI screening epidemiology analysis used to reproduce the CHO-facing report graph data derived from `APY_CHO_epidemiology_outputs.xlsx` and `APY_CHO_epidemiology_outputs_v2.xlsx`.

## Command

Run from the repository root:

```powershell
python paper/sa_health_report/releases/apy_epi_report_v1/run_locked_apy_epi_report.py --output-dir paper/sa_health_report/releases/apy_epi_report_v1/artifacts/run1
```

The verification workflow runs the same command twice with separate output directories and compares numerical CSV/XLSX content.

## Backend

The locked analysis uses the Python APY backend (`python_apy_v9_port`). MATLAB is not used for this release run.

## Frozen Model State

The release was prepared from git commit `70d842ef1b51175354943bd4250bea954165e560` on branch `codex/lock-apy-epi-report-run-v1`.

## Interpretation Limitations

Outputs support planning and sequencing, not exclusion from care. They represent direct person-level prevention benefits only and exclude downstream transmission benefits. The 1% LTBI and prevalence-update scenarios are sensitivity analyses, not accepted APY base-case estimates unless reviewed locally. Stopping and prevalence-update outputs are review triggers, not automatic stop rules.

Generated replicate-level artifacts are intentionally excluded from git under `artifacts/`; checksums are recorded in committed manifests after each locked run.
