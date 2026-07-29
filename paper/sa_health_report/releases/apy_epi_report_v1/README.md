# APY Epidemiology Report Release v1

This directory freezes the APY targeted LTBI screening epidemiology analysis used to reproduce the CHO-facing report graph data derived from `APY_CHO_epidemiology_outputs.xlsx` and `APY_CHO_epidemiology_outputs_v2.xlsx`.

## Command

Run from the repository root:

```powershell
python paper/sa_health_report/releases/apy_epi_report_v1/run_locked_apy_epi_report.py --output-dir paper/sa_health_report/releases/apy_epi_report_v1/artifacts/run1
```

The committed verification used:

```powershell
python paper/sa_health_report/releases/apy_epi_report_v1/run_locked_apy_epi_report.py --output-dir paper/sa_health_report/releases/apy_epi_report_v1/artifacts/run3 --compare-existing
python paper/sa_health_report/releases/apy_epi_report_v1/run_locked_apy_epi_report.py --output-dir paper/sa_health_report/releases/apy_epi_report_v1/artifacts/run4 --compare-existing
python paper/sa_health_report/releases/apy_epi_report_v1/verify_locked_release.py --run1 paper/sa_health_report/releases/apy_epi_report_v1/artifacts/run3 --run2 paper/sa_health_report/releases/apy_epi_report_v1/artifacts/run4
```

The verification workflow runs the same command twice with separate output directories and compares numerical CSV/XLSX content. Timestamps are stored in JSON manifests, not numerical CSV outputs.

## Backend

The locked analysis uses the Python APY backend (`python_apy_v9_port`). MATLAB is not used for this release run.

## Frozen Model State

The locked analysis was run from git commit `cd7727e94b5ca577db42c5de71291617d9675c21` on branch `codex/lock-apy-epi-report-run-v1`.

The first clean locked run recorded its environment at `2026-07-29T11:10:38.976789+00:00` and its run manifest at `2026-07-29T11:59:21.172180+00:00`.

## Interpretation Limitations

Outputs support planning and sequencing, not exclusion from care. They represent direct person-level prevention benefits only and exclude downstream transmission benefits. The 1% LTBI and prevalence-update scenarios are sensitivity analyses, not accepted APY base-case estimates unless reviewed locally. Stopping and prevalence-update outputs are review triggers, not automatic stop rules.

Generated replicate-level artifacts are intentionally excluded from git under `artifacts/`; checksums are recorded in committed manifests after each locked run.
