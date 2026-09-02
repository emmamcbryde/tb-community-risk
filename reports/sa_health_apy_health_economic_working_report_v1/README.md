# SA Health APY health-economic working report

This directory contains the editable Word report, PDF verification copy,
report-ready CSV tables, figures and a concise reproducibility manifest.
Large event-ledger and economics CSVs are stored in
`technical_event_ledger_and_economics_csv.zip` to keep the committed report
package portable.

Formatting reference: `paper/sa_health_report/APY_SA_Health_LTBI_model_report_draft.docx`
and the locked epidemiology-only release under
`paper/sa_health_report/releases/apy_epi_report_v1/`.

Reproduce from the repository root:

```powershell
python scripts/build_sa_health_reference_package.py
python scripts/build_sa_health_word_report.py
```

The report is a working analysis for planning and refinement. It is not final
policy evidence and does not claim cost-effectiveness.

Code commit: `12f697aa7caccad9a9290a7c7eef2dad54257e46`
Configuration hash: `a84f98039b269502d17f9c6e4b1e97fb6e9d35a211294bb8bccb637f7af7c9f1`
Economics configuration hash: `57ce96656c71ac1d2e1a4e8fbc77578d9e292ed02737e926529ed97c80e4ce67`
Evidence registry hash: `72a9410cff2f41b7b252e28a87527d3c9683cfa810b2a83179fa460635f14e43`
