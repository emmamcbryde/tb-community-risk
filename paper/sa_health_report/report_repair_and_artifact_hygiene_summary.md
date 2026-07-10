# Report Repair and Artifact Hygiene Summary

## Broken file

- `paper/sa_health_report/APY_SA_Health_LTBI_model_report_draft.docx`
- Size: 293,497 bytes
- ZIP status: valid ZIP
- `ZipFile.testzip()`: `None`
- `word/document.xml`: present
- XML parse result: failed with `ParseError: unbound prefix: line 1, column 60306`
- `python-docx`: not available in the current Python environment, so the document could not be opened/resaved with `python-docx`
- LibreOffice/soffice/Word command-line conversion: not available on PATH, so PDF conversion was not run

## Likely cause

The generated report package reused section/page XML from a source Word document while writing a minimal `document.xml` root. The copied section XML included namespace-prefixed attributes that were not declared on the new minimal root element. Microsoft Word rejects this as invalid WordprocessingML.

The issue does not appear to be caused by inserted figures or charts. The rebuilt report keeps figure content as highlighted placeholders and keeps Excel chart workbooks separate.

## Last known-good source used

- `APY_LTBI_health_economics_report_updated_post_tb.docx`
- Size: 303,167 bytes
- ZIP status: valid ZIP
- `ZipFile.testzip()`: `None`
- `word/document.xml`: present
- XML parse result: OK

## Rebuilt file

- `paper/sa_health_report/generated/APY_SA_Health_LTBI_model_report_REBUILT.docx`
- Size: 293,559 bytes
- ZIP status: valid ZIP
- `ZipFile.testzip()`: `None`
- `word/document.xml`: present
- XML parse result: OK

The corrupt generated draft was copied, not moved or deleted, to:

- `paper/sa_health_report/archive/APY_SA_Health_LTBI_model_report_draft.docx`

## Figure handling

No figure objects were inserted into the rebuilt Word report. The report uses yellow-highlighted placeholders for figures and keeps chart/workbook artefacts separate under:

- `paper/sa_health_report/figures/`
- `paper/sa_health_report/excel/`

Recommended final layout workflow: open the formula-linked Excel workbook, copy final charts manually, and paste into Word as editable Excel chart objects after content review.

## Ignore-rule changes

`.gitignore` was updated to ignore generated report artefacts:

- `/paper/sa_health_report/generated/`
- `/paper/sa_health_report/archive/`
- `/paper/sa_health_report/*.docx`
- `/paper/sa_health_report/*.pdf`
- `/paper/sa_health_report/*.tmp`
- `/paper/sa_health_report/figures/`
- `/paper/sa_health_report/tables/`
- `/paper/sa_health_report/excel/`
- `/paper/excel/*.xlsx`
- `/paper/figures/`
- `/outputs/`

Markdown documentation remains intended to be trackable via:

- `!/paper/**/README.md`
- `!/paper/**/*.md`

## Generated files already tracked

The following generated artefacts are already tracked and will not be affected by `.gitignore` until explicitly removed from the index:

- `outputs/apy_paper_daly_icer/apy_cost_outlay_audit.csv`
- `outputs/apy_paper_daly_icer/apy_cost_outlay_audit_notes.md`
- `outputs/apy_paper_daly_icer/apy_cost_saving_break_even.csv`
- `outputs/apy_paper_daly_icer/apy_daly_icer_costs.csv`
- `outputs/apy_paper_daly_icer/apy_daly_icer_health_outcomes.csv`
- `outputs/apy_paper_daly_icer/apy_daly_icer_notes.md`
- `outputs/apy_paper_daly_icer/apy_daly_icer_sensitivity.csv`
- `outputs/apy_paper_daly_icer/apy_daly_icer_strategy_summary.csv`
- `outputs/apy_paper_daly_icer/apy_daly_qaly_source_audit.csv`
- `outputs/apy_paper_daly_icer/apy_daly_qaly_source_audit_notes.md`
- `outputs/apy_paper_daly_icer/apy_post_tb_sequelae_notes.md`
- `outputs/apy_paper_daly_icer/apy_post_tb_sequelae_scenarios.csv`
- `outputs/apy_paper_daly_icer_quick/apy_daly_icer_costs.csv`
- `outputs/apy_paper_daly_icer_quick/apy_daly_icer_health_outcomes.csv`
- `outputs/apy_paper_daly_icer_quick/apy_daly_icer_notes.md`
- `outputs/apy_paper_daly_icer_quick/apy_daly_icer_sensitivity.csv`
- `outputs/apy_paper_daly_icer_quick/apy_daly_icer_strategy_summary.csv`
- `paper/excel/APY_economics_decision_tree_and_chart_data.xlsx`
- `paper/excel/APY_economics_decision_tree_formula_linked.xlsx`
- `paper/excel/APY_economics_decision_tree_formula_linked_post_tb.xlsx`
- `paper/excel/APY_economics_decision_tree_formula_linked_v3.xlsx`
- `paper/excel/APY_economics_decision_tree_formula_linked_v4.xlsx`
- `paper/excel/APY_economics_decision_tree_formula_linked_v5.xlsx`
- `paper/figures/cost_per_daly_averted_by_strategy.png`
- `paper/figures/cost_per_qaly_gained_by_strategy.png`
- `paper/figures/net_cost_vs_dalys_averted.png`
- `paper/figures/net_cost_vs_qalys_gained.png`
- `paper/figures/prevent_strategy_daly_sensitivity.png`

Tracked Markdown files in `paper/` may reasonably remain tracked:

- `paper/excel/APY_economics_chart_fix_notes.md`
- `paper/excel/APY_economics_formula_linking_notes.md`
- `paper/word_excel_chart_instructions.md`
- `paper/word_report_cost_saving_explanation.md`
- `paper/word_report_update_summary.md`

## Recommended index cleanup commands

Do not run these without explicit approval:

```powershell
git rm --cached -r outputs
git rm --cached paper/excel/*.xlsx
git rm --cached -r paper/figures
```

After index cleanup, keep lightweight Markdown documentation tracked.
