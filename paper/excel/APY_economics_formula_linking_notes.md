# APY Economics Formula Linking Notes

Input workbook inspected: `C:\Users\emmas\Documents\GITHUB\tb-community-risk\APY_economics_decision_tree_and_editable_charts_v2.xlsx`. The requested `paper/excel/APY_economics_decision_tree_and_editable_charts_v2.xlsx` path was used if present; otherwise the same-named workbook in the repository root was used.

## Sheets Updated

- `Model_Outputs`: static imported APY outputs from strategy, health, costs, and sensitivity CSVs.
- `Economic_Calculations`: formula-linked economics, DALYs, QALYs, ICERs, NMB, and break-even calculations.
- `Strategy_Summary`: readable formula-linked summary table.
- `Chart_Data`: formula-linked chart source data.
- `Editable_Charts`: native Excel charts pointing to `Chart_Data`.
- `Formula_Audit`: formula checks and notes on assumptions requiring APY reruns.

## Imported APY Model Outputs

The following remain imported stochastic model outputs and are not recalculated in Excel: strategy, N, nReps, screen coverage, screen window, follow-up horizon, screened count, TPT starts/completions, false-positive treated count, infections cured/protected, active TB cases prevented, baseline/intervention active TB by horizon, and relative reduction. To change these values, rerun APY and refresh `Model_Outputs`.

## Formula-Linked Values

Costs, DALYs, QALYs, ICERs, net monetary benefit, break-even TB cases, additional cost margin, and chart values are formula-linked to `Inputs_Assumptions` plus imported model outputs.

## Assumptions That Update Charts

Test costs, regimen costs, active TB disease cost, programme setup/running, false-positive incremental cost, return-for-results cost, clinic review cost, active TB exclusion/CXR/lab cost, travel/outreach/staff support cost, DALY assumptions, utility assumptions, LTBI treatment utility decrement, and WTP thresholds update `Economic_Calculations`, `Strategy_Summary`, `Chart_Data`, and charts.

## Formula Status

`Strategy_Summary` and `Chart_Data` are formula-linked. Open the workbook in Excel and allow recalculation; Python writes formulas but does not calculate cached values.
