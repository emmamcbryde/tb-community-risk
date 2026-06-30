# APY Economics Chart Fix Notes

## What Was Wrong

The previous workbook used formula-looking strings as scatter series titles, so Excel displayed legend labels such as `=Chart_Data!$A$2`. The cost-effectiveness plane charts also used one data point per strategy without reliably visible marker formatting, which could make the plot area appear blank.

## What Changed

- Rebuilt `Chart_Data` with explicit formula-linked columns A-J.
- Recreated all charts as native Excel charts linked to `Chart_Data`.
- Set scatter chart series names as plain text: Random, LTBI-targeted, Cure-targeted, Prevention-targeted.
- Set scatter charts to one visible marker series per strategy and removed connecting lines where supported.
- Added formula audit checks for chart data formulas, chart source notes, unresolved costs, and model-output rerun requirements.

## Editable Charts

The `Editable_Charts` sheet contains native Excel charts for cost per DALY, cost per QALY, net cost versus DALYs, net cost versus QALYs, active TB cases prevented, and break-even cost-saving margin. Copy them into Word using Paste Special as Microsoft Excel Chart Object or Paste Link as Microsoft Excel Chart Object.

## Remaining Limitations

Python writes formulas and chart definitions but does not calculate Excel's cached formula values. Open the workbook in Excel and allow recalculation before copying charts. Model-output values still require rerunning APY and refreshing `Model_Outputs`; Excel does not recreate the stochastic ABM.
