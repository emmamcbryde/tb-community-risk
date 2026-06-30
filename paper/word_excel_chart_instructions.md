# Using Editable Excel Charts In Word

The workbook `paper/excel/APY_economics_decision_tree_and_chart_data.xlsx` contains native Excel charts on the `Editable_Charts` sheet.

## Embedded Chart Workflow

1. Open the Excel workbook.
2. Go to `Editable_Charts`.
3. Select the chart.
4. Copy.
5. In Word, use `Home -> Paste -> Paste Special`.
6. Choose `Microsoft Excel Chart Object`.

Embedded charts travel with the Word document and can be edited by double-clicking the chart in Word. They do not update automatically if the Excel workbook later changes.

## Linked Chart Workflow

1. Open the Excel workbook.
2. Select and copy the chart.
3. In Word, use `Home -> Paste -> Paste Special`.
4. Choose `Paste Link -> Microsoft Excel Chart Object`.

Linked charts update when the source workbook changes, but the Word document depends on the workbook path remaining available. Use linked charts only if the workbook will be circulated with the report or stored in a stable shared location.

## Recommendation

For manuscript circulation, embedded Excel chart objects are usually safer. For active analysis drafts, linked charts are useful while tables and assumptions are still changing.
