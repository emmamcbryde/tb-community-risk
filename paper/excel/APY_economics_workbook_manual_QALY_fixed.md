# APY Economics Workbook Manual: QALY-Fixed Version

Workbook: `APY_economics_decision_tree_formula_linked_QALY_fixed.xlsx`

This workbook is for SA Health / APY Lands review of APY LTBI screening economics. It is a report and decision-support workbook. It does not change APY epidemiological model equations.

## Sheet Guide

| Sheet | Purpose |
|---|---|
| `README` | High-level workbook orientation and recalculation notes. |
| `Inputs_Assumptions` | Editable local costs, DALY inputs, QALY inputs, mortality assumptions, and willingness-to-pay thresholds. |
| `Model_Outputs` | Imported APY model outputs. These values should be changed only by rerunning/exporting APY outputs. |
| `Economic_Calculations` | Formula-linked economics, DALYs, mortality-adjusted QALYs, ICERs, NMB, and break-even calculations. |
| `Strategy_Summary` | Main review table for strategies, using Dale/Bauer mortality-adjusted QALYs as the primary QALY result. |
| `Cost_Saving_Check` | Cost-saving and break-even audit. Break-even TB cases are counts, not dollars. |
| `Decision_Tree` | Decision-pathway summary for screened/treated groups. |
| `Sensitivity` | Scenario/sensitivity rows from the report outputs. |
| `Chart_Data` | Formula-linked chart source data. Charts should reference this sheet rather than hard-coded values. |
| `Editable_Charts` | Native Excel charts that update after recalculation. |
| `Formula_Audit` | Formula and unit checks for key calculation links. |
| `Post_TB_Assumptions` | Editable post-TB sequelae assumptions and evidence notes. |
| `Post_TB_Scenarios` | Post-TB DALY/QALY/cost scenario calculations. |
| `Post_TB_Chart_Data` | Formula-linked post-TB chart source data. |
| `Post_TB_Audit` | Post-TB formula and unresolved-assumption checks. |

## Locally Determined Parameters

The following inputs require local review before final policy interpretation:

- Programme setup cost.
- Programme running cost.
- Outreach, travel, staff time, and follow-up support costs.
- Return-for-results costs.
- Clinical review and active TB exclusion costs.
- False-positive incremental costs.
- APY-specific TB case fatality risk.
- APY-specific years of life lost per TB death.
- Whether mortality-adjusted QALYs should use APY-specific life expectancy rather than the current YLL times healthy utility approach.
- APY-specific post-TB lung disease prevalence.
- Post-TB severity distribution, care costs, and duration.

## DALY Versus QALY Frameworks

DALYs and QALYs are different health-outcome frameworks and should not be mixed without clear labelling.

DALYs in this workbook:

- `yldAverted = nPreventedActiveTB * active TB disability weight * active TB duration`.
- `yllAverted = tbDeathsPrevented * YLL per TB death`.
- `dalysAverted = yldAverted + yllAverted`.

Primary QALYs in this workbook:

- Use Dale/Bauer compatible utility inputs for active TB morbidity.
- Include mortality benefits when DALYs include YLL.
- Use `QALYs gained, Dale/Bauer mortality-adjusted` as the primary QALY denominator for cost per QALY and QALY NMB.

GBD-aligned QALYs:

- Are retained as a sensitivity column.
- Use `yldAverted + mortalityQALYsGained - treatmentDecrementQALYLoss - SAEQALYLoss`.
- Should not be treated as the same construct as the Dale/Bauer utility-decrement analysis.

## Mortality-Adjusted QALY Logic

The key formulas in `Economic_Calculations` are:

- `tbDeathsPrevented = nPreventedActiveTB * TB case fatality risk`.
- `QALYs gained from deaths prevented = tbDeathsPrevented * YLL per TB death * healthy/LTBI utility`.
- `QALYs gained, Dale/Bauer morbidity-only = nPreventedActiveTB * (healthy/LTBI utility - active TB utility) * active TB duration`.
- `treatmentDecrementQALYLoss = nTotalCoursesStarted * LTBI treatment utility decrement * treatment duration`.
- `SAEQALYLoss = 0` in the current workbook unless a local SAE/ADR utility-loss linkage is specified.
- `QALYs gained, Dale/Bauer mortality-adjusted = Dale/Bauer morbidity-only QALYs + mortality QALYs - treatmentDecrementQALYLoss - SAEQALYLoss`.
- `QALYs gained, GBD-aligned mortality-adjusted = yldAverted + mortality QALYs - treatmentDecrementQALYLoss - SAEQALYLoss`.

Primary ICER and NMB formulas use Dale/Bauer mortality-adjusted QALYs:

- `costPerQALYGained = netCostVsBaseline / QALYs gained, Dale/Bauer mortality-adjusted`.
- `nmbQALY_45000 = QALYs gained, Dale/Bauer mortality-adjusted * WTP low - netCostVsBaseline`.
- `nmbQALY_75000 = QALYs gained, Dale/Bauer mortality-adjusted * WTP high - netCostVsBaseline`.

## Break-Even Calculations

Break-even TB cases are counts. Cost-saving margin columns are dollars.

- `breakEvenTBCasesPrevented = totalProgramCost / activeTBDiseaseCostPerCase`.
- `additionalCostMarginBeforeNotCostSaving = IF(netCostVsBaseline < 0, -netCostVsBaseline, 0)`.
- `additionalCostPerScreenedBeforeNotCostSaving = additionalCostMarginBeforeNotCostSaving / nScreened`.
- `additionalCostPerTPTStartBeforeNotCostSaving = additionalCostMarginBeforeNotCostSaving / nTotalCoursesStarted`.

Do not use break-even TB cases as a dollar value in charts.

## Double-Counting Risks

Avoid double counting when adding local costs or post-TB assumptions:

- Do not include active TB disease costs both as programme costs and as disease costs averted.
- Do not add post-TB care costs to programme delivery costs if they are already included in active TB/post-TB disease-cost assumptions.
- Do not combine DALYs and QALYs in a single ICER denominator.
- Do not use both Dale/Bauer utility decrement and GBD disability weight as primary QALY morbidity in the same result. One is primary; the other is sensitivity.
- Do not treat mortality QALYs and YLL as separate health gains in the same QALY denominator; mortality QALYs are the QALY translation of deaths prevented.

## Inputs That Require Rerunning APY

Changing the following requires rerunning the APY model or replacing `Model_Outputs`:

- Screening strategy.
- Test type if it changes modelled test performance/pathways rather than only unit cost.
- Regimen if it changes treatment effect, completion, adverse event, or duration assumptions.
- Cohort size or population/risk-factor inputs.
- LTBI prevalence, active TB prevalence, or risk-factor prevalence.
- Screening coverage, screening window, or follow-up horizon.
- Treatment start probability, completion probability, adverse reaction stop probability, and regimen efficacy.
- Any epidemiological parameter that changes active TB cases prevented or treatment-pathway counts.

Changing the following does not require rerunning APY, but does require workbook recalculation:

- Unit costs.
- Programme setup/running/delivery costs.
- Active TB disease cost per case.
- DALY disability weights and YLL assumptions.
- QALY utilities and WTP thresholds.
- Post-TB scenario assumptions.

## Review Notes

Open the workbook in Excel and allow recalculation before reviewing numeric outputs. Python writes formulas but does not calculate Excel cached formula values.
