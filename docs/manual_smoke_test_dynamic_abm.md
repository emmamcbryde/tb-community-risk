# Manual Smoke Test: Dynamic + APY ABM Integration

This smoke test checks the Streamlit integration between the Python dynamic model and the APY v9 MATLAB ABM workflow.

## Start the app

```powershell
streamlit run streamlit_app.py
```

## Dynamic workflow

1. Open **Dynamic Model > Dynamic Workflow**.
2. Run **1) Calibrate**.
   - If using dynamic calibration mode **Two-epoch beta: historical + recent 10 years**, confirm calibration completes with separate historical and recent transmission beta fits. The projection uses the recent beta.
   - This mode may be useful when recent incidence is falling, but it is not guaranteed to solve all declining-epidemic calibration issues.
3. Run **2) Simulate**.
4. Confirm the simulation table and projection chart render.
5. Confirm `st.session_state["dynamic_results_bundle"]` exists.
   - Use the **Dynamic + ABM Compare** page and open **Dynamic bundle debug**.
   - The debug panel should show `contractVersion = dynamic_output_contract_v9`, headline metrics, and projection annual row count.

## APY ABM workflow

1. Open **APY v9 ABM > Scenario**.
2. Click **Load backend defaults**.
3. Open **APY v9 ABM > Run Model**.
4. Run the APY model using the default Python APY backend. MATLAB is only required if the user selects and runs the MATLAB backend.
5. Confirm `st.session_state["results_bundle"]` exists.
   - Use the **Dynamic + ABM Compare** page and open **APY bundle debug**.
   - The debug panel should show top-level APY bundle keys, headline metrics, and metadata.
   - After a successful APY run with do-nothing/natural-history outputs included, confirm `technical.dynamicComparison.available = true`.
   - Confirm the APY debug panel lists dynamic-comparison metric names.

## Integrated comparison page

Open **Integrated workflow > Dynamic + ABM Compare** and check:

1. If the dynamic bundle is missing, the page warns that Dynamic Model must be calibrated and simulated.
2. If the APY bundle is missing, the page warns that APY Run Model must be run.
3. If both bundles are present, the page shows:
   - latest result summaries;
   - collapsed Dynamic/APY debug panels;
   - comparable metrics where both models expose numeric aligned values;
   - warnings for missing or structurally non-comparable metrics;
   - a comparison CSV download.
4. The page must not start MATLAB or require MATLAB unless the user returns to the APY Run Model page and runs the APY backend.

## Expected limitation

Some comparison rows may be unavailable until the APY results bundle exposes dynamic-style cumulative baseline/intervention active TB case metrics. This is expected and should be shown as a warning, not hidden.

## Confirmed APY dynamic-comparison block

The APY bundle debug panel should show `technical.dynamicComparison.available = true` when do-nothing/natural-history outputs are included. A confirmed MATLAB smoke test produced:

```text
technical.dynamicComparison.available = true
source = doNothing.derived
population = 1500
followHorizon = 20
cumulative_baseline_active_tb_cases = 37
cumulative_intervention_active_tb_cases = 25
cumulative_cases_averted = 12
relative_reduction_cumulative_active_tb_cases ~= 0.3196
```

Exact values may differ if APY scenario settings, seed, population, `nReps`, or follow-up horizon differ.
