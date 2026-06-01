# Manual Smoke Test: Experimental Python APY Backend

MATLAB remains the APY v9 reference backend. The Python APY backend is
experimental. Python economics is newly ported; attributable-risk add-ons
remain outside the Python backend.

1. Start the app:

   ```powershell
   streamlit run streamlit_app.py
   ```

2. Open the Scenario page.

3. Set `APY backend` to `Python APY v9 port (experimental)`.

4. Click `Load backend defaults`.

5. In `Advanced / run controls`, set a small smoke-test run:

   - `N = 100` or `500`
   - `Simulation replicates = 5` or `20`
   - `Random seed = 1`

6. Validate the scenario.

7. Open Run Model and run the Python APY model.

8. Open Results and confirm:

   - metadata shows backend `python`;
   - key metrics are displayed;
   - `technical.dynamicComparison` is available from `doNothing.derived`.

9. Open Dynamic + ABM Compare and confirm it can read the latest APY bundle.

10. Open Economics.

11. Click `Load economics defaults` and confirm blank/optional assumptions
    render without errors.

12. Click `Load KWAB150 preset`.

13. Confirm AUD 2019 Australia values populate the visible economics controls.

14. Click `Run economics`.

15. Confirm:

    - economics summary rows render;
    - economics status renders;
    - economics summary CSV download appears;
    - no MATLAB session is required for the Python backend economics run.

16. Return to Scenario and switch `APY backend` back to `MATLAB v9 reference`.

17. If MATLAB is available, confirm the MATLAB backend still loads defaults,
    validates, and runs as before.
