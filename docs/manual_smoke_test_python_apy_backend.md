# Manual Smoke Test: Experimental Python APY Backend

MATLAB remains the APY v9 reference backend. The Python APY backend is
experimental and currently supports a partial economics subset while excluding
attributable-risk add-ons.

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

10. Open Economics and confirm it explains that Python APY supports a partial
    economics subset and that the full economics model remains MATLAB-backed.

11. Return to Scenario and switch `APY backend` back to `MATLAB v9 reference`.

12. If MATLAB is available, confirm the MATLAB backend still loads defaults,
    validates, and runs as before.
