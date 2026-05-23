# Manual Smoke Test: Python APY Backend

Use this smoke path for the normal Python-first APY workflow. MATLAB remains
the APY v9 reference engine, but it is only needed here for optional reference
checks or for generating fixtures.

Current conservative status:

- Python APY is the normal backend for this smoke test.
- Economics support is partial in Python; full parity with MATLAB economics is
  not assumed.
- Targeting and dynamic/ABM compare support are partial.
- Attributable-risk support and chart/export parity are partial.

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

7. Open Run Model and run the Python APY model. MATLAB is not required for this
   normal Python APY smoke run.

8. Open Results and confirm:

   - metadata shows backend `python`;
   - key metrics are displayed;
   - any attributable-risk, chart, or export outputs are treated as partial
     parity checks, not complete parity evidence.

9. Open Dynamic + ABM Compare and confirm it can read the latest APY bundle.
   Treat this as a partial compare smoke check.

10. Open Economics and confirm it explains that Python APY supports a partial
    economics subset; do not treat economics parity as complete.

11. Optional/reference only: return to Scenario and switch `APY backend` to
    `MATLAB v9 reference`.

12. If MATLAB is available, optionally confirm the MATLAB backend still loads
    defaults, validates, and runs as before, or use it to generate reference
    fixtures.
