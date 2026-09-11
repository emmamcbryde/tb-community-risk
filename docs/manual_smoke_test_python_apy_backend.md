# Manual Smoke Test: APY Screening Workflow

This smoke test checks the standard non-dynamic APY screening workflow used
by the Streamlit application.

1. Start the app:

   ```powershell
   streamlit run streamlit_app.py
   ```

2. Open **Set up**.

3. Click **Use default parameters**.

4. Review **Current working defaults** and expand **Age distribution and risk
   factors**.

5. Click **Review or change parameters**.

6. Confirm the parameter tabs are shown with:

   - `Parameter`
   - `Value used by model`
   - `Unit`
   - `Source`

7. Optionally edit a low-risk setup parameter. Economic-only assumptions are
   edited later on **Health Economics** and should not require rerunning
   screening outcomes.

8. If screening, treatment, demographic or epidemiological settings are
   changed, confirm the page marks previous results as stale.

9. Open **Run Analysis** and run the APY screening model.

10. Open **Results** and confirm headline screening and active-TB outcomes
    render.

11. Open **Health Economics**.

12. Confirm the page leads with analysis status, headline economic results,
    cost breakdown and economic scenario comparison.

13. Expand **View or change economic assumptions**, edit one low-risk cost
    in the `Value used by model` column, then click **Recalculate economics
    using current screening outcomes**. Confirm the override summary appears
    and screening outcomes are unchanged.

14. Download the assumptions or economics summary export and confirm it opens.

15. Open **Evidence & Assumptions** and confirm unresolved evidence remains
    visible.
