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

7. Optionally edit a low-risk field such as an economic unit cost, validate
   the parameters, and apply the change. Economic-only edits should not
   require rerunning screening outcomes.

8. If screening, treatment, demographic or epidemiological settings are
   changed, confirm the page marks previous results as stale.

9. Open **Run Analysis** and run the APY screening model.

10. Open **Results** and confirm headline screening and active-TB outcomes
    render.

11. Open **Health Economics**.

12. Confirm the assumptions workspace, annual costs, cost categories and
    summary results render without errors.

13. Recalculate health economics after an economic-only edit and confirm the
    screening outcomes are unchanged.

14. Download the workbook and confirm it opens.

15. Open **Evidence & Assumptions** and confirm unresolved evidence remains
    visible.
