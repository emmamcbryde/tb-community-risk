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

5. In **Analysis settings**, choose either:

   - `Quick preview: 100 repetitions` for a fast check; or
   - `SA Health reference: 2,000 repetitions` with seed `1` for the
     reference configuration.

6. Click **Review or change parameters**.

7. Confirm the parameter tabs are shown with:

   - `Parameter`
   - `Value used by model`
   - `Unit`
   - `Source`

8. Optionally edit a low-risk setup parameter. Economic-only assumptions are
   edited later on **Health Economics** and should not require rerunning
   screening outcomes.

9. If screening, treatment, demographic or epidemiological settings are
   changed, confirm the page marks previous results as stale.

10. Open **Run Analysis**. Confirm the page shows analysis type,
    repetitions, seed and preview/reference status, then run the analysis.

11. Open **Results** and confirm key metrics, the detailed summary and the
    per-100-person summary render. Use **Export results** for downloads.

12. Open **Health Economics**.

13. Confirm the page leads with analysis status, headline economic results,
    cost breakdown and economic scenario comparison.

14. Expand **View or change economic assumptions**, edit one low-risk cost
    in the `Value used by model` column, then click **Recalculate economics
    using current screening outcomes**. Confirm the override summary appears
    and screening outcomes are unchanged.

15. Download the assumptions or economics summary export and confirm it opens.

16. Open **Evidence & Assumptions** and confirm unresolved evidence remains
    visible. Expand **Caveats & technical information** for calibration,
    uncertainty and reproducibility details.
