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

5. In **Analysis settings**, choose the analysis mode:

   - `Deterministic expected-value analysis` for a rapid exploratory
     calculation with no repetitions or random seed; or
   - `Stochastic individual-based analysis` for repeated simulated
     populations.

6. For stochastic analysis, choose either:

   - `Quick preview: 100 repetitions` for a fast check; or
   - `SA Health reference: 2,000 repetitions` with seed `1` for the
     reference configuration.

   Repetitions and seed do not apply to deterministic expected-value runs.

7. Click **Review or change parameters**.

8. Confirm the parameter tabs are shown with:

   - `Parameter`
   - `Value used by model`
   - `Unit`
   - `Source`

9. Optionally edit a low-risk setup parameter. Economic-only assumptions are
   edited later on **Health Economics** and should not require rerunning
   screening outcomes.

10. If screening, treatment, demographic or epidemiological settings are
   changed, confirm the page marks previous results as stale.

11. Open **Run Analysis**. Confirm the page shows the selected analysis type
    and preview/reference status. Repetitions and seed should appear only for
    stochastic analysis. Run Analysis performs a final safety validation
    automatically when you click **Run analysis**; there is no separate
    validation step.

12. Open **Results** and confirm key metrics, the detailed summary and the
    per-100-person summary render. Use **Export results** for downloads.

13. Open **Health Economics**.

14. Confirm the page leads with analysis status, headline economic results,
    cost breakdown and economic scenario comparison.

15. Expand **View or change economic assumptions**, edit one low-risk cost
    in the `Value used by model` column, then click **Recalculate economics
    using current screening outcomes**. Confirm the override summary appears
    and screening outcomes are unchanged.

16. Download the assumptions or economics summary export and confirm it opens.

17. Open **Evidence & Assumptions** and confirm unresolved evidence remains
    visible. Expand **Caveats & technical information** for calibration,
    analysis-mode, uncertainty and reproducibility details.
