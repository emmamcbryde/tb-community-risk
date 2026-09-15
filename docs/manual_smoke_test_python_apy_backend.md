# Manual Smoke Test: APY Screening Workflow

This smoke test checks the standard non-dynamic APY screening workflow used
by the Streamlit application.

1. Start the app:

   ```powershell
   streamlit run streamlit_app.py
   ```

2. Open **Set up**.

3. Confirm the page has initialized the APY / SA Health defaults without a
   separate load button. Use **Restore APY defaults** only if you need to
   discard changes.

4. Review **Parameters**. Expand **View age distribution** and **View risk
   factors** if you need to check the resolved APY profile.

5. In **Analysis settings**, choose the analysis mode:

   - `Expected outcomes — single deterministic run` for a rapid exploratory
     calculation with no repetitions or random seed; or
   - `Simulated community variation — multiple stochastic runs` for repeated
     simulated populations.

6. For stochastic analysis, choose either:

   - `Quick preview: 100 repetitions` for a fast check; or
   - `SA Health reference: 2,000 repetitions` with seed `1` for the
     reference configuration.

   Repetitions and seed do not apply to deterministic expected-value runs.

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

10. Open **Run Analysis**. Confirm the page shows the selected analysis type
    and preview/reference status. Repetitions and seed should appear only for
    stochastic analysis. Run Analysis performs a final safety validation
    automatically when you click **Run analysis**; there is no separate
    validation step.

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
    analysis-mode, uncertainty and reproducibility details.
