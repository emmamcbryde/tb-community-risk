# Appendix C. Using the Streamlit and Python APY model

## 1. Purpose and scope

The Streamlit APY workflow supports public-health planning for LTBI screening in the APY Lands. It helps users compare screening strategies, update epidemiological assumptions, run the APY model, review epidemiological outputs, optionally run economics, and download a consolidated Excel workbook.

The model is for planning and sequencing screening activity. It must not be used to deny screening, clinical review, or treatment to an individual person.

## 2. How to open the Streamlit application

In a deployed environment, open the Streamlit application URL provided for the APY model. The main workflow is through the Streamlit pages, not through editing MATLAB scripts or JSON files.

Technical local note: in a local development checkout, the application is started with Streamlit in the usual way for this repository. This appendix is not intended as an installation manual.

## 3. Selecting the Python APY backend

Open the `Scenario` page and use the `APY backend` selector.

Choose `Python APY v9 port (experimental)` to run the Python APY backend without MATLAB. The page continues to label this backend as experimental. MATLAB remains the reference implementation for parity testing and validation provenance.

## 4. Loading the APY defaults

On the `Scenario` page, select `Load backend defaults`. This loads the current APY default configuration, including the default cohort size, screening settings, test and regimen assumptions, calibration targets, and risk-factor inputs from `default_data.csv` and `default_age_distribution.csv` where applicable.

## 5. Changing the assumed prevalence of MTB infection/LTBI

In `Epidemiological assumptions`, use `Use default APY MTB infection prevalence` to choose whether the APY default is used.

When this box is selected, the internal configuration stores `ltbiPrevalence = None`, which tells the Python APY model to use the default APY calibration target.

When this box is deselected, enter `Assumed prevalence of MTB infection (LTBI)` as a percentage. The interface displays a percentage, but the model configuration stores a fraction. For example, 1% is stored as `0.01`.

This value is a calibration target for the simulated cohort. It is not the observed yield from the screening programme. Changing it recalibrates the probability of infection and can change expected test yield, false-positive burden, infections cured/protected, and active TB cases prevented.

## 6. Changing active-TB prevalence if required

The same section includes `Use default APY active-TB prevalence`. When selected, the internal configuration stores `activeTBPrevalence = None` and uses the APY default active-TB calibration target.

When deselected, enter the active-TB prevalence calibration target as a percentage. This is also stored internally as a fraction.

Custom active-TB prevalence should be used cautiously and reviewed epidemiologically, especially when paired with very low LTBI prevalence assumptions.

## 7. Optionally changing risk-factor prevalences

Open `Risk-factor prevalence overrides` on the `Scenario` page.

For close contact, marijuana use, renal disease, diabetes, smoking, chronic lung disease, and harmful alcohol/other drug use, the default option uses values from `default_data.csv`.

For each factor, users can select either a single overall prevalence or three age-group prevalences: 0-4 years, 5-14 years, and >=15 years. Values are entered as percentages and stored as fractions.

Female sex and BCG vaccination are available as advanced overall-prevalence overrides.

To return to model defaults, set the factor source back to `Use default`.

## 8. Selecting screening coverage, test, regimen and targeting strategy

The Scenario page exposes the routine planning controls:

- Scenario label.
- Screening coverage.
- Test type: IGRA or TST.
- Regimen.
- Screening strategy.
- Cohort size.
- Number of replicates.
- Random seed.
- Screening window.
- Follow-up horizon.

Routine public-health scenarios should be configured through these controls rather than by editing the `Full config` JSON.

## 9. Validating and running the model

After editing a scenario, use `Validate current config` on the Scenario page or `Validate scenario` on the Run Model page. Validation errors are shown with user-facing labels and the underlying config field.

On the `Run Model` page, select `Run Python APY model` or `Validate and run Python APY model`. The Python backend runs without MATLAB.

## 10. Recognising stale results

If epidemiological assumptions or run controls are changed after a model run, the app marks existing results as stale. Stale results should not be used as current outputs.

When results are stale, the Results page warns the user and disables the consolidated Excel workbook download until the model is rerun.

## 11. Interpreting headline epidemiological outputs

The Results page shows key metrics and summary results from the APY results bundle. Common outputs include people screened, test-positive non-active people, TPT starts and completions, false-positive treatment starts, infections cured/protected, active TB cases prevented, and number needed to screen or treat.

Active TB cases prevented are direct person-level effects in the simulated cohort. Downstream transmission effects are not included.

## 12. Running economics if required

Economics can be run from the `Economics` page after an APY model run. The economics workflow uses the selected backend outputs and current economics assumptions.

Economics is optional. If economics has not been run, the Results workbook includes an Economics sheet stating `Economics not run`; it does not replace missing economics outputs with zero.

## 13. Downloading and understanding the Excel workbook

On the Results page, after a current model run, select `Download consolidated APY results workbook`.

The workbook includes:

- `Read_me`: scenario name, purpose, interpretation guardrail, direct-effect limitation, backend warning, and stale-result status.
- `Scenario_inputs`: backend, scenario label, cohort size, replicates, seed, screening settings, strategy, test, regimen, LTBI and active-TB prevalence inputs, and treatment-cascade assumptions.
- `Risk_factor_inputs`: default or custom source for each risk factor and any entered prevalence values.
- `Headline_results`: key metric rows with human-readable labels.
- `Summary_results`: available model summary rows.
- `Natural_history`: do-nothing and intervention comparison fields where available.
- `Technical_metadata`: backend, model version, contract version, Git provenance where available, calibration information, source-data filenames, stale status, and current limitations.
- `Economics`: economics summary rows if economics has been run, or a clear `Economics not run` statement.

## 14. Saving a scenario for later reuse

Use `Save scenario` on the Scenario page. The saved JSON includes the current scenario config and economics config if available.

Custom LTBI prevalence, active-TB prevalence, and risk-factor overrides are preserved. Loading the scenario restores those values into the Streamlit controls.

## 15. Comparing alternative scenarios

Use the comparison workflow after running scenarios that represent the alternatives of interest. If the base scenario changes after comparison outputs are created, the app marks comparison outputs as stale.

## 16. Model limitations and appropriate use

The APY model is an individual-level planning model for screening strategy comparison. It is not a diagnostic tool, not an operational eligibility rule, and not a tool for denying care.

The Python backend currently supports the tested epidemiological and economics workflow. It remains experimental. It should not be described as fully validated or production-ready beyond the validation evidence in the repository.

The model reports direct person-level prevention. It does not include downstream transmission benefits in these APY outputs.

## 17. Python-MATLAB validation provenance

The current user workflow is Streamlit. Epidemiological and economics runs can be completed using the Python backend without MATLAB. MATLAB is retained as the reference implementation for parity testing. Progression attributable-risk add-ons remain outside the current Python backend.

Validation provenance is documented separately in the repository, including Python-versus-MATLAB distributional validation material. Exact individual-level equality is not expected because random number streams differ.

## 18. Troubleshooting

If a control is changed and results become stale, rerun the model before downloading or interpreting outputs.

If validation reports an input error, review the displayed user-facing field label and the config field shown beside it.

If a very low LTBI prevalence is paired with the default active-TB target and calibration fails, review whether the assumed LTBI and active-TB calibration targets are epidemiologically coherent.

If the Excel download is unavailable, check whether the current results are stale or whether the model has not yet been run.
