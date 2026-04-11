# APY LTBI Screening Model (v9)

This folder contains the MATLAB-based APY screening and preventive-treatment model.

## Main entry point

The main user-facing file is:

`run_tb_screening_user_options_v9.m`

This is the file most users should open and edit.

## Core engine

The main model engine is:

`tb_screening_mc_model_v9.m`

This file contains the core simulation logic.

## Default input files

The default APY input files are:

- `default_data.csv`
- `default_age_distribution.csv` (or `.xlsx`, if still retained)

These provide the baseline parameter values used in the default scenario.

## Main output folder

All generated outputs should be written to:

`abm/output/`

This folder is intended for generated CSVs, figures, and temporary outputs, and should be excluded from Git tracking.

## Common scripts

### User-facing run
- `run_tb_screening_user_options_v9.m`

### Example run
- `run_tb_screening_example_v9.m`

### Main convenience workflow
- `v9main.m`

### Add-ons
- `run_tb_screening_compare_strategies_v9.m`
- `run_tb_screening_targeted_gradient_v9.m`
- `run_tb_screening_targeting_profile_v9.m`
- `run_tb_screening_targeting_optima_v9.m`
- `run_tb_screening_igra_charts_v9.m`
- `run_tb_screening_igra_addons_v9.m`
- `run_tb_screening_do_nothing_v9.m`
- `run_tb_screening_reactivation_attributable_v9.m`
- `run_tb_screening_natural_history_addons_v9.m`

## Typical workflow

1. Open MATLAB in the `abm/` folder
2. Open `run_tb_screening_user_options_v9.m`
3. Edit the USER INPUTS section as needed
4. Run the script
5. Review results and outputs in `abm/output/`

## Notes

- The v9 model is intended for planning and prioritisation.
- It is not intended to refuse screening or treatment to any individual.
- Targeting outputs should be interpreted as sequencing tools, not eligibility rules.