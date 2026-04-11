# AGENTS.md — APY MATLAB ABM (v9)

This file gives repo-specific instructions for AI coding agents working inside `abm/`.

## Scope

This folder contains the MATLAB-based APY LTBI screening and preventive-treatment model.

Work in this folder should focus on:
- the v9 APY model
- user-facing usability improvements
- safe input handling
- output management
- reporting / chart generation
- maintenance of the current validated workflow

Do **not** treat older model versions as the active codebase unless explicitly asked.

---

## Canonical starting point

### Main user-facing entry point
`run_tb_screening_user_options_v9.m`

This is the main file a real user should open and edit.

### Core engine
`tb_screening_mc_model_v9.m`

This is the reference model engine.

### Convenience / orchestration script
`v9main.m`

This is an analyst helper script, not the canonical minimal user entry point.

---

## Default data files

The default APY data files are:

- `default_data.csv`
- `default_age_distribution.csv`

If backward compatibility is still needed, `.xlsx` support may remain temporarily, but CSV is preferred as the canonical format.

Do not hard-code personal machine paths.

---

## Output handling

All generated outputs should go into:

`abm/output/`

Rules:
- do not write generated CSVs or PNGs into the main `abm/` folder
- preserve `.gitignore` behaviour for `abm/output/`
- create the output directory if needed
- keep filenames stable unless explicitly asked to rename them

---

## Editing rules

When editing this codebase:

1. Prefer **minimal edits** over rewrites.
2. Preserve compatibility with:
   - `run_tb_screening_user_options_v9.m`
   - current v9 add-on scripts
3. Do not create new model versions unless explicitly asked.
4. Do not rename the core engine or user-facing runner unless explicitly asked.
5. Keep APY defaults intact unless the requested task is specifically to revise them.
6. If changing data filenames, update all references consistently.
7. If changing output locations, update all writing scripts consistently.
8. Preserve the planning purpose of the model:
   - this model supports sequencing and prioritisation
   - it should not be framed as a gatekeeping tool for denying care

---

## Validation expectations

After making code changes:

- search for remaining hard-coded filenames or paths
- search for output-writing calls (`writetable`, `saveas`, `exportgraphics`, etc.)
- confirm outputs still go to `abm/output/`
- check that user-facing inputs still flow into the engine correctly
- if MATLAB is available, prefer lightweight validation commands only
- do not run expensive long simulations unless explicitly asked

Suggested lightweight checks:
- inspect the canonical entry point
- inspect output paths
- run a small example only if requested

---

## User-facing design intent

For APY v9, keep the tool as user-friendly as possible.

### Simple mode should prioritise:
- LTBI prevalence
- active TB prevalence
- age distribution
- risk-factor prevalence
- test choice
- treatment choice
- cascade of care
- screening strategy
- screening coverage

### Advanced mode can include:
- exact-age inputs
- disease OR overrides
- test-performance overrides
- calibration details
- seed / replicate count
- alternative scenario controls

If asked to simplify the interface, hide complexity rather than removing capability.

---

## Documentation expectations

If user-facing workflow changes, update the relevant README or help text.

Important distinction:
- `README.md` = for users
- `AGENTS.md` = for coding agents

If a new script becomes the main entry point, document that explicitly.

---

## Non-goals unless explicitly requested

Do not:
- merge this MATLAB model into the Streamlit system
- redesign the epidemiologic logic
- remove archived or older versions permanently
- change default epidemiologic assumptions for reasons not requested by the user
- add broad new features just because they seem useful

---

## Preferred working style

When given a task:
1. identify the smallest set of files that need to change
2. inspect first
3. summarise what will be changed
4. edit conservatively
5. validate
6. summarise the exact files changed


## UI mode rules

The APY v9 interface uses two modes:
- Simple mode
- Advanced mode

Do not infer what these mean from context.
Follow `SIMPLE_ADVANCED_SPEC.md` as the source of truth.

If a field is not explicitly listed for Simple mode, put it in Advanced mode unless told otherwise.