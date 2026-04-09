# AGENTS.md

## Project scope
This repo contains the APY v9 MATLAB ABM in `abm/`.

## Editing rules
- Treat `tb_screening_mc_model_v9.m` as the reference engine.
- Do not create new model versions unless explicitly asked.
- Prefer minimal edits over rewrites.
- Keep outputs in `abm/output/`.
- Do not hard-code local machine paths.
- Preserve compatibility with `run_tb_screening_user_options_v9.m`.

## Validation
After code changes:
- search for syntax issues
- check output paths
- if MATLAB is available, run lightweight validation commands only

## User-facing principle
This model is for planning and sequencing, not for denying care.