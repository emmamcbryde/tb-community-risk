# Country/subnational TB scientific-review application

An isolated **internal scientific-review workspace**, not a production application
or an externally validated country model. For planning and sequencing, not denying
care. The principal operational route is observed/user-specified TBI. Production
incidence inference remains blocked.

From the repository root, in the project's Python environment:

```powershell
python -m streamlit run country_model_review_ui/app.py --server.address 127.0.0.1 --server.enableCORS true --server.enableXsrfProtection true --browser.gatherUsageStats false
```

Open the localhost URL printed by Streamlit (normally `http://127.0.0.1:8501`).
An existing development instance may use another explicitly selected port.
The protective flags apply only to this local launch; no deployed configuration
or existing entry point is changed. The application works with the installed
Streamlit 1.32 API and the repository-pinned 1.39 runtime.

`app.py` lives below its own directory and uses a fixed seven-view allow-list in
`views/`. There is **no adjacent `pages/` directory**. Root legacy pages,
SA Health, explicit-TBI, no-transmission, Starsim and economics navigation are
not registered or imported. Do not launch this through an existing application.

## Quick review

1. Load national frozen demography or supply a custom/subnational population CSV.
   Only 2024 is supported by the retained baseline contracts. Female/male records
   stay separate. First UN WPP validation is a slower, explicit, cached operation.
2. Choose observed/user TBI. An overall total stays total-only unless age/recency
   derivation is explicitly selected. Direct age/sex and recent/remote profiles
   preserve supplied values. Blank values are not zero.
3. Supply risk prevalences/definitions, or explicitly choose an RR 1 reference.
   Missing/excluded/not-applicable/unreviewed states remain distinct. Legacy
   effects require an explicit diagnostic-sensitivity selection.
4. Inspect natural-history assumptions and gates. Adult source curves are
   sensitivity-only; pediatric curves and remote tail are not central defaults.
5. Run the deterministic review. Missing risk/progression inputs do not erase
   a valid TBI profile: they produce unresolved/blocked downstream outputs.
6. Inspect population-weighted results, uncertainties and assumptions. Material
   edits mark results stale. Prepare full-precision JSON/CSV exports explicitly.
7. Review the 3D evidence and unresolved decisions. Its empirical checks were
   internal synthesis/model checking, not external scientific validation.

**Restore working example** explicitly loads invented inputs: a custom adult
1000-person cohort,20% overall TBI, constant diagnostic pressure shape and RR 1.
The subgroup example is also invented. Neither is an Australian estimate.
No fictional example is used automatically by incidence inference.

## Uploads

See [CSV examples](examples/) and the full
[architecture/input specification](../docs/country_subnational_model/milestone3e_ui.md).
CSV probabilities are proportions 0–1; UI fields labelled `%` use 0–100. Sources
and review status are required. Blank upper age means open-ended, not an invented
terminal age. Population files require female/male counts and matching geography.
No unlogged interpolation, source relabelling or clipping occurs.

Use the Population page's custom ID `custom:synthetic_review` with the provided
example CSVs. The risk example is deliberately heavy smoking, not current smoking.
Direct total profiles have unavailable recency unless the explicit hybrid policy
is selected. A supplied direct 100+ total may resolve TBI prevalence but does not
supply an infection-time distribution or progression curve.

The incidence fork requires a separately uploaded, source-bound native review
bundle, explicit enabling, every scientific gate, then unique-root/refinement
clearance. UI selections or a hand-set session flag cannot authorize it. Actual
country configurations remain unapproved. See the [gate contract](../docs/country_subnational_model/milestone3e_ui.md#diagnostic-inversion-bundle).

## Validation and limitations

```powershell
python -m pytest tests/country_tb/test_review_ui_model.py tests/country_tb/test_review_ui_apptest.py -q --durations=8 -p no:tmpdir
python scripts/validate_country_tb_review_ui.py --output abm/output/country_tb_m3e_rendered --prospective
```

The second command records AppTest element/table evidence and timings, not browser
screenshots. Generated outputs stay under `abm/output/`. JSON retains native
80-digit decimals/exact rational pairs and binary 64 forward values; presentation
rounding never changes them. Export preparation is lazy.

Instantaneous rates and optional one-year first-episode projections are distinct.
The latter uses the unchanged 3B operator and a 50000-support-node preflight budget.
Oversized requests are blocked without dropping people. No long-horizon country
claim, stochastic UI, interventions, screening, treatment, economics, DALYs/QALYs
or dynamic transmission is implemented. See [validation](../docs/country_subnational_model/milestone3e_validation.md)
for exact results, runtime versions, browser availability and remaining limits.
