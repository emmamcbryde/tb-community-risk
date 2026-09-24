# General community TB application - Milestone 1

Branch: `feature/general-community-tb-model` (based on the frozen release tag
`sa-health-apy-he-v1.0.0`, commit `03cc16e4a52a55e10019dab16bd4f599571c4b90`).

This milestone adds the architecture and first interface for a setting-neutral
community TB screening decision-support application. It does not change the
validated engines, economic or DALY formulas, the dynamic model, the MATLAB
reference files or the frozen release artifacts.

## Entry points

| Entry point | Purpose |
| --- | --- |
| `general_app.py` | General application (new). Six-step workflow in `general_pages/`. |
| `streamlit_app.py` | Existing setting-specific workflow and research pages (unchanged). |
| Tag `sa-health-apy-he-v1.0.0` | Frozen release; deployed separately and never modified. |

General workflow: Set up population -> Define intervention -> Run analysis ->
Results -> Health economics -> Evidence and technical information.

## General defaults

* Profile label: **Demonstration working defaults**, with a warning that inputs are
  not evidence for any particular country.
* Population size 10,000; stochastic analysis with 1,000 simulated populations;
  deterministic expected-value preview available.
* Screening test (IGRA), regimen (3HP), coverage, targeting and every other model
  value are taken unchanged from the validated working defaults.
* One action, **Restore demonstration defaults**, resets the profile,
  intervention and analysis settings and clears all user-defined values.
* Run time: roughly 0.8-0.9 s per simulated population of 10,000 on a typical
  laptop, so the default 1,000-simulation analysis takes about 15 minutes.

## Population-profile contract (`engine/profiles/population_profile.py`)

Schema `population_profile_v1`; JSON round-trip and deterministic SHA-256 hashing
(canonical sorted JSON). Fields: profile ID and name; location (name, kind,
ISO3); population size; age distribution (bands); LTBI prevalence; incidence
(source, snapshot ID, provenance, series of year/estimate/lower/upper, derived
data year range, measure = estimated TB disease incidence); trend (method,
settings, optional user annual percentage change); risk factors; data vintage;
notes; demonstration flag.

Every numeric input is a `ProfileValue` with separate `state`
(`value` / `missing` / `not_applicable` / `excluded`), `provenance`
(`bundled` / `user_defined`) and `reviewStatus`
(`not_reviewed` / `reviewed` / `not_required`). Zero is `state=value, value=0`;
missing is `state=missing, value=null`; these are never conflated.

Risk factors carry: enabled, prevalence, effect estimate, effect-measure type
(RR, HR or OR, retained exactly), unit, evidence source, review status, notes,
the engine key (if any) and the list of user-modified fields. A profile with no
risk factors is valid and runnable.

### Engine mapping (`engine/profiles/engine_mapping.py`)

* Bundled values are passed to the engine as "use engine default", so an
  unchanged demonstration profile reproduces the validated engine inputs exactly
  (tested).
* User-defined prevalence is applied uniformly across age groups (existing engine
  override behaviour).
* **Internal assumption:** the engine multiplies the hazard of progression from
  infection to disease by the effect estimate, whatever its declared type. ORs and
  RRs are not converted. This is shown on the Evidence and technical information page.
* Disabled or excluded risk factors, and profiles without risk factors, run with
  prevalence 0 for that factor (no stratification).
* The attached incidence series is descriptive only in this milestone.
* Each run records `generalProfileLink` (profile ID, schema, hash, incidence
  snapshot ID and data checksum), so saved results identify their exact inputs.

## Data boundary

```
gtbreport2025 checkout / WHO CSV  --(maintainer runs scripts/import_who_incidence.py)-->
    data/who_incidence/<snapshot>.csv + <snapshot>_manifest.json   (versioned, checksummed)
        --> engine/who_incidence/snapshot.py (offline load + checksum + validation)
            --> Streamlit general application (no network)
```

The application never reads the neighbouring repository, calls WHO services or
imports external code at runtime. Tests do not need MATLAB, internet access or a
`gtbreport2025` checkout.

Manifest fields: schema version, snapshot ID and kind, completeness flag, source
dataset, report year, source and public download URLs, upstream repository,
commit and source-file checksums, extraction date, transformation version, column
list, data file name/checksum/row count, countries, year range, validation
status and warnings, licence status and citation.

Validation (`engine/who_incidence/schema.py`) catches duplicate country-year rows,
invalid ISO3 codes, invalid years, non-numeric and negative values, lower bounds
above estimates, upper bounds below estimates, missing years (warning),
implausible year-on-year discontinuities (warning, ratio > 2), schema changes and
incompatible data vintages.

The bundled snapshot `who-gtb2025-fixture-v1` is a **test fixture** of six
countries (AUS, BRA, GBR, IDN, PHL, ZAF; 2000-2024; 150 rows) extracted from
`inc_mort/analysis/est.rda` at gtbreport2025 commit `666088c`. It is not the
complete WHO dataset, and the interface says so.

## Audit summary: gtbreport2025

* Fork `emmamcbryde/gtbreport2025` of upstream `GTB-TME/gtbreport2025`; `main` at
  `666088c` ("First commit of the public version of the 2025 global TB report
  repository", 12 Nov 2025). Report year 2025; estimates cover 2000-2024.
* No licence file in the repository. WHO states that data use is subject to WHO's
  terms and conditions and data policy. The licence status is recorded as
  `to_be_confirmed`. Only estimate values are extracted; no code is copied.
* Contents mix source data (R `.rda` extracts of the WHO database snapshot
  2025-07-30, UN population, surveys, external indicators), generated outputs
  (`inc_mort/analysis/*.rda`, `disaggregation/output/*`), estimation code (R, Stan,
  MATLAB) and report code (R Markdown).
* Country incidence: `inc_mort/analysis/est.rda` (3.4 MB, 5,347 rows x 203
  columns, 217 ISO3 countries/areas, 2000-2024). Fields `inc`, `inc.lo`, `inc.hi`,
  `inc.sd` (per 100,000), `inc.num`/`inc.lo.num`/`inc.hi.num`, `pop`,
  `source.inc`, HIV-split `inc.h*`/`inc.nh*`, and method intermediates.
* Population: `data/gtb/other/pop.rda` (UNPD by sex and age band, 1950-2050).
  Country names/codes: `data/gtb/other/cty.rda` (ISO2, ISO3, numeric, WHO region).
* Risk factors: `inc_mort/analysis/attributable_cases.rda` (prevalence, PAF and
  attributable incidence for HIV, diabetes, alcohol, smoking, undernutrition) and
  `rf.rda` (attributable incident numbers). Age/sex disaggregation:
  `disaggregation/output/db_estimates_country_all.Rdata`.
* Data dictionary: `dic.rda` (WHO database dictionary) plus `doc/who_tb_dcf2025.pdf`.
* Reproducibility: the pipeline depends on `import/load_gtb.R` and WHO database
  access, which are not in the public repository. Generated outputs can be read
  without credentials, but they cannot be regenerated from the public repository.
* Vintage: `est.rda` is the full retrospective series from the 2025 round. WHO
  revises all past years each round (`estimates2024/old.rda` holds the 2024-round
  series for comparison), so vintages must not be mixed.
* Smallest authoritative inputs for country incidence: either the public WHO CSV
  `https://extranet.who.int/tme/generateCSV.asp?ds=estimates` (plus
  `ds=dictionary`), or `inc_mort/analysis/est.rda` + `data/gtb/other/cty.rda` at a
  pinned commit.

## Next milestone: WHO import

1. Obtain the licence position for redistributing WHO estimate extracts. Decide
   whether the production snapshot stores all countries or is built on demand by
   maintainers.
2. Choose the production source. The recommendation is the public WHO CSV, as the
   current-round series, with `est.rda` at a pinned commit as a cross-check. Record
   both checksums.
3. Run `scripts/import_who_incidence.py who-csv --kind production`. Commit the
   snapshot under `data/who_incidence/snapshots/`, then add a selector for
   multiple bundled vintages. Vintages are never combined.
4. Add ISO3 validation against the WHO country list (`cty`) and a regression
   comparison between the CSV and `est.rda` values for the same round.
5. Extend profiles with UNPD age distributions (optional), clearly labelled.

## Scientific decisions still required

* Which trend method (log-linear recent, penalised spline or state-space), which
  fitting period, and how to handle 2020-2022 COVID-era disruption.
* How WHO bounds are propagated. They are not symmetric and not a posterior
  distribution.
* The linkage model from estimated disease incidence to infection pressure or
  force of infection. It must account for progression, risk-factor prevalence
  and case detection. The disease-incidence slope is not to be used as the
  infection-pressure slope.
* How country incidence is used in calibration of the individual-based engine
  (currently calibrated to bundled demonstration targets) and in the dynamic
  model.
* Whether ORs from the literature should be converted before being applied as
  progression hazard multipliers, and by which method (baseline-risk dependent).
* Country-specific LTBI prevalence, age distribution, test accuracy, treatment
  cascade and cost inputs, and their evidence review workflow.
