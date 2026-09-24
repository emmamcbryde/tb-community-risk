# Population-profile schema (`population_profile_v2`)

Status: **implemented and validated**. Definitions:
`engine/profiles/population_profile.py`. Version 1 profiles are migrated on load.

## Value states

Every numeric input is a `ProfileValue` with separate fields:

* `value`: a number, or null;
* `state`: `value`, `missing`, `not_applicable` or `excluded`;
* `provenance`: `bundled` (demonstration), `who_snapshot`, `local_upload` or
  `user_defined`;
* `reviewStatus`: `not_reviewed`, `reviewed` or `not_required`;
* `unit`, `source` and `notes`.

A true zero is recorded as `state=value, value=0`, and a missing value as
`state=missing, value=null`. The two are never conflated.

## Top-level fields

| Field | Notes |
| --- | --- |
| `schemaVersion` | `population_profile_v2` |
| `profileId`, `name` | Country application prefixes the ISO3 code and country name |
| `location` | `name`, `kind` (`demonstration`, `country`, `territory`, `subnational`), `iso3`, `whoRegion`, `nationalPopulation`, `nationalPopulationYear`, `nationalPopulationSource`. The national population is context only, not the simulated population. |
| `populationSize` | Simulated population (default 10,000) |
| `ageDistribution`, `ageDistributionSource` | Broad age bands |
| `ltbiPrevalence` | Demonstration calibration target unless replaced |
| `incidence` | `source`, `snapshotId`, `provenance`, `measure` (`estimated_tb_disease_incidence`), `unit`, `dataYearRange`, `series` [{year, estimate, lower, upper}], `sourceDetail` (snapshot ID, report year, access date, data SHA-256, citation; or file name, file SHA-256, content hash), `dataHash`, `notes` |
| `trend` | `method` plus the trend settings (window, period, COVID handling, excluded years, seed, draws) and an optional user-specified annual change |
| `riskFactors` | See `docs/risk_factor_schema.md` |
| `dataVintage`, `notes`, `demonstration` | |

## Hashing and round-trip

* `profile_hash()` is the SHA-256 of the canonical, sorted, compact JSON.
* `PopulationProfile.from_json(p.to_json()) == p` holds for all profiles.
* `incidence.dataHash` is checked on load.

## What country data may change

Applying WHO country data changes only `location`, `incidence`, `profileId`,
`name`, `dataVintage` and `trend`. It never changes LTBI prevalence, age
distribution, test accuracy, the treatment cascade, costs, DALYs or risk factors.
User-supplied incidence or a subnational location can be replaced only after an
explicit choice.

## Link recorded with each analysis (`generalProfileLink`)

Each analysis records:

* the profile ID, name, schema and hash;
* the location and ISO3 code;
* the incidence provenance, snapshot ID, snapshot data SHA-256 and incidence
  data hash;
* the data vintage;
* `incidenceUsedByEngine: false`.
