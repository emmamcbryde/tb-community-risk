# Immutable scientific validation inputs

These are source-preserving **internal validation snapshots**, not country
parameter defaults, empirical latent-TBI truth or externally validated estimates.
See [3D protocol](../../docs/country_subnational_model/milestone3d_protocol.md),
[evidence](../../docs/country_subnational_model/milestone3d_evidence.md) and
[results](../../docs/country_subnational_model/milestone3d_empirical.md).

`snapshots/tbi-validation-fe05d6d1ab8c/` preserves96 selected observations:

- Innes AL et al. Public Health Action2023,doi:10.5588/pha.23.0020,
  [PMC10446657](https://pmc.ncbi.nlm.nih.gov/articles/PMC10446657/),CC BY4.0.
- Miramontes R et al. PLoS One2015,doi:10.1371/journal.pone.0140881,
  [PMC4633161](https://pmc.ncbi.nlm.nih.gov/articles/PMC4633161/),public domain.
- Laemmle-Ruff I et al. Emerg Infect Dis2022,doi:10.3201/eid2809.212426,
  [PMC9423895](https://pmc.ncbi.nlm.nih.gov/articles/PMC9423895/),public-domain
  source with attribution retained.

Source licences/URLs, full citations, access date2026-10-30, table/cell locators,
raw strings and original XML SHA-256 are in `sources.json` and
`source_observations.json`. `model_inputs.json` contains exact count ratios or
published weighted percentages divided by100, never reconstructed counts.
`manifest.json` checksums all three files and records the scientific hash.
No values are digitized from plots. Counts, rounded percentages, survey-weighted
estimates, missingness, overlapping marginal views and source ambiguities are
preserved. The loader verifies bytes and rederives transformed inputs offline.

`mortality/abs-life2022-2024-6785449a5200/` contains selected ABS Table9Australia
values and the1-qx transformation. Attribution: Commonwealth of Australia,
Australian Bureau of Statistics, *Life expectancy2022–2024*;
[source](https://www.abs.gov.au/statistics/people/population/life-expectancy/2022-2024),
[CC BY4.0 rights](https://www.abs.gov.au/website-privacy-copyright-and-disclaimer).
The workbook checksum, original cells,202 raw rows and200 finite survival
intervals are retained. Raw age100 does not expand a demographic100+ cell.
All-cause-vs-competing-cause mortality and country adoption remain unapproved.
There is no Australian mortality fallback for another geography.

Re-extraction requires separately acquired primary files matching recorded
source hashes; ordinary validation never downloads them. From repository root:

```
python scripts/freeze_country_tb_validation_3d.py --primary abm/output/country_tb_m3d_audit_20261030/primary
python scripts/freeze_country_tb_mortality_3d.py --workbook abm/output/country_tb_m3d_audit_20261030/primary/ABS_lifetables_2022_2024.xlsx
python scripts/validate_country_tb_empirical_3d.py --output abm/output/country_tb_m3d_validation
```

The extraction scripts use standard-library XML parsing and read-only openpyxl
for the workbook. Identical reimport is a no-op; any changed payload fails rather
than replacing a snapshot. New source versions require new immutable directories.
Generated reports remain under `abm/output/`. Read canonical JSON explicitly as
UTF-8 on Windows. Production inference remains blocked.
