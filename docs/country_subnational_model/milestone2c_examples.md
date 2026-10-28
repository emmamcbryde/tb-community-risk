# Milestone 2C diagnostic results

These are **diagnostic software scenarios, not country TBI estimates**.
No real TBI anchor or assay-performance value is adopted. The country example
uses the same invented supplied alpha **0.005/year**, common to both sexes,
2024.5 baseline time, a five-year recent window, the frozen WHO 2025 incidence
snapshot, and frozen ABS/WPP 2024 demographic weights.

Run `python scripts/preview_country_tb_calibration.py`. Exports are immutable,
hash-identified folders under `abm/output/country_tb_m2c_preview` by default.
This recovery used `abm/output/country_tb_m2c_reinfection_preview`; intermediate
exports remain preserved. Current scenario metadata is additionally published
in `final-*` folders in that directory. No output replaces frozen inputs.

## Central diagnostic profiles and unresolved ages

| Quantity | Australia | Viet Nam |
| --- | ---: | ---: |
| Selected demography | ABS June 2024 estimate | WPP 2024 Medium projection |
| Female population | 13,698,889 | 51,526,284 |
| Male population | 13,505,920 | 49,461,398 |
| Complete derived population | 27,204,809 | 100,987,682 |
| Age-covered population | 27,198,221 | 100,980,543 |
| Unestimated 100+ population | 6,588 | 7,139 |
| Population excluded from age estimation | 0.0242163% | 0.00706918% |
| Covered any TBI | 16.907234% | 20.092000% |
| Covered recent | 2.270459% | 2.313659% |
| Covered remote-only | 14.636775% | 17.778341% |
| Covered prior remote plus recent (subset of recent) | 0.349814% | 0.436243% |
| Full-population point prevalence | null | null |
| Derived minus published both-sex population | 0 | -4 |

The subset is not added to total TBI. Both-sex published rounding differences
remain visible. Every female/male age cell and all five display bands
reconcile; derived-both bands are a separate view, never additional people.
The unresolved tail stays in population totals and in conservative missing-age
bounds. These bounds are not confidence intervals.

Central profile scientific SHA-256:

- Australia: `f58d14a5460679da38f10a53c048bb2fe46f17b0b97afc7cc016a3d68ddd24a0`.
- Viet Nam: `33c00a23d8a62c41cf7385e71438d25fc1ba4ecd9887f6d7b56ccfb80426fd46`.

National shape reconstruction ignores infection acquired elsewhere. Migration
could materially change these profiles, especially Australia's. These figures
must not be interpreted as empirical national burdens.

## WHO trajectory and backcast propagation

Each WHO trajectory is smoothed separately and normalized to its own 2024
level. The alternative historical policy has reviewed diagnostic slope
-0.01/year, annual calendar-grid levels, bounds [0,100], and a 1900 start.
It is an assumed scenario, not an estimated historical trend.

| WHO trajectory / backcast | Australia covered TBI | Viet Nam covered TBI |
| --- | ---: | ---: |
| Lower / constant earliest | 16.808249% | 14.253900% |
| Central / constant earliest | 16.907234% | 20.092000% |
| Upper / constant earliest | 16.888127% | 24.577117% |
| Lower / bounded backward trend | 18.256452% | 15.087474% |
| Central / bounded backward trend | 18.365815% | 21.467331% |
| Upper / bounded backward trend | 18.347374% | 26.248119% |

These are source-data/structural sensitivity trajectories, **not TBI
confidence intervals**. Normalization explains why Australia's WHO upper
trajectory does not produce the highest profile. Alpha is supplied, not
re-fitted in these country scenarios. Scale/assay/anchor parameter scenarios
are separate invented demonstrations below.

## Raw and smoothed incidence through 2020-2022

Incidence units remain WHO cases per 100,000 person-years; these levels are
not infection hazards. Original lower/upper values remain in every export.
The smoother is five-year centred local-linear OLS of positive log incidence,
with truncated three-point endpoint windows, normalized at 2024.

| Year | AUS raw | AUS smoothed | AUS proxy | VNM raw | VNM smoothed | VNM proxy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2019 | 6.8 | 6.714944 | 0.974022 | 176 | 177.478940 | 0.966274 |
| 2020 | 7.2 | 6.497873 | 0.942536 | 171 | 175.549362 | 0.955768 |
| 2021 | 6.5 | 6.436754 | 0.933670 | 171 | 176.124318 | 0.958898 |
| 2022 | 5.6 | 6.455575 | 0.936400 | 178 | 177.309117 | 0.965349 |
| 2023 | 6.2 | 6.370684 | 0.924086 | 185 | 180.952124 | 0.985183 |
| 2024 | 6.9 | 6.894035 | 1 | 182 | 183.673598 | 1 |

The 2020-2022 log residuals are AUS +0.102606, +0.009778, -0.142178 and VNM
-0.026257, -0.029527, +0.003889. The smoother visibly attenuates the Australian
2022 dip. That is descriptive smoothing, not a correction for ascertainment,
COVID disruption or true exposure. Residuals, slopes and exact window members
remain in `shape_diagnostics.csv` for review.

## Reinfection effects

Removing the repeat-event restriction changes recent/remote classification
without changing any-TBI prevalence under the no-clearance assumptions. In
the country diagnostics, recent prevalence increases by 0.349814 percentage
points (Australia) and 0.436243 points (Viet Nam), exactly the auditable
prior-remote-plus-recent subset. With only one event permitted as a comparison,
recent would be 1.920645% and 1.877416%, respectively; remote-only would be
higher by the same amounts. This comparison is not another production mode.

The following constant-pressure example isolates infection-time effects:
alpha*q=0.02/year, exact age 40, recent window 5 years. The conditional mean
and 5th/50th/95th percentiles refer to time since infection among infected
people. Allowing repeat events shifts the most-recent-event distribution
towards baseline. Exact numeric comparison is recorded in the recovery
export-verification receipt; no progression risk is computed.

| Constant-pressure diagnostic | Repeated Poisson events | One event permitted, comparison only |
| --- | ---: | ---: |
| Recent prevalence | 9.516258% | 4.725634% |
| Remote-only prevalence | 45.550845% | 50.341470% |
| Prior remote plus recent | 4.790624% | 0 |
| Any TBI | 55.067104% | 55.067104% |
| Mean time since most recent/only event, years | 17.361351 | 22.638649 |
| 5th percentile, years | 1.395985 | 2.973647 |
| Median, years | 16.102326 | 23.897674 |
| 95th percentile, years | 37.026353 | 38.604015 |

## Synthetic calibration and assay/anchor uncertainty

Invented log-linear incidence, constant earliest backcast and perfect
synthetic observations at female age 5, male age 20 and female age 65 recover
alpha **0.013/year** to better than 1e-18 with explicitly selected WLS.
The invented population is 190,018, of which 190,000 is age-covered and 18
remains in unestimated 100+ cells. A separate binomial test recovers a known
scale from actual invented counts 200/1,000 at age 20 for both sexes.

The following single-anchor demonstration invents IGRA Se/Sp and positivity
0.20 at exact age 20, with equal sex weights. Anchor bounds 0.15-0.25 are an
invented diagnostic range. They are not a real survey or performance study.

| Scenario | Se | Sp | Anchor | Recovered alpha/year |
| --- | ---: | ---: | ---: | ---: |
| Central | 0.85 | 0.96 | 0.20 | 0.009364620 |
| Assay alternative A | 0.75 | 0.98 | 0.20 | 0.012048291 |
| Assay alternative B | 0.95 | 0.92 | 0.20 | 0.006315936 |
| Anchor lower | 0.85 | 0.96 | 0.15 | 0.006210993 |
| Anchor upper | 0.85 | 0.96 | 0.25 | 0.012770796 |

All five choices propagate through the complete synthetic age-sex profile,
including unresolved population. No pseudo-count is manufactured. Parameter,
source-data, structural-scenario and explicit-RNG validation uncertainty remain
separate in the scenario contract.

## PRK blocked example

All selected WHO incidence observations in 2000-2024 are missing. The
structured unavailable artifact identifies all 25 offending years and rejects
shape generation. No country TBI profile is generated, and no zero,
neighbouring-country incidence or arbitrary positive epsilon is substituted.
PRK's frozen 2024 population remains available: 26,498,835 derived people,
with a +1 difference from published both-sex totals retained from Milestone 2B.
An independently reviewed supplied history would be a separate explicit mode.

## Candidate evidence for later review; no values adopted

- [Marks et al., Ca Mau adult population survey](https://pubmed.ncbi.nlm.nih.gov/29471900/):
  a cluster-sampled provincial adult QuantiFERON Gold In-Tube survey. Review
  original age/sex weights, survey time, threshold, design uncertainty and
  provincial-to-national transportability before any use.
- [First national tuberculin survey in Viet Nam](https://pubmed.ncbi.nlm.nih.gov/23676155/):
  a 2006-2007 national stratified cluster survey candidate. Review original
  target ages, TST thresholds, BCG/NTM effects, weights and historical pressure.
- [Dalmau et al., Queensland testing analysis](https://ojs.cdi.cdc.gov.au/index.php/cdi/article/view/137):
  selected TST/IGRA testing data from 2016-2020. This is a potential source for
  reviewing observation/selection processes, not a representative national
  prevalence anchor. A population survey or explicit selection model would
  still be needed.

These primary reports were located for review only. Neither prevalence values
nor assay-performance defaults from them are entered into a contract or fit.
