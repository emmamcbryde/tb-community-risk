# Milestone 1 Scenario and Economics Foundation

APY is now represented as the `apy_demonstration` population preset in the same
versioned scenario contract used for future populations. The generic preset
table is `data/population_presets.json`; no APY-only branch is required to load
the demonstration scenario.

Scenario records include identity/version, population size, age-profile source,
LTBI prevalence assumptions, risk-factor assumptions, targeting criteria,
eligible and screened counts or proportions, screening window, follow-up
horizon, comparator, intervention, sources, notes, and the direct-effects scope
statement.

The default comparator is current practice with no additional systematic LTBI
screening. It includes ordinary background clinical care where already
represented. The default intervention is targeted LTBI screening and preventive
treatment. IGRA is the demonstration test, TST remains selectable, and
sensitivity and specificity are stored separately from coverage and treatment
parameters. The demonstration regimen is 3HP, with existing alternative regimen
labels retained.

Economics use an Australian health-system perspective, AUD target currency, and
2025-26 target price year. Existing local cost inputs are preserved as original
source values. Their source price years are unresolved unless explicitly
verified, so they are not silently treated as 2025-26 values.

Cost-year conversion is separate from future discounting:

```text
converted_cost = original_cost * target_year_index / source_year_index
```

Inflation index values are version-controlled in `data/inflation_indices.csv`
and are never retrieved live during a run. Missing source years or missing index
values produce unresolved warnings instead of fabricated costs. Downstream
economics calculations only use converted target-year costs when conversion
status is valid. Model-year costs remain in constant target-year prices; future
inflation is not applied.

The editable discounting structure supports 3% and 0% annual rates. Discounting
is reserved for future costs and health outcomes, not historical source-price
normalisation.

The primary health outcome is DALYs averted. Natural outputs remain visible
alongside DALY configuration, including TB cases prevented, infections
effectively treated, false-positive treatments, and adverse events where the
model already supports them. Missing DALY or mortality inputs should remain
explicitly unresolved rather than being invented.

The willingness-to-pay benchmark is stored as an editable cost per DALY averted.
The demonstration concept is 1 x Australian GDP per capita per DALY averted, but
no live GDP lookup is performed and no unverified value is supplied. Treat it as
an illustrative benchmark, not an official Australian funding threshold.

Current scope statement:

```text
The current static and individual-based analyses estimate direct benefits,
harms and costs for screened and treated individuals. Transmission-mediated
benefits are not yet included.
```
