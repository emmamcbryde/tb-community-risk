# WHO incidence import validation

Result: PASSED

## Fatal (0)
None.

## Warning (150)
- `crosscheck_unpublished` x150: PRK 2000 incidence_per_100k: present in the report-round analysis output but not published in the public WHO dataset; the public (blank) value is kept.

## Expected missingness (1)
- `no_estimates` x1: PRK: population published but no incidence estimates.

## Incomplete series (8)
- `short_series` x8: ANT: series covers 2000-2009 rather than 2000-2024.

## Info (2)
- `crosscheck_passed` x1: All 37279 values published in both sources agree within published rounding precision.
- `unused_columns` x1: 34 published column(s) are not part of this snapshot contract.
