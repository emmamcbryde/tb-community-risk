# MATLAB Reference Outputs

MATLAB APY v9 remains the reference implementation while the Python-native APY ABM port is developed.

Reference outputs in this folder are intended for Python validation. They should be compact, portable, and committed only when they are useful as stable test fixtures.

The Python APY port must match MATLAB reference outputs within documented tolerances.
Exact individual-level equality is not expected because MATLAB and NumPy random streams differ.

Validation should focus on:

- summary statistics;
- medians and uncertainty intervals;
- stable output-contract fields;
- `technical.dynamicComparison` values where available;
- compact deterministic parameter snapshots where available;
- scenario configuration compatibility.

Generated reference outputs should not be treated as general model outputs.
Keep ordinary generated APY files in `abm/output/`, and only commit reference files here when they are intentionally small and useful for tests.
