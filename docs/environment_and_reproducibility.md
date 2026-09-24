# Environment and reproducibility

Status: **implemented and validated** (Milestone 2).

## Supported configuration

| Item | Supported |
| --- | --- |
| Python | 3.12.x (tested 3.12.12) |
| Streamlit | 1.39.0 (pinned in `requirements.txt`) to 1.54.0; both tested |
| Other packages | exactly as pinned in `requirements.txt` (pandas 2.2.2, numpy 1.26.4, altair 5.3.0, matplotlib 3.9.2; scipy and openpyxl unpinned) |
| MATLAB | not required |
| Internet | not required by the application or tests |

`requirements.txt` and `runtime.txt` are part of the frozen release and are not
edited on this branch. `runtime.txt` contains the literal text `runtime.txt`
rather than a Python version; Streamlit Community Cloud takes its Python version
from the app settings, so this is recorded rather than changed.

Check any environment with:

```bash
python scripts/check_environment.py
```

It exits non-zero for an unsupported Python or Streamlit version and warns about
version drift and about x86-64 Python emulated on ARM64.

## Clean environment (Windows, uv)

```bash
uv venv --python 3.12 .venv-tb
uv pip install --python .venv-tb/Scripts/python.exe -r requirements.txt pytest
.venv-tb/Scripts/python.exe scripts/check_environment.py
```

If uv cannot download a Python build because of a corporate TLS proxy, add
`--native-tls` or point `--python` at an existing 3.12 interpreter. On macOS or
Linux use `.venv-tb/bin/python`.

Conda alternative:

```bash
conda create -n tbmodel python=3.12
conda activate tbmodel
pip install -r requirements.txt pytest
```

## Tests

```bash
python -m pytest -q                                    # full suite (about 2.5 h on the reference machine)
python -m pytest -q tests/test_general_*.py            # general application (a few minutes)
python -m pytest -q tests/test_frozen_release_integrity.py tests/test_calibration_memoisation.py
```

The general-application tests need neither MATLAB, internet access nor a
`gtbreport2025` checkout.

## Investigation of the historical 90-second page timeouts

Two existing rendered-page tests timed out after 90 s in Milestone 1:

* `test_apy_infection_history.py::...test_supported_deterministic_run_survives_results_and_economics_navigation`
* `test_sa_health_reference_package.py::...test_rendered_health_economics_widgets_recalculate_without_changing_health`

Findings:

1. **The Streamlit version is not the cause.** Both tests were run in clean Python
   3.12.12 environments with Streamlit 1.39.0 (exact `requirements.txt`) and 1.54.0.
   Test 1 timed out in both (101 s and 99 s). Test 2 passed in both (366 s and 371 s
   wall time, including class set-up) but failed in the `tbmodel` conda environment.
   `pip check` reported no broken requirements in any environment.
2. **Where the time goes.** Profiling the first deterministic preview in a fresh process
   (N = 50) showed 164.9 s under the profiler, of which 141.7 s was
   `resolve_calibration_for_config`. `calibrate_from_config` was called **twice with
   identical inputs**: once directly and once inside
   `build_reference_calibration_artifact`. A second preview in the same process took
   14.8 s. None of the following contributed: workbook construction (lazy), chart
   rendering, frozen-artifact loading, session-state recursion or stochastic
   execution. No rendered test starts the 1,000-run analysis.
3. **Why this machine is slow.** The reference machine has a Snapdragon X Elite (ARM64)
   CPU. The pinned wheels (pandas 2.2.2, numpy 1.26.4) exist only for x86-64 on
   Windows, so Python runs under x86-64 emulation. The pure-Python calibration loop
   takes about 38 s per call there. A native ARM64 environment could not be built,
   because pandas 2.x publishes no Windows ARM64 wheels.
4. **Correction.** `engine/apy/calibration_policy.py` now memoises the duplicated call,
   keyed by the hash of the normalised configuration and returning deep copies.
   `tests/test_calibration_memoisation.py` shows that calibration runs once and that
   the results match an uncached call exactly. No formula, timeout or assertion was
   changed.
5. **Result after the correction.** Test 1 passes in 64 s (Streamlit 1.39) and 62 s
   (Streamlit 1.54). Test 2 passes in the `tbmodel` environment in 290 s, of which
   139 s is the timed interaction.

On a native x86-64 machine the same code is expected to be several times faster.
The 90 s limits remain tight on emulated ARM64 hardware: `check_environment.py`
warns about this, and heavy background load on the machine can still cause
timeouts.

## Test-isolation finding (full suite, clean Streamlit 1.39 environment)

`tests/test_dynamic_abm_compare_page_helpers.py`, part of the frozen release,
replaces `sys.modules["streamlit"]` with a mock at module level and never restores
it. Any later test in the same process that imports Streamlit at run time,
including `AppTest._run`, which reads `st.secrets`, then gets the mock. In the first
full run, all 16 general-interface tests that sort after it failed with
`module 'streamlit' has no attribute 'secrets'`. They pass when run alone.

`tests/test_general_app_interface.py` now restores the real Streamlit module
before each render; no assertions were changed. The existing test file was not
edited, because it belongs to the frozen release. The recommended future fix is
to wrap its module replacement in `patch.dict(sys.modules, ...)`.

## Windows cleanup warnings

* Streamlit may log `missing ScriptRunContext` and `use_container_width` deprecation
  messages during AppTest runs (Streamlit 1.54). They are harmless.
* The existing Word-report test leaves `reports/test_tmp/` and `outputs/test_tmp/`.
  They are untracked and can be deleted.
* Temporary directories created by tests can occasionally fail to delete on Windows
  while a file handle is open. Retry or delete them manually.

## Expected page timing on the reference machine (emulated x86-64)

| Operation | Time |
| --- | --- |
| General Set up population page (production snapshot, first render) | about 2-4 s |
| Snapshot load with full re-validation (first time in a process) | 0.9 s, cached afterwards |
| Trend fit with 1,000 propagated draws | < 0.1 s (smooth), about 0.1-1 s (log-linear, first call includes import) |
| Deterministic preview, first in a process (calibration included) | about 55-70 s at 10,000 people |
| Deterministic preview, later in the same process | about 15-25 s |
| Stochastic analysis, 1,000 simulated populations of 10,000 | see `docs/performance_benchmark.md` |
