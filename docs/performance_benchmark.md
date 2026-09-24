# Performance benchmark

Status: **measured**. Reference machine:
* Snapdragon X Elite X1E78100, 12 cores, Windows 11;
* Python 3.12.12 x86-64 running under emulation;
* Streamlit 1.54, pandas 2.2.2, numpy 1.26.4.

Native x86-64 machines are expected to be several times faster for the pure-Python
parts (calibration).

## Light operations (median of 5)

| Operation | Time |
| --- | --- |
| Snapshot load with full re-validation (first time in a process) | 553 ms |
| Snapshot load (cached) | 0.2 ms |
| Country selector options (217 areas) | 45 ms |
| Log-linear trend with 1,000 propagated draws | 5 ms |
| Log-linear trend without propagation | 3 ms |
| Smooth trend with 1,000 propagated draws | 13 ms |
| Incidence chart specification | 112 ms |
| Profile hash / JSON round trip | 0.5 ms / 2.9 ms |
| Engine configuration build | 17 ms |

Propagation is vectorised: both estimators are linear in log incidence, so 1,000
refits are a single matrix product. The results are identical to refitting each
draw in a loop.

## Analyses (population 10,000, demonstration profile)

| Analysis | Time |
| --- | --- |
| Deterministic expected-value preview, first in a process (calibration included) | about 55-70 s |
| Deterministic preview, later in the same process | about 15-25 s |
| Stochastic, 20 simulated populations | calibration 37 s + 4.6 s per population when profiled (overhead included) |
| **Stochastic default, 1,000 simulated populations** | **534 s (8.9 min) including calibration; 0.53 s per population** |
| Health economics for 1,000 simulated populations | 175 s (2.9 min); a cost-only recalculation takes the same time (180 s) and does not rerun epidemiology |

The 1,000-population run was measured while other work was running on the machine,
so these figures are conservative. The interface estimate uses 0.55 s per population
of 10,000 plus 40 s for calibration, and about 0.18 s per population for economics.
Before this milestone the estimate assumed 0.9 s per population (about 15 minutes),
which was too pessimistic.

## Where the time goes (20-simulation run under the profiler; relative shares)

| Component | Time under the profiler |
| --- | --- |
| Reference calibration (once per process per epidemiological configuration) | about 37 s unprofiled |
| Event-record construction (`_append_agent_based_ledger_rows`, binning) | about 0.4 s per simulated population |
| Individual simulation (`simulate_one_cohort`) | about 0.2 s per simulated population |
| Result serialisation (`to_json_like`) | about 0.4 s per simulated population |
| Health economics from a completed run | about 0.2 s per simulated population |

## Caching and reuse

* **Calibration** is cached in-process, keyed by its epidemiological inputs. The
  duplicated calibration call in the full reference policy is now memoised
  (outputs identical; `tests/test_calibration_memoisation.py`). Changing only the
  test, regimen, coverage or targeting does not recalibrate.
* **Identical configurations** reuse completed results (last 3 per session),
  matched by epidemiological configuration hash.
* **Cost-only changes** recalculate economics from the completed run; epidemiology
  is never rerun.
* **Downloads** (provenance package and workbook) are built only on request.

## Not optimised in this milestone

These would need equivalence tests before changing:
* vectorising the pure-Python calibration loops;
* lighter result serialisation;
* running simulated populations in parallel.

The 1,000-simulation default was kept as requested. The interface shows an
estimated run time and asks for confirmation before long runs.
