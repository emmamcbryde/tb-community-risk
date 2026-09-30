"""Run the Starsim TB natural-history feasibility prototype (Milestone 1).

Demonstration values only: not calibrated, not a clinical or country model, not
for policy use. Requires the separate prototype environment
(``requirements-starsim-prototype.txt``).

Examples:
    python scripts/run_starsim_tb_prototype.py
    python scripts/run_starsim_tb_prototype.py --seed 7 --json
    python scripts/run_starsim_tb_prototype.py --out-dir outputs/starsim_tb_m1 --csv
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.starsim_tb.parameters import DEMONSTRATION_VALUES_NOTICE, demonstration_parameters  # noqa: E402
from engine.starsim_tb.prototype import results_table, run_prototype, summary  # noqa: E402

DEFAULT_OUT_DIR = ROOT / "outputs" / "starsim_tb_m1"
RESULTS_FORMAT = "starsim_tb_m1_prototype_results/1"  # Not a production result contract


def _peak_memory_mb() -> float | None:
    try:
        import psutil
    except ImportError:
        return None
    info = psutil.Process().memory_info()
    peak = getattr(info, "peak_wset", None)  # Windows peak working set
    if peak is None:
        try:
            import resource

            usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            peak = usage * (1 if sys.platform == "darwin" else 1024)
        except ImportError:
            peak = info.rss
    return round(peak / 2**20, 1)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, default=None, help="random seed (default: demonstration seed)")
    parser.add_argument("--n-agents", type=int, default=None)
    parser.add_argument("--stop-year", type=float, default=None)
    parser.add_argument("--json", action="store_true", help="write machine-readable JSON results")
    parser.add_argument("--csv", action="store_true", help="write the per-timestep table as CSV")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args(argv)

    changes = {}
    if args.seed is not None:
        changes["rand_seed"] = args.seed
    if args.n_agents is not None:
        changes["n_agents"] = args.n_agents
    if args.stop_year is not None:
        changes["stop_year"] = args.stop_year
    pars = demonstration_parameters(**changes)

    run = run_prototype(pars)
    metadata = dict(run.metadata, peak_process_memory_mb=_peak_memory_mb())
    summ = summary(run)
    table = results_table(run)

    print(f"NOTICE: {DEMONSTRATION_VALUES_NOTICE}")
    print(f"Starsim {metadata['starsim_version']} | Python {metadata['python_version']} | {metadata['platform']}")
    print(f"Commit {metadata['git_commit']} | parameter hash {metadata['parameter_hash_sha256'][:16]} | seed {pars.rand_seed}")
    print(f"{pars.n_agents} agents, {pars.start_year:g}-{pars.stop_year:g}, dt = {pars.timestep_years:.6g} y ({pars.n_timesteps} steps)")
    print(f"Runtime {run.runtime_seconds:.2f} s | peak process memory {metadata['peak_process_memory_mb']} MB")
    print(json.dumps(summ, indent=2))

    if args.json or args.csv:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        stem = f"starsim_tb_m1_seed{pars.rand_seed}"
        if args.json:
            path = args.out_dir / f"{stem}.json"
            payload = {"format": RESULTS_FORMAT, "metadata": metadata, "summary": summ, "timeseries": table}
            path.write_text(json.dumps(payload, indent=1, allow_nan=False), encoding="utf-8")
            print(f"Wrote {path}")
        if args.csv:
            path = args.out_dir / f"{stem}_timeseries.csv"
            with path.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(table.keys())
                writer.writerows(zip(*table.values()))
            print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
