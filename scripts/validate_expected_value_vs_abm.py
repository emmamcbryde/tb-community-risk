from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from engine.apy.expected_value import run_expected_value
from engine.apy.runner import run_replicates


EVENTS = [
    "recent_ltbi_at_baseline",
    "remote_ltbi_at_baseline",
    "recent_latent_at_screen",
    "remote_latent_at_screen",
    "screened",
    "true_positive_latent",
    "true_positive_recent",
    "true_positive_remote",
    "false_positive",
    "tpt_started_total",
    "tpt_completed_total",
    "infection_effectively_treated_total",
    "infection_effectively_treated_recent",
    "infection_effectively_treated_remote",
    "active_tb_cases_prevented",
    "active_tb_cases_prevented_recent",
    "active_tb_cases_prevented_remote",
    "active_tb_cases",
]


SCENARIOS = {
    "random_simple": {
        "N": 700,
        "nReps": 60,
        "seed": 2401,
        "screeningStrategy": "random",
        "screenCoverage": 0.35,
        "baselineRecentLTBIProportion": 0.35,
        "ltbiStateAssumptions": {
            "baselineRecentLTBIProportion": 0.35,
            "recentToRemoteTransitionRatePerYear": 0.2,
            "status": "configured",
        },
    },
    "apy_prevent_targeted": {
        "N": 700,
        "nReps": 60,
        "seed": 2402,
        "screeningStrategy": "prevent",
        "screenCoverage": 0.30,
        "baselineRecentLTBIProportion": 0.35,
        "ltbiStateAssumptions": {
            "baselineRecentLTBIProportion": 0.35,
            "recentToRemoteTransitionRatePerYear": 0.2,
            "status": "configured",
        },
    },
}


def run_validation(configs: dict[str, dict[str, Any]] | None = None) -> pd.DataFrame:
    rows = []
    for scenario_name, config in (configs or SCENARIOS).items():
        ev = _wide(run_expected_value(config, quadrature_points=48)["eventLedger"]["replicateTotals"])
        abm = _wide(run_replicates(config, keep_example_cohort=False)["eventLedger"]["replicateTotals"])
        for arm in ["comparator", "intervention"]:
            ev_row = ev[ev["arm"] == arm].iloc[0]
            abm_arm = abm[abm["arm"] == arm]
            for event in EVENTS:
                if arm == "comparator" and event not in {"active_tb_cases"}:
                    continue
                values = abm_arm[event].astype(float)
                mean = float(values.mean())
                mcse = float(values.std(ddof=1) / (len(values) ** 0.5)) if len(values) > 1 else 0.0
                floor = 0.25
                tolerance = max(3.0 * mcse, floor)
                delta = float(ev_row[event]) - mean
                rows.append(
                    {
                        "scenario": scenario_name,
                        "arm": arm,
                        "eventName": event,
                        "expectedValue": float(ev_row[event]),
                        "abmMean": mean,
                        "abmMCSE": mcse,
                        "delta": delta,
                        "tolerance": tolerance,
                        "withinTolerance": abs(delta) <= tolerance,
                    }
                )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true", help="Emit JSON records instead of a table.")
    args = parser.parse_args()
    results = run_validation()
    if args.json:
        print(json.dumps(results.to_dict(orient="records"), indent=2))
    else:
        print(results.to_string(index=False))
    return 0 if bool(results["withinTolerance"].all()) else 1


def _wide(totals: pd.DataFrame) -> pd.DataFrame:
    frame = totals.copy()
    id_cols = ["replicateId", "arm"]
    return frame.pivot_table(
        index=id_cols,
        columns="eventName",
        values="value",
        aggfunc="first",
    ).reset_index()


if __name__ == "__main__":
    raise SystemExit(main())
