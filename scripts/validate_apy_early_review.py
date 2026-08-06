from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.apy.config import build_default_config
from engine.apy.early_review import run_early_screening_review
from engine.apy.ltbi_state import apply_ltbi_state_edit


def main() -> None:
    cfg = apply_ltbi_state_edit(
        build_default_config(),
        baseline_recent_proportion=0.35,
        transition_rate_per_year=0.2,
        source="Synthetic validation fixture, not APY reference evidence",
        status="configured_reviewed",
        notes="Used only for software validation.",
    )
    cfg.update({"N": 30, "screeningStrategy": "random"})
    common = {
        "reviewId": "synthetic_validation",
        "screenedToDate": 6,
        "plannedTotalScreened": 12,
        "eligiblePopulation": 30,
        "prior": {
            "type": "discrete_grid",
            "weights": [0.5, 0.5],
            "source": "Synthetic validation prior, not APY reference evidence",
        },
        "prevalenceGrid": [0.0, 0.1],
    }
    low = run_early_screening_review(cfg, None, {**common, "testPositiveToDate": 0})
    high = run_early_screening_review(cfg, None, {**common, "testPositiveToDate": 5})
    assert low["validation"]["isValid"], low["validation"]
    assert high["validation"]["isValid"], high["validation"]
    assert abs(sum(row["posteriorWeight"] for row in low["priorPosteriorTable"]) - 1.0) < 1e-9
    assert low["posterior"]["summary"]["mean"] < high["posterior"]["summary"]["mean"]
    print("Synthetic early-review validation passed.")
    print("Low-yield posterior mean:", low["posterior"]["summary"]["mean"])
    print("High-yield posterior mean:", high["posterior"]["summary"]["mean"])


if __name__ == "__main__":
    main()
