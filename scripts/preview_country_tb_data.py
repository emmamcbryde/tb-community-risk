"""Offline developer data preview; absent from all deployed app navigation."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.country_tb.data_profile import baseline_data_profile, load_candidates, write_preview
from engine.country_tb.incidence import WHO_2025, load_who_snapshot


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--countries", nargs="+", default=["AUS", "VNM", "PRK"])
    parser.add_argument("--snapshot", default=WHO_2025)
    parser.add_argument("--output", type=Path, default=ROOT / "abm/output/country_tb_m2a_preview")
    args = parser.parse_args(argv)
    snapshot = load_who_snapshot(args.snapshot)
    candidates = load_candidates()
    for country in args.countries:
        profile = baseline_data_profile(country, snapshot, candidates)
        print(write_preview(profile, args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
