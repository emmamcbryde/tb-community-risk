"""Maintainer helper: export cross-check columns from a local Global TB Report checkout.

Reads inc_mort/analysis/est.rda (read-only) and writes a CSV used only to
cross-check the public WHO estimates during import. Requires the optional
``pyreadr`` package; neither the application nor its tests need it. The export
is not committed to this repository.

    python scripts/export_gtbreport_crosscheck.py ../gtbreport2025 gtbreport2025_est_export.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

COLUMNS = ["iso3", "year", "inc", "inc.lo", "inc.hi", "inc.num", "inc.lo.num", "inc.hi.num", "pop"]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("gtbreport_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)
    try:
        import pyreadr
    except ImportError:
        print("This helper needs the optional 'pyreadr' package.")
        return 1
    est = pyreadr.read_r(str(args.gtbreport_root / "inc_mort" / "analysis" / "est.rda"))["est"]
    missing = [column for column in COLUMNS if column not in est.columns]
    if missing:
        print(f"est.rda schema changed; missing {missing}.")
        return 1
    est[COLUMNS].sort_values(["iso3", "year"]).to_csv(args.output, index=False, float_format="%.17g", lineterminator="\n")
    print(f"Wrote {len(est)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
