"""Closed-form SYNTHETIC backcast comparison; no country data or engine imports."""
from __future__ import annotations

import argparse
import csv
import io
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AGES = (0, 1, 4, 10, 24, 40, 65, 85)


def exponential_integral(start, end, slope):
    """Integral of 0.01 exp[-slope(t-2000)] with explicit finite endpoints."""
    if end <= start:
        return 0.0
    if slope == 0:
        return .01 * (end - start)
    return .01 * (math.exp(-slope * (start - 2000)) - math.exp(-slope * (end - 2000))) / slope


def cumulative(age, policy, *, multiplier=1.0):
    birth = 2024 - age
    observed = exponential_integral(max(2000, birth), 2024, .02)
    if birth >= 2000:
        return observed
    if policy == "constant_earliest":
        earlier = .01 * (2000 - birth)
    elif policy in {"bounded_early_trend", "anchor_conditional_slope"}:
        slope = .02 if policy == "bounded_early_trend" else .01
        if .01 * math.exp(slope * (2000 - birth)) > .1:
            raise ValueError("Synthetic bound exceeded; reject instead of capping.")
        earlier = exponential_integral(birth, 2000, slope)
    elif policy == "supplied_history":
        earlier = sum(rate * max(0, min(2000, hi) - max(birth, lo))
                      for lo, hi, rate in ((1900, 1960, .005), (1960, 1980, .02), (1980, 2000, .012)))
        if birth < 1900:
            raise ValueError("Missing supplied history.")
    elif policy == "uncertain_earlier_level":
        earlier = .01 * multiplier * (2000 - birth)
    else:
        raise ValueError("Unknown backcast.")
    return earlier + observed


def examples():
    rows = []
    for policy in ("constant_earliest", "bounded_early_trend", "supplied_history",
                   "anchor_conditional_slope", "uncertain_earlier_level"):
        for age in AGES:
            hazards = [cumulative(age, policy, multiplier=m) for m in (.5, 1, 2)] if policy == "uncertain_earlier_level" else [cumulative(age, policy)]
            probabilities = [-math.expm1(-h) for h in hazards]
            weights = (.25, .5, .25) if len(hazards) == 3 else (1,)
            total = sum(w * p for w, p in zip(weights, probabilities))
            hr = exponential_integral(2024 - min(age, 2), 2024, .02)
            recent = sum(w * math.exp(-(h - hr)) * -math.expm1(-hr) for w, h in zip(weights, hazards))
            rows.append({"policy": policy, "age": age, "total_ever_infected": total,
                         "recent_first_infected": recent, "remote_only": total - recent,
                         "scenario_min": min(probabilities), "scenario_max": max(probabilities)})
    return rows


def table():
    rows = examples()
    lines = ["| Backcast | Age 4 | Age 24 | Age 40 | Age 65 | Age 85 |",
             "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for policy in dict.fromkeys(r["policy"] for r in rows):
        values = [r["total_ever_infected"] for r in rows if r["policy"] == policy and r["age"] in (4, 24, 40, 65, 85)]
        lines.append("| " + policy + " | " + " | ".join(f"{100*p:.3f}%" for p in values) + " |")
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if args.output:
        output = args.output.resolve()
        if not output.is_relative_to((ROOT / "abm/output").resolve()):
            raise ValueError("Preview output must stay under abm/output.")
        output.parent.mkdir(parents=True, exist_ok=True)
        stream = io.StringIO(newline="")
        rows = examples()
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        output.write_text(stream.getvalue(), encoding="utf-8")
    print(table())


if __name__ == "__main__":
    main()
