from __future__ import annotations

import argparse
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from engine.apy.distributional_validation import (  # noqa: E402
    run_reference_suite_distributional_validation,
)


DEFAULT_REFERENCE_ROOT = Path("validation") / "matlab_reference"
DEFAULT_SUITE_FILE = DEFAULT_REFERENCE_ROOT / "scenario_suite_v1.json"
DEFAULT_OUTPUT_DIR = (
    Path("validation") / "output" / "apy_python_matlab_distributional_validation"
)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    scenario_ids = _split_scenario_ids(args.scenario_id)
    config_overrides = {"nReps": int(args.quick)} if args.quick else None

    result = run_reference_suite_distributional_validation(
        reference_root=args.reference_root,
        suite_file=args.suite_file,
        scenario_ids=scenario_ids,
        config_overrides=config_overrides,
    )
    write_outputs(result, Path(args.output_dir))
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Python APY vs MATLAB reference distributional validation."
    )
    parser.add_argument(
        "--reference-root",
        default=str(DEFAULT_REFERENCE_ROOT),
        help="Directory containing MATLAB reference scenario fixtures.",
    )
    parser.add_argument(
        "--suite-file",
        default=str(DEFAULT_SUITE_FILE),
        help="Scenario suite JSON file.",
    )
    parser.add_argument(
        "--scenario-id",
        action="append",
        default=[],
        help="Scenario ID to run. Can be repeated or comma-separated.",
    )
    parser.add_argument(
        "--quick",
        nargs="?",
        const=20,
        type=int,
        help="Override nReps for faster local checks. Defaults to 20 if no value is supplied.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory to write validation outputs.",
    )
    return parser.parse_args(argv)


def write_outputs(result: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    scenario_rows = result["scenarioRows"]
    metric_rows = result["metricRows"]
    scenario_rows.to_csv(output_dir / "scenario_summary.csv", index=False)
    metric_rows.to_csv(output_dir / "metric_validation_rows.csv", index=False)
    (output_dir / "validation_report.md").write_text(
        build_markdown_report(scenario_rows, metric_rows),
        encoding="utf-8",
    )


def build_markdown_report(scenario_rows, metric_rows) -> str:
    lines = [
        "# APY Python vs MATLAB Distributional Validation",
        "",
        "This report is diagnostic. Broad tolerances are used because MATLAB and NumPy random streams differ.",
        "",
        "## Scenario Summary",
        "",
    ]
    if scenario_rows.empty:
        lines.append("No scenarios were evaluated.")
    else:
        lines.extend(_markdown_table(scenario_rows))

    lines.extend(["", "## Metric Rows", ""])
    if metric_rows.empty:
        lines.append("No metric rows were produced.")
    else:
        display_cols = [
            "scenario_id",
            "Metric",
            "PythonMedian",
            "MatlabMedian",
            "AbsoluteDifference",
            "RelativeDifference",
            "Tolerance",
            "Pass",
        ]
        lines.extend(_markdown_table(metric_rows[display_cols]))
    lines.append("")
    return "\n".join(lines)


def _markdown_table(df) -> list[str]:
    columns = list(df.columns)
    rows = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for record in df.to_dict(orient="records"):
        rows.append(
            "| "
            + " | ".join(_format_cell(record.get(column)) for column in columns)
            + " |"
        )
    return rows


def _format_cell(value) -> str:
    if value != value:
        return ""
    return str(value).replace("|", "\\|")


def _split_scenario_ids(values: list[str]) -> list[str] | None:
    scenario_ids = []
    for value in values:
        scenario_ids.extend(
            part.strip() for part in str(value).split(",") if part.strip()
        )
    return scenario_ids or None


if __name__ == "__main__":
    raise SystemExit(main())
