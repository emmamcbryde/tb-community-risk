from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def load_reference_bundle(path: str | Path) -> dict[str, Any]:
    return _load_reference_json(path, "MATLAB reference bundle")


def load_reference_manifest(path: str | Path) -> dict[str, Any]:
    return _load_reference_json(path, "MATLAB reference manifest")


def load_reference_scenario_config(path: str | Path) -> dict[str, Any]:
    return _load_reference_json(path, "MATLAB reference scenario config")


def load_reference_summary(path: str | Path) -> pd.DataFrame:
    reference_path = _require_reference_file(path, "MATLAB reference summary CSV")
    return pd.read_csv(reference_path)


def load_reference_dynamic_comparison(path: str | Path) -> dict[str, Any]:
    return _load_reference_json(path, "MATLAB dynamic-comparison reference")


def load_reference_scenario_suite(path: str | Path) -> list[dict[str, Any]]:
    reference_path = _require_reference_file(path, "MATLAB reference scenario suite")
    with reference_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list):
        raise ValueError(
            f"MATLAB reference scenario suite must contain a JSON list: {reference_path}"
        )
    for index, item in enumerate(data):
        if not isinstance(item, dict) or "scenario_id" not in item:
            raise ValueError(
                "Each scenario suite entry must be an object with scenario_id "
                f"(entry {index})"
            )
    return data


def list_reference_scenario_dirs(reference_root: str | Path) -> list[Path]:
    root = Path(reference_root)
    if not root.is_dir():
        raise FileNotFoundError(f"MATLAB reference root not found: {root}")
    return sorted(
        path for path in root.iterdir()
        if path.is_dir() and (path / "manifest.json").is_file()
    )


def load_reference_dir(path: str | Path) -> dict[str, Any]:
    reference_dir = Path(path)
    if not reference_dir.is_dir():
        raise FileNotFoundError(f"MATLAB reference directory not found: {reference_dir}")

    return {
        "manifest": load_reference_manifest(reference_dir / "manifest.json"),
        "scenario_config": load_reference_scenario_config(
            reference_dir / "scenario_config.json"
        ),
        "dynamic_comparison": load_reference_dynamic_comparison(
            reference_dir / "matlab_dynamic_comparison.json"
        ),
        "summary": load_reference_summary(reference_dir / "matlab_summary.csv"),
    }


def _load_reference_json(path: str | Path, label: str) -> dict[str, Any]:
    reference_path = _require_reference_file(path, label)
    with reference_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"{label} must contain a JSON object: {reference_path}")
    return data


def _require_reference_file(path: str | Path, label: str) -> Path:
    reference_path = Path(path)
    if not reference_path.is_file():
        raise FileNotFoundError(f"{label} not found: {reference_path}")
    return reference_path
