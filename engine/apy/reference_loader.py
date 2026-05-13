from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def load_reference_bundle(path: str | Path) -> dict[str, Any]:
    return _load_reference_json(path, "MATLAB reference bundle")


def load_reference_summary(path: str | Path) -> pd.DataFrame:
    reference_path = _require_reference_file(path, "MATLAB reference summary CSV")
    return pd.read_csv(reference_path)


def load_reference_dynamic_comparison(path: str | Path) -> dict[str, Any]:
    return _load_reference_json(path, "MATLAB dynamic-comparison reference")


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
