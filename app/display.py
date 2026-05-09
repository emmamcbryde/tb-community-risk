from __future__ import annotations

import json
import re
from datetime import datetime
from typing import Any

import pandas as pd


def arrow_safe_dataframe(data: Any) -> pd.DataFrame:
    """Return a display-only DataFrame with Arrow-safe object columns."""
    frame = pd.DataFrame(data).copy()
    for column in frame.columns:
        if frame[column].dtype == "object":
            frame[column] = frame[column].map(display_string)
    return frame


def display_string(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True)
    return str(value)


def safe_download_stem(label: object, suffix: str) -> str:
    label_text = str(label or "scenario").strip() or "scenario"
    safe_label = re.sub(r"[^A-Za-z0-9._-]+", "_", label_text).strip("._-")
    safe_label = safe_label or "scenario"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{safe_label}_{timestamp}_{suffix}"


def economics_summary_csv(economics_results: dict[str, Any]) -> bytes:
    rows = economics_results.get("summaryRows") or []
    return pd.DataFrame(rows).to_csv(index=False).encode("utf-8")


def economics_assumptions_json(config: dict[str, Any]) -> str:
    return json.dumps(config, indent=2, sort_keys=True)
