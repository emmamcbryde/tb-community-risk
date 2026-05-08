from __future__ import annotations

import json
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
