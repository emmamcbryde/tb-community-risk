from __future__ import annotations

from collections.abc import Mapping, Sequence
import json
import math
from numbers import Number
from typing import Any


def to_json_like(value: Any) -> Any:
    """Convert MATLAB/Python values into Streamlit-safe JSON-like objects."""
    if value is None:
        return None

    if _is_pandas_dataframe(value):
        return _pandas_dataframe_to_rows(value)

    if _is_pandas_series(value):
        return to_json_like(value.to_dict())

    if _is_numpy_array(value):
        return to_json_like(value.tolist())

    if _is_numpy_value(value):
        return to_json_like(value.item())

    if _is_matlab_table(value):
        return _matlab_table_to_rows(value)

    if _is_matlab_array(value):
        return _matlab_array_to_json_like(value)

    if isinstance(value, Mapping):
        return {str(key): to_json_like(item) for key, item in value.items()}

    if hasattr(value, "_fieldnames"):
        return {
            str(field): to_json_like(getattr(value, field))
            for field in value._fieldnames
        }

    if isinstance(value, bool):
        return value

    if isinstance(value, str):
        return value

    if isinstance(value, Number):
        return _number_to_json_like(value)

    if hasattr(value, "isoformat"):
        return value.isoformat()

    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")

    if hasattr(value, "tolist"):
        return to_json_like(value.tolist())

    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [to_json_like(item) for item in value]

    return value


def ensure_json_safe(value: Any) -> Any:
    """Normalize a value and verify that plain json.dumps can serialize it."""
    normalized = to_json_like(normalize_for_matlab(value))
    json.dumps(normalized)
    return normalized


def normalize_for_matlab(value: Any) -> Any:
    """Remove Python-only containers before sending payloads to MATLAB."""
    if isinstance(value, Mapping):
        out: dict[str, Any] = {}
        for key, item in value.items():
            key_text = str(key)
            if key_text == "ageDistributionTable":
                if _is_empty_table_like(item):
                    out[key_text] = []
                else:
                    out["ageDistributionTableRows"] = _table_like_to_rows(item)
                    out[key_text] = []
                continue
            out[key_text] = normalize_for_matlab(item)
        return out

    if _is_pandas_dataframe(value):
        return _pandas_dataframe_to_rows(value)

    if _is_pandas_series(value):
        return to_json_like(value.to_dict())

    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [normalize_for_matlab(item) for item in value]

    return value


def to_matlab_arg(engine: Any, value: Any) -> Any:
    """Convert JSON-like Python values to MATLAB values via jsondecode."""
    value = ensure_json_safe(value)
    if isinstance(value, Mapping) or _is_json_sequence(value):
        return engine.jsondecode(json.dumps(value), nargout=1)
    return value


def _number_to_json_like(value: Number) -> int | float:
    if isinstance(value, bool):
        return value
    number = float(value)
    if math.isnan(number) or math.isinf(number):
        return number
    if number.is_integer():
        return int(number)
    return number


def _is_numpy_value(value: Any) -> bool:
    return value.__class__.__module__.startswith("numpy") and hasattr(value, "item")


def _is_numpy_array(value: Any) -> bool:
    return value.__class__.__module__.startswith("numpy") and hasattr(value, "tolist")


def _is_pandas_dataframe(value: Any) -> bool:
    return (
        value.__class__.__module__.startswith("pandas")
        and value.__class__.__name__ == "DataFrame"
    )


def _is_pandas_series(value: Any) -> bool:
    return (
        value.__class__.__module__.startswith("pandas")
        and value.__class__.__name__ == "Series"
    )


def _is_empty_table_like(value: Any) -> bool:
    if value is None:
        return True
    if _is_pandas_dataframe(value) or _is_pandas_series(value):
        return bool(getattr(value, "empty", False))
    if _is_matlab_table(value):
        return True
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return len(value) == 0
    return False


def _table_like_to_rows(value: Any) -> list[dict[str, Any]]:
    if _is_pandas_dataframe(value):
        return _pandas_dataframe_to_rows(value)
    if _is_pandas_series(value):
        return [to_json_like(value.to_dict())]
    return to_json_like(value)


def _pandas_dataframe_to_rows(value: Any) -> list[dict[str, Any]]:
    if bool(getattr(value, "empty", False)):
        return []
    rows = value.to_dict(orient="records")
    return [to_json_like(row) for row in rows]


def _is_matlab_array(value: Any) -> bool:
    module = value.__class__.__module__
    if not module.startswith("matlab"):
        return False
    return hasattr(value, "size") or hasattr(value, "_data")


def _matlab_array_to_json_like(value: Any) -> Any:
    size = tuple(int(dim) for dim in getattr(value, "size", ()))
    flat = _matlab_flat_values(value)

    if not flat or any(dim == 0 for dim in size):
        return []

    converted = [_matlab_scalar_to_json_like(item) for item in flat]

    if size in ((1, 1), (1,)):
        return converted[0]

    if len(size) == 2:
        rows, cols = size
        if rows == 1 or cols == 1:
            return converted
        return [
            [converted[row + col * rows] for col in range(cols)]
            for row in range(rows)
        ]

    return converted


def _matlab_flat_values(value: Any) -> list[Any]:
    data = getattr(value, "_data", None)
    if data is not None:
        return list(data)
    try:
        return list(value)
    except TypeError:
        return [value]


def _matlab_scalar_to_json_like(value: Any) -> Any:
    if isinstance(value, bool):
        return value
    if isinstance(value, Number):
        return _number_to_json_like(value)
    return to_json_like(value)


def _is_matlab_table(value: Any) -> bool:
    return value.__class__.__name__.lower() == "table"


def _matlab_table_to_rows(value: Any) -> list[dict[str, Any]]:
    # MATLAB Engine exposes tables opaquely in many environments. MATLAB-side
    # bundle fields already include row structs, so this is a conservative stub.
    return []


def _is_json_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))
