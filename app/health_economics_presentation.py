from __future__ import annotations

import math
from typing import Any


def classify_icer_result(incremental_cost: Any, dalys_averted: Any) -> dict[str, Any]:
    inc = _number(incremental_cost)
    dalys = _number(dalys_averted)
    if inc is None or dalys is None or abs(dalys) < 1e-12:
        return {"classification": "ICER not calculable", "icer": None}
    if inc < 0 and dalys > 0:
        return {"classification": "Dominant", "icer": None}
    if inc > 0 and dalys > 0:
        return {"classification": "Higher cost with health gain", "icer": inc / dalys}
    if inc > 0 and dalys <= 0:
        return {"classification": "Dominated", "icer": None}
    if inc < 0 and dalys < 0:
        return {"classification": "Trade-off", "icer": None}
    return {"classification": "No material difference", "icer": None}


def _number(value: Any) -> float | None:
    if value in (None, "", []):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None
