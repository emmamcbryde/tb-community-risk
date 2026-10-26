"""Presentation-only formatting; never use these strings as model inputs."""
from __future__ import annotations

from decimal import Decimal, ROUND_HALF_UP
from numbers import Real
from typing import Iterable

import pandas as pd


def format_number(value: Real, *, decimals: int = 2, minimum_decimals: int = 0) -> str:
    text = f"{float(value):,.{decimals}f}"
    if decimals > minimum_decimals:
        whole, fraction = text.split(".")
        fraction = fraction.rstrip("0").ljust(minimum_decimals, "0")
        text = whole + ("." + fraction if fraction else "")
    return "0" if text == "-0" else text


def format_percent(proportion: Real, *, decimals: int = 1, trim: bool = False) -> str:
    return format_number(
        100.0 * float(proportion), decimals=decimals,
        minimum_decimals=0 if trim else decimals,
    ) + "%"


def percentage_input_format(proportion: Real, *, baseline: bool = False) -> str:
    # Preserve a meaningful second decimal for small baseline active-TB inputs.
    percent = 100.0 * float(proportion)
    return "%.2f" if baseline and abs(percent - round(percent, 1)) > 1e-10 else "%.1f"


def format_aud(value: Real | Decimal, *, per_unit: str | None = None) -> str:
    """Future currency contract only; this module creates no economic outputs.

    Round halves away from zero; an authorised future arithmetic ICER can pass
    per_unit="DALY" to use the same whole-dollar presentation rule.
    """
    rounded = Decimal(str(value)).quantize(Decimal("1"), rounding=ROUND_HALF_UP)
    sign = "-" if rounded < 0 else ""
    text = f"{sign}AUD {abs(rounded):,.0f}"
    return text + (f" per {per_unit}" if per_unit else "")


def display_table(
    rows: pd.DataFrame | Iterable[dict], *,
    percentage_columns: Iterable[str] = (),
    integer_columns: Iterable[str] = (),
) -> pd.DataFrame:
    """Return a formatted copy, retaining the source's full precision."""
    table = pd.DataFrame(rows).copy(deep=True)
    percentages, integers = set(percentage_columns), set(integer_columns)
    for column in table.columns:
        def present(value):
            if not isinstance(value, Real) or isinstance(value, bool):
                return value
            if pd.isna(value):
                return "N/A"
            if column in percentages:
                return format_percent(value)
            return format_number(value, decimals=0 if column in integers else 2)
        table[column] = table[column].map(present)
    return table
