from __future__ import annotations

import math

import pandas as pd


def empirical_quantile(values, p: float) -> float:
    x = sorted(float(v) for v in values if not pd.isna(v))
    if not x:
        return math.nan
    n = len(x)
    pos = 1 + (n - 1) * p
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return x[lo - 1]
    return x[lo - 1] + (pos - lo) * (x[hi - 1] - x[lo - 1])


def summarise_numeric_rows(raw_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for column in raw_df.columns:
        series = raw_df[column]
        if not (pd.api.types.is_numeric_dtype(series) or pd.api.types.is_bool_dtype(series)):
            continue
        values = [float(v) for v in series.dropna().tolist()]
        if values:
            median = float(pd.Series(values).median())
            low95 = empirical_quantile(values, 0.025)
            high95 = empirical_quantile(values, 0.975)
        else:
            median = low95 = high95 = math.nan
        rows.append(
            {
                "Metric": column,
                "Median": median,
                "Low95": low95,
                "High95": high95,
            }
        )
    return pd.DataFrame(rows, columns=["Metric", "Median", "Low95", "High95"])
