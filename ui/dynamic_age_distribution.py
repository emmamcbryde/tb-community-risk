from __future__ import annotations

import pandas as pd


def default_five_year_age_distribution() -> pd.DataFrame:
    rows = []
    for start in range(0, 100, 5):
        rows.append(
            {
                "AgeGroup": f"{start}-{start + 4}",
                "AgeStart": start,
                "AgeEnd": start + 4,
                "Proportion": 5 / 101,
            }
        )
    rows.append(
        {
            "AgeGroup": "100+",
            "AgeStart": 100,
            "AgeEnd": 100,
            "Proportion": 1 / 101,
        }
    )
    return pd.DataFrame(rows)


def default_single_year_age_distribution() -> pd.DataFrame:
    return pd.DataFrame(
        {"AgeGroup": list(range(0, 101)), "Proportion": [1 / 101] * 101}
    )


def _looks_like_five_year_table(df: pd.DataFrame) -> bool:
    return {"AgeStart", "AgeEnd", "Proportion"}.issubset(df.columns)


def normalise_age_distribution(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or "Proportion" not in df.columns:
        return default_single_year_age_distribution()

    out = df.copy()
    out["Proportion"] = pd.to_numeric(out["Proportion"], errors="coerce").fillna(0.0)
    out["Proportion"] = out["Proportion"].clip(lower=0.0)

    total = float(out["Proportion"].sum())
    if total > 0:
        out["Proportion"] = out["Proportion"] / total
        return out

    if _looks_like_five_year_table(out):
        return default_five_year_age_distribution()
    return default_single_year_age_distribution()


def expand_five_year_age_distribution(
    df: pd.DataFrame, population: float
) -> pd.DataFrame:
    df = normalise_age_distribution(df)
    if not _looks_like_five_year_table(df):
        return pd.DataFrame(
            {
                "age": range(0, 101),
                "population": [float(population) / 101] * 101,
            }
        )

    rows = []
    for _, row in df.iterrows():
        start = int(row["AgeStart"])
        end = int(row["AgeEnd"])
        ages = list(range(start, end + 1))
        if not ages:
            continue
        per_age_population = float(row["Proportion"]) * float(population) / len(ages)
        for age in ages:
            rows.append({"age": age, "population": per_age_population})
    return pd.DataFrame(rows)
