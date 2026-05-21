from __future__ import annotations

import json

from app.display import (
    economics_assumptions_rows,
    economics_summary_csv,
    economics_summary_rows,
)


def test_economics_assumptions_rows_have_stable_order_and_labels() -> None:
    rows = economics_assumptions_rows(
        {
            "metadata": {
                "currencyCode": "AUD",
                "priceYear": 2024,
                "locationLabel": "APY Lands",
            },
            "costs": {
                "test": {"IGRA": 65, "TST": 35},
                "regimen": {
                    "x3HP": 120,
                    "x4R": 110,
                    "x3HR": 100,
                    "x6H": 90,
                    "x9H": 80,
                },
                "falsePositiveIncrementalPerPerson": 12,
                "activeTBDiseasePerCase": 10000,
                "programSetupTotal": 5000,
                "programRunningTotal": 7000,
            },
        }
    )

    assert [row["field"] for row in rows] == [
        "currencyCode",
        "priceYear",
        "locationLabel",
        "test.IGRA",
        "test.TST",
        "regimen.3HP",
        "regimen.4R",
        "regimen.3HR",
        "regimen.6H",
        "regimen.9H",
        "falsePositiveIncrementalPerPerson",
        "activeTBDiseasePerCase",
        "programSetupTotal",
        "programRunningTotal",
    ]
    assert all(list(row) == ["field", "value"] for row in rows)


def test_economics_summary_rows_have_stable_column_order() -> None:
    rows = economics_summary_rows(
        {
            "summaryRows": [
                {
                    "metric": "totalImplementedCost",
                    "label": "Total implemented cost",
                    "value": 123.45,
                    "category": "directCost",
                    "status": "implemented",
                    "includedInTotal": True,
                    "source": "calculated",
                },
                {
                    "metric": "costPerTBCasePrevented",
                    "label": "Cost per TB case prevented",
                    "value": 67.89,
                    "status": "implemented",
                    "note": "denominator supplied",
                },
            ]
        }
    )

    assert [list(row) for row in rows] == [
        [
            "label",
            "value",
            "category",
            "status",
            "includedInTotal",
            "metric",
            "source",
        ],
        ["label", "value", "status", "metric", "note"],
    ]


def test_economics_summary_csv_uses_deterministic_summary_rows() -> None:
    payload = {
        "summaryRows": [
            {
                "metric": "totalImplementedCost",
                "value": 123.45,
                "label": "Total implemented cost",
            },
            {
                "value": 67.89,
                "label": "Cost per TB case prevented",
                "metric": "costPerTBCasePrevented",
            },
        ]
    }

    first = economics_summary_csv(payload)
    second = economics_summary_csv(payload["summaryRows"])

    assert first == second
    assert first.decode("utf-8").splitlines() == [
        "label,value,metric",
        "Total implemented cost,123.45,totalImplementedCost",
        "Cost per TB case prevented,67.89,costPerTBCasePrevented",
    ]


def test_economics_display_rows_are_json_serialisable() -> None:
    assumptions = economics_assumptions_rows(
        {
            "metadata": {"currencyCode": "AUD", "priceYear": []},
            "costs": {"test": {"IGRA": []}},
        }
    )
    summary = economics_summary_rows(
        [
            {
                "label": "Total implemented cost",
                "value": 123.45,
                "includedInTotal": True,
                "metric": "totalImplementedCost",
            }
        ]
    )

    json.dumps({"assumptions": assumptions, "summary": summary}, sort_keys=True)
