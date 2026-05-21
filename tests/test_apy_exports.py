from __future__ import annotations

import json
import math
import unittest

import numpy as np

from engine.apy.exports import (
    UNSUPPORTED_CHART_EXPORT_MESSAGE,
    UNSUPPORTED_HEADLINE_KEY_METRICS_MESSAGE,
    UNSUPPORTED_HEADLINE_SUMMARY_MESSAGE,
    UNSUPPORTED_SUMMARY_EXPORT_MESSAGE,
    chart_numeric_series,
    headline_display_tables,
    json_export_payload,
    summary_csv_export,
    summary_rows_for_csv,
)


SUMMARY_COLUMNS_FOR_TEST = ["Metric", "Median", "Low95", "High95"]


class ApyExportTests(unittest.TestCase):
    def test_json_export_payload_is_json_serialisable(self) -> None:
        bundle = {
            "metadata": {"scenarioLabel": "Example", "replicates": np.int64(3)},
            "headline": {
                "summaryRows": [
                    {"Metric": "nScreened", "Median": np.float64(10.5), "Low95": math.nan}
                ]
            },
        }

        payload = json_export_payload(bundle)

        self.assertEqual(payload["metadata"]["replicates"], 3)
        self.assertIsNone(payload["headline"]["summaryRows"][0]["Low95"])
        json.dumps(payload)

    def test_summary_rows_for_csv_are_stable_and_csv_ready(self) -> None:
        bundle = {
            "headline": {
                "summaryRows": [
                    {
                        "High95": 13.0,
                        "Metric": "nScreened",
                        "Median": 10.0,
                        "Low95": 7.0,
                        "Notes": "screened",
                    },
                    {
                        "Metric": "nPreventedActiveTB",
                        "Median": 2.0,
                        "Low95": 1.0,
                        "High95": 3.0,
                        "Notes": "prevented",
                    },
                ]
            }
        }

        rows = summary_rows_for_csv(bundle)

        self.assertEqual(
            list(rows[0].keys()),
            ["Metric", "Median", "Low95", "High95", "Notes"],
        )
        self.assertEqual(rows[0]["Metric"], "nScreened")
        self.assertEqual(rows[1]["Median"], 2.0)

    def test_summary_csv_export_is_explicitly_unsupported_without_rows(self) -> None:
        status = summary_csv_export({"headline": {}})

        self.assertFalse(status["available"])
        self.assertEqual(status["message"], UNSUPPORTED_SUMMARY_EXPORT_MESSAGE)
        self.assertEqual(status["rows"], [])

    def test_headline_display_tables_are_stable_and_json_safe(self) -> None:
        bundle = {
            "headline": {
                "keyMetricsRows": [
                    {
                        "High95": np.float64(13.0),
                        "Metric": "nScreened",
                        "Median": np.float64(10.0),
                        "Low95": math.nan,
                        "Notes": "screened",
                    }
                ],
                "summaryRows": [
                    {
                        "Metric": "nPreventedActiveTB",
                        "Median": 2.0,
                        "Low95": 1.0,
                        "High95": 3.0,
                    }
                ],
            }
        }

        tables = headline_display_tables(bundle)

        self.assertEqual(
            tables["keyMetricsRows"]["columns"],
            ["Metric", "Median", "Low95", "High95", "Notes"],
        )
        self.assertTrue(tables["keyMetricsRows"]["available"])
        self.assertIsNone(tables["keyMetricsRows"]["rows"][0]["Low95"])
        self.assertEqual(tables["summaryRows"]["columns"], SUMMARY_COLUMNS_FOR_TEST)
        json.dumps(tables)

    def test_headline_display_tables_are_explicitly_unavailable_without_rows(self) -> None:
        tables = headline_display_tables({"headline": {"keyMetricsRows": []}})

        self.assertFalse(tables["keyMetricsRows"]["available"])
        self.assertEqual(
            tables["keyMetricsRows"]["message"],
            UNSUPPORTED_HEADLINE_KEY_METRICS_MESSAGE,
        )
        self.assertEqual(tables["keyMetricsRows"]["columns"], [])
        self.assertEqual(tables["keyMetricsRows"]["rows"], [])
        self.assertFalse(tables["summaryRows"]["available"])
        self.assertEqual(
            tables["summaryRows"]["message"],
            UNSUPPORTED_HEADLINE_SUMMARY_MESSAGE,
        )

    def test_chart_numeric_series_uses_existing_plot_data_rows(self) -> None:
        bundle = {
            "charts": {
                "plotDataRows": [
                    {
                        "screeningStrategy": "cure",
                        "screenCoverage": 0.1,
                        "meanCuredInfection": 5.0,
                        "meanPreventedActiveTB": 1.5,
                    },
                    {
                        "screeningStrategy": "cure",
                        "screenCoverage": 0.2,
                        "meanCuredInfection": 8.0,
                        "meanPreventedActiveTB": 2.0,
                    },
                ]
            }
        }

        status = chart_numeric_series(bundle)

        self.assertTrue(status["available"])
        series_by_name = {series["name"]: series for series in status["series"]}
        self.assertEqual(
            set(series_by_name),
            {"screenCoverage", "meanCuredInfection", "meanPreventedActiveTB"},
        )
        self.assertEqual(series_by_name["meanCuredInfection"]["points"][0]["value"], 5.0)
        self.assertEqual(
            series_by_name["meanCuredInfection"]["points"][0]["row"],
            {"screeningStrategy": "cure"},
        )

    def test_chart_numeric_series_is_explicitly_unsupported_without_plot_data(self) -> None:
        status = chart_numeric_series({"charts": {"available": False}})

        self.assertFalse(status["available"])
        self.assertEqual(status["message"], UNSUPPORTED_CHART_EXPORT_MESSAGE)
        self.assertEqual(status["series"], [])


if __name__ == "__main__":
    unittest.main()
