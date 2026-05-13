from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import MagicMock


def load_page_module():
    fake_streamlit = types.ModuleType("streamlit")
    fake_streamlit.session_state = {
        "dynamic_results_bundle": {"metadata": {}},
        "results_bundle": {"metadata": {}, "technical": {}},
    }
    fake_streamlit.warning = MagicMock()
    fake_streamlit.stop = MagicMock(side_effect=RuntimeError("st.stop"))
    fake_streamlit.title = MagicMock()
    fake_streamlit.caption = MagicMock()
    fake_streamlit.subheader = MagicMock()
    fake_streamlit.info = MagicMock()
    fake_streamlit.dataframe = MagicMock()
    fake_streamlit.expander = MagicMock()
    fake_streamlit.download_button = MagicMock()
    fake_streamlit.altair_chart = MagicMock()
    fake_streamlit.write = MagicMock()
    fake_streamlit.cache_resource = lambda **_kwargs: (lambda func: func)

    class DummyContext:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    fake_streamlit.expander.return_value = DummyContext()
    sys.modules["streamlit"] = fake_streamlit

    path = Path(__file__).resolve().parents[1] / "pages" / "7_Dynamic_ABM_Compare.py"
    spec = importlib.util.spec_from_file_location("dynamic_abm_compare_page_for_tests", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class DynamicAbmComparePageHelperTests(unittest.TestCase):
    def test_comparison_summary_rows(self) -> None:
        module = load_page_module()
        rows = [
            {"metric": "horizon", "dynamic_value": 20, "abm_value": 20, "comparable": True},
            {"metric": "population", "dynamic_value": 1500, "abm_value": 1500, "comparable": True},
            {"metric": "cumulative cases averted", "dynamic_value": 10, "abm_value": 12, "comparable": True},
        ]
        abm = {
            "technical": {
                "dynamicComparison": {
                    "available": True,
                    "source": "doNothing.derived",
                }
            }
        }

        summary = module.comparison_summary_rows({"metadata": {}}, abm, rows)
        by_field = {row["field"]: row["value"] for row in summary}

        self.assertEqual(by_field["Comparable metrics available"], "3 of 6")
        self.assertEqual(by_field["APY dynamic-comparison available"], True)
        self.assertEqual(by_field["APY dynamic-comparison source"], "doNothing.derived")
        self.assertEqual(by_field["Population match"], "yes")
        self.assertEqual(by_field["Horizon match"], "yes")


if __name__ == "__main__":
    unittest.main()
