import unittest
from pathlib import Path

import pandas as pd

import ui.static_ui as static_ui


class _FakeSidebar:
    def __init__(self, run_button=False):
        self.run_button = run_button

    def number_input(self, label, *args, **kwargs):
        if "value" in kwargs:
            return kwargs["value"]
        if args:
            return args[-1]
        return 0

    def slider(self, label, *args, **kwargs):
        if "value" in kwargs:
            return kwargs["value"]
        if len(args) >= 3:
            return args[2]
        return args[-1] if args else 0

    def selectbox(self, label, options, *args, **kwargs):
        index = int(kwargs.get("index", 0))
        return options[index]

    def radio(self, label, options, *args, **kwargs):
        if "Age Distribution" in label or "Choose method" in label:
            return "Default global"
        index = int(kwargs.get("index", 0))
        return options[index]

    def text_input(self, label, value="", *args, **kwargs):
        return value

    def file_uploader(self, *args, **kwargs):
        return None

    def subheader(self, *args, **kwargs):
        return None

    def button(self, label, *args, **kwargs):
        return self.run_button and label == "Run Static Projection"


class _FakeStreamlit:
    def __init__(self, run_button=False, model_choice=None):
        self.sidebar = _FakeSidebar(run_button=run_button)
        self.errors = []
        self.titles = []
        self.radio_options = []
        self.model_choice = model_choice
        self.dataframes = []

    def title(self, text):
        self.titles.append(text)

    def radio(self, label, options, *args, **kwargs):
        self.radio_options.append(list(options))
        return self.model_choice or options[0]

    def subheader(self, *args, **kwargs):
        return None

    def warning(self, *args, **kwargs):
        return None

    def info(self, *args, **kwargs):
        return None

    def success(self, *args, **kwargs):
        return None

    def error(self, message):
        self.errors.append(str(message))

    def dataframe(self, data, *args, **kwargs):
        self.dataframes.append(data)

    def altair_chart(self, *args, **kwargs):
        return None

    def line_chart(self, *args, **kwargs):
        return None

    def write(self, *args, **kwargs):
        return None


class LegacyStaticUITests(unittest.TestCase):
    def test_static_run_path_uses_static_parameters_without_name_error(self):
        fake_st = _FakeStreamlit(run_button=True)
        original_st = static_ui.st
        try:
            static_ui.st = fake_st
            static_ui.render_static_ui()
        finally:
            static_ui.st = original_st

        self.assertEqual(fake_st.errors, [])
        self.assertTrue(any(isinstance(df, pd.DataFrame) for df in fake_st.dataframes))

    def test_static_ui_no_longer_references_removed_dynamic_defaults(self):
        source = static_ui.__loader__.get_source(static_ui.__name__)
        self.assertNotIn("load_dynamic_parameters", source)
        self.assertNotIn("immune_pct", source)

    def test_legacy_app_labels_static_page_as_research_and_development(self):
        fake_st = _FakeStreamlit(model_choice="Legacy static approximation (R&D)")
        calls = []
        app_source = Path("ui/app.py").read_text(encoding="utf-8")
        app_source = app_source.replace("import streamlit as st", "st = injected_st")
        app_source = app_source.replace(
            "from ui.dynamic_ui import render_dynamic_ui",
            "render_dynamic_ui = injected_render_dynamic_ui",
        )
        app_source = app_source.replace(
            "from ui.static_ui import render_static_ui",
            "render_static_ui = injected_render_static_ui",
        )
        namespace = {
            "__file__": str(Path("ui/app.py")),
            "injected_st": fake_st,
            "injected_render_static_ui": lambda: calls.append("static"),
            "injected_render_dynamic_ui": lambda: calls.append("dynamic"),
        }
        exec(compile(app_source, "ui/app.py", "exec"), namespace)

        self.assertEqual(calls, ["static"])
        self.assertIn(
            ["Dynamic transmission model (R&D)", "Legacy static approximation (R&D)"],
            fake_st.radio_options,
        )


if __name__ == "__main__":
    unittest.main()
