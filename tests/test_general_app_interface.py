from __future__ import annotations

from pathlib import Path
import io
import re
import socket
import unittest
from unittest.mock import patch
import zipfile

from streamlit.testing.v1 import AppTest

from app.general.terminology import (
    DETERMINISTIC_LABEL,
    RESTORE_DEFAULTS_LABEL,
    STOCHASTIC_LABEL,
    USER_DEFINED_MARK,
    WORKFLOW_PAGES,
    forbidden_matches,
)
from engine.profiles.demonstration import build_demonstration_profile
from engine.profiles.engine_mapping import build_engine_config, epidemiological_config_hash
from engine.profiles.population_profile import with_population_size
from engine.who_incidence.trend import TrendMethod


ROOT = Path(__file__).resolve().parents[1]
PAGES = [ROOT / path for path, _ in WORKFLOW_PAGES]
SETUP = ROOT / "general_pages" / "1_Set_up_population.py"
GENERAL_SOURCES = [
    ROOT / "general_app.py",
    *PAGES,
    *sorted((ROOT / "app" / "general").glob("*.py")),
    *sorted((ROOT / "engine" / "profiles").glob("*.py")),
    *sorted((ROOT / "engine" / "who_incidence").glob("*.py")),
]
# Maintainer-only importer: may name the upstream repository URL, never a local checkout path.
MAINTAINER_ONLY = {"who_import.py"}
ELEMENT_KINDS = (
    "title", "header", "subheader", "markdown", "caption", "info", "warning", "success", "error", "text",
    "button", "selectbox", "radio", "number_input", "slider", "checkbox", "expander", "code", "metric", "multiselect",
)


def _no_network(*args, **kwargs):
    raise AssertionError(f"Network access attempted during page rendering: {args!r}")


def _visible_text(app: AppTest) -> list[str]:
    texts: list[str] = []
    for kind in ELEMENT_KINDS:
        for element in getattr(app, kind, []):
            for attr in ("value", "label", "body", "help"):
                value = getattr(element, attr, None)
                if isinstance(value, str):
                    texts.append(value)
            options = getattr(element, "options", None)
            if options:
                texts.extend(str(option) for option in options)
    for frame in getattr(app, "dataframe", []):
        value = frame.value
        texts.extend(str(column) for column in value.columns)
        texts.extend(str(cell) for cell in value.astype(str).to_numpy().ravel())
    return texts


def _render(path: Path, session: dict | None = None) -> AppTest:
    app = AppTest.from_file(str(path), default_timeout=300)
    for key, value in (session or {}).items():
        app.session_state[key] = value
    with patch.object(socket.socket, "connect", _no_network), patch("socket.create_connection", _no_network):
        app.run(timeout=300)
    return app


def _run(element) -> None:
    with patch.object(socket.socket, "connect", _no_network):
        element.run(timeout=300)


def _country_option(app: AppTest, iso3: str) -> str:
    return next(option for option in app.selectbox(key="general_country_candidate").options if option.endswith(f"({iso3})"))


class GeneralSourceGuardTests(unittest.TestCase):
    def test_general_sources_have_no_setting_specific_user_text_or_network(self) -> None:
        literal = re.compile(r"(['\"])(?:(?!\1).)*\1")
        for path in GENERAL_SOURCES:
            source = path.read_text(encoding="utf-8")
            with self.subTest(file=path.name):
                self.assertNotRegex(source, r"import (requests|urllib|http\.client|socket)\b")
                self.assertNotRegex(source, r"\.\./gtbreport2025|GITHUB[\\/]gtbreport2025")
                if path.name not in MAINTAINER_ONLY:
                    self.assertNotIn("gtbreport2025", source)
                if path in PAGES or path.name == "general_app.py":
                    for match in literal.finditer(source):
                        self.assertEqual(forbidden_matches(match.group(0)), [], match.group(0))


class GeneralPageRenderingTests(unittest.TestCase):
    """Render each ordinary page offline and inspect the visible text."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.rendered = {path.name: _render(path) for path in PAGES}

    def test_pages_render_without_exceptions_or_forbidden_terms(self) -> None:
        for name, app in self.rendered.items():
            with self.subTest(page=name):
                self.assertFalse(app.exception, [e.value for e in app.exception])
                found = [m for text in _visible_text(app) for m in forbidden_matches(text)]
                self.assertEqual(found, [])

    def test_population_defaults_to_10000_with_demonstration_warning(self) -> None:
        app = self.rendered["1_Set_up_population.py"]
        self.assertEqual(app.number_input(key="general_population_input").value, 10_000)
        warnings = " ".join(item.value for item in app.warning)
        self.assertIn("demonstration working defaults, not evidence for any particular country", warnings)
        labels = [button.label for button in app.button]
        self.assertEqual(labels.count(RESTORE_DEFAULTS_LABEL), 1)
        self.assertEqual(app.selectbox(key="general_country_candidate").value, "Select a country or area")

    def test_stochastic_default_is_1000_and_not_launched_without_confirmation(self) -> None:
        app = self.rendered["3_Run_analysis.py"]
        radio = app.radio(key="general_analysis_type")
        self.assertEqual(radio.value, "agent_based")
        self.assertIn(STOCHASTIC_LABEL, radio.options)
        self.assertIn(DETERMINISTIC_LABEL, radio.options)
        self.assertEqual(app.number_input(key="general_n_sims").value, 1_000)
        run = next(button for button in app.button if button.label == "Run analysis")
        self.assertTrue(run.disabled)
        self.assertIn("I understand this analysis takes about", app.checkbox(key="general_confirm_long_run").label)
        self.assertIsNone(app.session_state["general_results_bundle"])

    def test_deterministic_preview_available(self) -> None:
        app = _render(ROOT / "general_pages" / "3_Run_analysis.py")
        _run(app.radio(key="general_analysis_type").set_value("expected_value"))
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["general_analysis"]["analysisMethod"], "expected_value")
        self.assertNotIn("general_n_sims", [item.key for item in app.number_input])
        run = next(button for button in app.button if button.label == "Run analysis")
        self.assertFalse(run.disabled)

    def test_downloads_are_lazy(self) -> None:
        app = _render(ROOT / "general_pages" / "6_Evidence_and_technical_information.py")
        self.assertNotIn("general_provenance_package", app.session_state)
        _run(next(b for b in app.button if b.label == "Prepare provenance package").click())
        package = app.session_state["general_provenance_package"]
        names = set(zipfile.ZipFile(io.BytesIO(package["bytes"])).namelist())
        for expected in ("population_profile.json", "risk_factors.csv", "override_audit.csv", "environment.json", "LIMITATIONS.md", "model_configuration.json"):
            self.assertIn(expected, names)
        self.assertTrue(package["name"].startswith("demonstration_stochastic-1000_no-incidence_"))


class CountryWorkflowTests(unittest.TestCase):
    def test_selection_previews_without_overwriting(self) -> None:
        app = _render(SETUP)
        before = app.session_state["general_profile"]
        _run(app.selectbox(key="general_country_candidate").set_value(_country_option(app, "AUS")))
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["general_profile"]["location"], before["location"])
        self.assertEqual(app.session_state["general_profile"]["incidence"], before["incidence"])
        self.assertIn("Apply selected country incidence data", [b.label for b in app.button])
        texts = " ".join(_visible_text(app))
        self.assertIn("Annual change in estimated incidence", texts)
        self.assertIn("Not changed by country data", texts)
        self.assertIn("They are not yet used to infer infection pressure or transmission", texts)

    def test_preview_survives_page_navigation(self) -> None:
        app = _render(SETUP)
        option = _country_option(app, "PHL")
        _run(app.selectbox(key="general_country_candidate").set_value(option))
        remembered = app.session_state["general_country_candidate_memory"]
        self.assertEqual(remembered, option)
        # Streamlit discards widget state of pages that are not rendered; a return visit keeps only
        # non-widget session keys, which a fresh render with the remembered value reproduces.
        returned = _render(SETUP, {"general_country_candidate_memory": remembered})
        self.assertEqual(returned.selectbox(key="general_country_candidate").value, option)

    def test_apply_changes_only_supported_fields_and_keeps_overrides(self) -> None:
        app = _render(SETUP)
        _run(app.number_input(key="general_population_input").set_value(5000))
        demo = app.session_state["general_profile"]
        _run(app.selectbox(key="general_country_candidate").set_value(_country_option(app, "ZAF")))
        _run(next(b for b in app.button if b.label == "Apply selected country incidence data").click())
        self.assertFalse(app.exception)
        after = app.session_state["general_profile"]
        self.assertEqual(after["location"]["iso3"], "ZAF")
        self.assertEqual(after["incidence"]["provenance"], "who_snapshot")
        self.assertTrue(after["incidence"]["snapshotId"].startswith("who-gtb2025-incidence-"))
        for key in ("ltbiPrevalence", "ageDistribution", "riskFactors", "populationSize"):
            self.assertEqual(after[key], demo[key], key)
        self.assertEqual(after["populationSize"]["provenance"], "user_defined")
        self.assertIn("incidence data are applied to this profile", " ".join(item.value for item in app.success))

    def test_trend_controls_change_only_trend_settings(self) -> None:
        app = _render(SETUP)
        _run(app.selectbox(key="general_country_candidate").set_value(_country_option(app, "IDN")))
        _run(app.selectbox(key="general_trend_window").set_value(5))
        self.assertEqual(app.session_state["general_trend_settings"]["window_years"], 5)
        _run(app.radio(key="general_trend_method").set_value(TrendMethod.PENALISED_SPLINE))
        self.assertEqual(app.session_state["general_trend_settings"]["method"], "penalised_spline")
        _run(app.selectbox(key="general_trend_covid").set_value("exclude"))
        self.assertEqual(app.session_state["general_trend_settings"]["covid_handling"], "exclude")
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["general_profile"]["incidence"]["provenance"], "bundled")

    def test_restore_clears_country_and_overrides(self) -> None:
        app = _render(SETUP)
        _run(app.number_input(key="general_population_input").set_value(4000))
        _run(app.selectbox(key="general_country_candidate").set_value(_country_option(app, "AUS")))
        _run(next(b for b in app.button if b.label == "Apply selected country incidence data").click())
        self.assertIn(USER_DEFINED_MARK, " ".join(_visible_text(app)))
        _run(next(b for b in app.button if b.label == RESTORE_DEFAULTS_LABEL).click())
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["general_profile"], build_demonstration_profile().to_dict())
        self.assertEqual(app.number_input(key="general_population_input").value, 10_000)
        self.assertEqual(app.selectbox(key="general_country_candidate").value, "Select a country or area")
        self.assertNotIn(USER_DEFINED_MARK + " values in this profile", " ".join(_visible_text(app)))

    def test_conflict_requires_choice_for_local_incidence(self) -> None:
        from engine.profiles.local_incidence import apply_local_incidence, parse_local_incidence
        from engine.profiles.country import USE_NEW

        upload = parse_local_incidence(
            b"location,year,measure,incidence_per_100k,lower,upper,source\nNorth,2022,estimated_incidence,50,40,60,Local survey\n"
            b"North,2023,estimated_incidence,48,39,58,Local survey\n",
            filename="north.csv",
        )
        local = apply_local_incidence(build_demonstration_profile(), upload, resolutions={"incidence": USE_NEW})
        app = _render(SETUP, {"general_profile": local.to_dict()})
        _run(app.selectbox(key="general_country_candidate").set_value(_country_option(app, "AUS")))
        apply = next(b for b in app.button if b.label == "Apply selected country incidence data")
        self.assertTrue(apply.disabled)
        _run(app.radio(key="general_conflict_incidence").set_value("keep_current"))
        _run(app.radio(key="general_conflict_location").set_value("keep_current"))
        _run(next(b for b in app.button if b.label == "Apply selected country incidence data").click())
        after = app.session_state["general_profile"]
        self.assertEqual(after["incidence"]["provenance"], "local_upload")
        self.assertEqual(after["location"]["name"], "North")

    def test_country_with_no_estimates_and_short_series(self) -> None:
        from engine.who_incidence import snapshot as snapshot_module

        fixture = snapshot_module.find_manifests(snapshot_module.FIXTURE_DIR)[0]
        with patch.object(snapshot_module, "default_manifest_path", return_value=fixture):
            snapshot_module._load_cached.cache_clear()
            app = _render(SETUP)
            _run(app.selectbox(key="general_country_candidate").set_value(_country_option(app, "PRK")))
            self.assertIn("publishes no incidence estimates", " ".join(item.value for item in app.warning))
            self.assertNotIn("Apply selected country incidence data", [b.label for b in app.button])
            _run(app.selectbox(key="general_country_candidate").set_value(_country_option(app, "ANT")))
            self.assertIn("the series is incomplete", " ".join(item.value for item in app.warning))
            self.assertIn("offline example subset", " ".join(item.value for item in app.caption))
        snapshot_module._load_cached.cache_clear()

    def test_missing_snapshot_gives_clear_message(self) -> None:
        from engine.who_incidence import snapshot as snapshot_module

        with patch.object(snapshot_module, "default_manifest_path", side_effect=snapshot_module.SnapshotUnavailable("No WHO incidence snapshot is installed under data/who_incidence/.")):
            app = _render(SETUP)
        self.assertFalse(app.exception)
        self.assertIn("Country data are unavailable", " ".join(item.value for item in app.error))
        self.assertEqual(app.number_input(key="general_population_input").value, 10_000)


class GeneralEndToEndTests(unittest.TestCase):
    """One small offline run of a profile without risk factors, then results pages."""

    @classmethod
    def setUpClass(cls) -> None:
        from adapters.python_apy_backend import PythonApyBackend
        from engine.profiles.demonstration import demonstration_profile_without_risk_factors
        from engine.profiles.population_profile import with_population_size

        profile = with_population_size(demonstration_profile_without_risk_factors(), 200)
        cls.profile = profile
        cls.analysis = {"analysisMethod": "agent_based", "nReps": 2, "seed": 1}
        cls.config = build_engine_config(profile, analysis=cls.analysis)
        backend = PythonApyBackend(ROOT)
        with patch.object(socket.socket, "connect", _no_network):
            cls.bundle = backend.run_scenario_bundle(cls.config, validation_report=backend.validate_config(cls.config))
        det_config = build_engine_config(profile, analysis={"analysisMethod": "expected_value", "nReps": 1000, "seed": 1})
        cls.det_config = det_config
        cls.det_bundle = backend.run_scenario_bundle(det_config, validation_report=backend.validate_config(det_config))

    def _session(self, bundle, config, analysis):
        return {
            "general_profile": self.profile.to_dict(),
            "general_analysis": analysis,
            "general_results_bundle": bundle,
            "general_results_config": config,
            "general_results_epi_hash": epidemiological_config_hash(config),
        }

    def test_profile_without_risk_factors_runs(self) -> None:
        self.assertTrue(self.bundle["headline"]["keyMetricsRows"])
        link = self.bundle["technical"]["interfaceConfig"]["generalProfileLink"]
        self.assertEqual(link["profileHash"], self.profile.profile_hash())

    def test_results_current_then_stale_after_epidemiological_change(self) -> None:
        session = self._session(self.bundle, self.config, self.analysis)
        results = _render(ROOT / "general_pages" / "4_Results.py", session)
        self.assertFalse(results.exception)
        self.assertIn("Results are current", " ".join(item.value for item in results.success))
        self.assertEqual([m for t in _visible_text(results) for m in forbidden_matches(t)], [])
        self.assertIn("simulation interval", " ".join(item.value for item in results.caption))
        stale = dict(session)
        stale["general_profile"] = with_population_size(self.profile, 300).to_dict()
        results = _render(ROOT / "general_pages" / "4_Results.py", stale)
        self.assertIn("out of date", " ".join(item.value for item in results.warning))

    def test_deterministic_results_show_no_stochastic_interval(self) -> None:
        analysis = {"analysisMethod": "expected_value", "nReps": 1000, "seed": 1}
        results = _render(ROOT / "general_pages" / "4_Results.py", self._session(self.det_bundle, self.det_config, analysis))
        self.assertFalse(results.exception)
        cells = " ".join(_visible_text(results))
        self.assertIn("N/A", cells)
        economics = _render(ROOT / "general_pages" / "5_Health_economics.py", self._session(self.det_bundle, self.det_config, analysis))
        _run(next(b for b in economics.button if b.label == "Calculate health economics").click())
        frame = economics.dataframe[-1].value
        self.assertEqual(set(frame["Low 95%"]), {"N/A"})

    def test_economics_paired_cloud_and_cost_only_changes(self) -> None:
        from app.general.economics import apply_cost_edits, default_economics_config, icer_cloud_points, run_economics

        session = self._session(self.bundle, self.config, self.analysis)
        economics = _render(ROOT / "general_pages" / "5_Health_economics.py", session)
        with patch("adapters.python_apy_backend.PythonApyBackend.run_scenario_bundle", side_effect=AssertionError("epidemiology rerun")):
            _run(next(b for b in economics.button if b.label == "Calculate health economics").click())
        self.assertFalse(economics.exception)
        visible = _visible_text(economics)
        self.assertIn("DALYs averted", visible)
        text = " ".join(visible).lower()
        self.assertNotIn("reference", text)
        self.assertEqual([m for t in visible for m in forbidden_matches(t)], [])
        self.assertIn("Cost-effectiveness plane", visible)

        base = run_economics(self.bundle, default_economics_config())
        edited = run_economics(self.bundle, apply_cost_edits(default_economics_config(), {"test_igra": 500.0}))
        before, after = icer_cloud_points(base), icer_cloud_points(edited)
        self.assertEqual([p["dalysAverted"] for p in before], [p["dalysAverted"] for p in after])
        self.assertTrue(all(b["incrementalCost"] < a["incrementalCost"] for b, a in zip(before, after)))

    def test_effect_warnings_do_not_change_configuration(self) -> None:
        from engine.profiles.effect_measures import effect_warnings

        profile = build_demonstration_profile()
        before = build_engine_config(profile)
        warnings = effect_warnings(profile)
        after = build_engine_config(profile)
        self.assertTrue({w.code for w in warnings} >= {"or_as_hazard_multiplier", "multiplied_factors", "extreme_combined_multiplier"})
        self.assertEqual(epidemiological_config_hash(before), epidemiological_config_hash(after))


if __name__ == "__main__":
    unittest.main()
