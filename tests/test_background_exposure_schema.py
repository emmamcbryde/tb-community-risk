"""Background-exposure (catalytic) configuration schema: validation and isolation from calculations."""

from __future__ import annotations

import json
from pathlib import Path
import unittest

from engine.profiles.background_exposure import (
    BACKGROUND_EXPOSURE_SCHEMA_VERSION,
    HAZARD_UNIT,
    AgeBandHazard,
    BackgroundExposure,
    BackgroundExposureValidationError,
    ExposureMode,
    HazardYear,
    background_exposure_from_profile_payload,
    hazard_value,
    no_background_exposure,
)
from engine.profiles.population_profile import Provenance, ValueState


ROOT = Path(__file__).resolve().parents[1]


def _payload(mode: str, **extra) -> dict:
    return {"schemaVersion": BACKGROUND_EXPOSURE_SCHEMA_VERSION, "mode": mode, **extra}


def _hz(value, **kwargs) -> dict:
    return hazard_value(value, **kwargs).to_dict()


class BackgroundExposureSchemaTests(unittest.TestCase):
    def test_default_and_legacy_profiles_map_to_none(self) -> None:
        from engine.profiles.demonstration import build_demonstration_profile

        profile_payload = build_demonstration_profile().to_dict()
        self.assertNotIn("backgroundExposure", profile_payload)
        config = background_exposure_from_profile_payload(profile_payload)
        self.assertIs(config.mode, ExposureMode.NONE)
        self.assertTrue(config.hazard_is_identically_zero)
        self.assertEqual(config, no_background_exposure())
        self.assertEqual(background_exposure_from_profile_payload({"schemaVersion": "population_profile_v1"}), no_background_exposure())
        self.assertEqual(BackgroundExposure.from_dict(None), no_background_exposure())

    def test_round_trip_and_deterministic_hash(self) -> None:
        configs = [
            no_background_exposure(),
            BackgroundExposure.from_dict(_payload("constant", constantHazard=_hz(0.002), source="Survey", citation="Doe 2020")),
            BackgroundExposure.from_dict(
                _payload(
                    "constant",
                    ageSpecificHazards=[
                        {"ageLower": 15, "ageUpper": None, "hazard": _hz(0.004)},
                        {"ageLower": 0, "ageUpper": 15, "hazard": _hz(0.001)},
                    ],
                )
            ),
            BackgroundExposure.from_dict(
                _payload("time_series", timeSeries=[{"year": 2021, "hazard": _hz(0.003)}, {"year": 2020, "hazard": _hz(None)}])
            ),
        ]
        for config in configs:
            with self.subTest(mode=config.mode.value):
                restored = BackgroundExposure.from_json(config.to_json())
                self.assertEqual(restored, config)
                self.assertEqual(restored.exposure_hash(), config.exposure_hash())
                self.assertEqual(config.to_json(), restored.to_json())
        self.assertEqual(len({c.exposure_hash() for c in configs}), len(configs))

    def test_ordering_is_deterministic(self) -> None:
        rows_a = [{"year": 2020, "hazard": _hz(0.001)}, {"year": 2022, "hazard": _hz(0.002)}, {"year": 2021, "hazard": _hz(0.003)}]
        a = BackgroundExposure.from_dict(_payload("time_series", timeSeries=rows_a))
        b = BackgroundExposure.from_dict(_payload("time_series", timeSeries=list(reversed(rows_a))))
        self.assertEqual([row.year for row in a.time_series], [2020, 2021, 2022])
        self.assertEqual(a.exposure_hash(), b.exposure_hash())
        with self.assertRaises(BackgroundExposureValidationError):
            BackgroundExposure(
                mode=ExposureMode.TIME_SERIES,
                time_series=(HazardYear(2022, hazard_value(0.1)), HazardYear(2021, hazard_value(0.1))),
            )

    def test_zero_is_distinct_from_missing(self) -> None:
        config = BackgroundExposure.from_dict(
            _payload("time_series", timeSeries=[{"year": 2020, "hazard": _hz(0.0)}, {"year": 2021, "hazard": _hz(None)}])
        )
        zero, missing = config.time_series
        self.assertIs(zero.hazard.state, ValueState.VALUE)
        self.assertEqual(zero.hazard.value, 0.0)
        self.assertIs(missing.hazard.state, ValueState.MISSING)
        self.assertIsNone(missing.hazard.value)
        constant_zero = BackgroundExposure.from_dict(_payload("constant", constantHazard=_hz(0.0)))
        self.assertFalse(constant_zero.hazard_is_identically_zero)
        self.assertNotEqual(constant_zero.exposure_hash(), no_background_exposure().exposure_hash())
        with self.assertRaises(BackgroundExposureValidationError):
            BackgroundExposure.from_dict(_payload("constant", constantHazard=_hz(None)))
        with self.assertRaises(BackgroundExposureValidationError):
            BackgroundExposure.from_dict(_payload("time_series", timeSeries=[{"year": 2020, "hazard": _hz(None)}]))

    def test_invalid_values_rejected(self) -> None:
        wrong_unit = {**_hz(0.01), "unit": "per 100,000 population per year"}
        not_finite = [{**_hz(0.01), "value": v} for v in (float("nan"), float("inf"), float("-inf"))]
        cases = {
            "negative": _payload("constant", constantHazard=_hz(-0.001)),
            "wrong hazard unit": _payload("constant", constantHazard=wrong_unit),
            "wrong config unit": _payload("constant", constantHazard=_hz(0.01), unit="per 100,000 population per year"),
            "non-numeric": _payload("constant", constantHazard={**_hz(0.01), "value": "0.01"}),
            "boolean": _payload("constant", constantHazard={**_hz(0.01), "value": True}),
            "duplicate years": _payload("time_series", timeSeries=[{"year": 2020, "hazard": _hz(0.1)}, {"year": 2020, "hazard": _hz(0.2)}]),
            "bad year": _payload("time_series", timeSeries=[{"year": 20.5, "hazard": _hz(0.1)}]),
            "empty series": _payload("time_series", timeSeries=[]),
            "none with values": _payload("none", constantHazard=_hz(0.01)),
            "constant without hazard": _payload("constant"),
            "constant with both": _payload("constant", constantHazard=_hz(0.01), ageSpecificHazards=[{"ageLower": 0, "ageUpper": None, "hazard": _hz(0.01)}]),
            "age gap": _payload("constant", ageSpecificHazards=[{"ageLower": 0, "ageUpper": 10, "hazard": _hz(0.01)}, {"ageLower": 15, "ageUpper": None, "hazard": _hz(0.01)}]),
            "age not from zero": _payload("constant", ageSpecificHazards=[{"ageLower": 5, "ageUpper": None, "hazard": _hz(0.01)}]),
            "age closed": _payload("constant", ageSpecificHazards=[{"ageLower": 0, "ageUpper": 80, "hazard": _hz(0.01)}]),
            "unknown mode": _payload("seasonal"),
            "unknown version": {**_payload("none"), "schemaVersion": "background_exposure_v0"},
            "incidence quantity": _payload("constant", constantHazard=_hz(0.01), quantity="estimated_tb_disease_incidence"),
            "WHO provenance": _payload("constant", constantHazard=_hz(0.01, provenance=Provenance.WHO_SNAPSHOT)),
            "excluded state": _payload("constant", constantHazard={**_hz(0.01), "state": "excluded"}),
        }
        for index, value in enumerate(not_finite):
            cases[f"non-finite {index}"] = _payload("constant", constantHazard=value)
        for name, payload in cases.items():
            with self.subTest(case=name):
                with self.assertRaises(BackgroundExposureValidationError):
                    BackgroundExposure.from_dict(payload)
        with self.assertRaises(BackgroundExposureValidationError):
            BackgroundExposure(mode=ExposureMode.CONSTANT, constant_hazard=hazard_value(-1.0))

    def test_high_hazard_is_valid_with_a_non_blocking_warning(self) -> None:
        high = BackgroundExposure.from_dict(_payload("constant", constantHazard=_hz(5.0)))
        self.assertEqual(high.constant_hazard.value, 5.0)
        self.assertEqual(len(high.plausibility_warnings()), 1)
        self.assertIn("provisional", high.plausibility_warnings()[0])
        series = BackgroundExposure.from_dict(
            _payload("time_series", timeSeries=[{"year": 2020, "hazard": _hz(1.0)}, {"year": 2021, "hazard": _hz(1.5)}])
        )
        self.assertEqual(len(series.plausibility_warnings()), 1)
        self.assertEqual(BackgroundExposure.from_dict(_payload("constant", constantHazard=_hz(0.01))).plausibility_warnings(), ())
        self.assertEqual(no_background_exposure().plausibility_warnings(), ())

    def test_non_finite_json_rejected(self) -> None:
        text = json.dumps(_payload("constant", constantHazard=_hz(0.01))).replace("0.01", "NaN")
        with self.assertRaises(BackgroundExposureValidationError):
            BackgroundExposure.from_json(text)

    def test_user_defined_inputs_keep_provenance(self) -> None:
        config = BackgroundExposure.from_dict(_payload("constant", constantHazard=_hz(0.01, notes="local estimate")))
        restored = BackgroundExposure.from_json(config.to_json())
        self.assertIs(restored.constant_hazard.provenance, Provenance.USER_DEFINED)
        self.assertEqual(restored.constant_hazard.source, "User-defined")
        self.assertEqual(restored.constant_hazard.notes, "local estimate")
        self.assertEqual(restored.constant_hazard.unit, HAZARD_UNIT)
        band = AgeBandHazard(0, None, hazard_value(0.01, provenance=Provenance.LOCAL_UPLOAD, source="upload.csv"))
        self.assertIs(BackgroundExposure(mode=ExposureMode.CONSTANT, age_specific_hazards=(band,)).age_specific_hazards[0].hazard.provenance, Provenance.LOCAL_UPLOAD)


class BackgroundExposureIsolationTests(unittest.TestCase):
    """The schema is scaffolding: nothing that calculates or renders ordinary pages uses it."""

    def test_no_calculation_or_interface_module_imports_the_schema(self) -> None:
        scanned = [
            *sorted((ROOT / "engine").rglob("*.py")),
            *sorted((ROOT / "app").rglob("*.py")),
            *sorted((ROOT / "general_pages").glob("*.py")),
            *sorted((ROOT / "pages").glob("*.py")),
            *sorted((ROOT / "adapters").rglob("*.py")),
            ROOT / "general_app.py",
            ROOT / "streamlit_app.py",
        ]
        own = {ROOT / "engine" / "profiles" / name for name in ("background_exposure.py", "background_exposure_hazard.py")}
        for path in scanned:
            if path in own:
                continue
            with self.subTest(file=path.relative_to(ROOT).as_posix()):
                self.assertNotIn("background_exposure", path.read_text(encoding="utf-8"))

    def test_engine_configuration_and_hash_unchanged(self) -> None:
        from engine.profiles.demonstration import build_demonstration_profile
        from engine.profiles.engine_mapping import build_engine_config, epidemiological_config_hash

        profile = build_demonstration_profile()
        config = build_engine_config(profile)
        text = json.dumps(config, sort_keys=True, default=str)
        for term in ("backgroundExposure", "exogenous_infection_hazard", HAZARD_UNIT):
            self.assertNotIn(term, text)
        self.assertEqual(epidemiological_config_hash(config), epidemiological_config_hash(build_engine_config(profile)))

    def test_no_starsim_dependency(self) -> None:
        for name in ("requirements.txt", "runtime.txt"):
            path = ROOT / name
            if path.exists():
                with self.subTest(file=name):
                    self.assertNotIn("starsim", path.read_text(encoding="utf-8").lower())
        for path in [*sorted((ROOT / "engine").rglob("*.py")), *sorted((ROOT / "app").rglob("*.py"))]:
            with self.subTest(file=path.relative_to(ROOT).as_posix()):
                self.assertNotRegex(path.read_text(encoding="utf-8"), r"^\s*(import|from)\s+starsim")


if __name__ == "__main__":
    unittest.main()
