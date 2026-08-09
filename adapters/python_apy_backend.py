from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

from adapters.backend import JsonDict
from adapters.serialization import to_json_like
from engine.apy.config import build_default_config
from engine.apy.economics import (
    build_default_economics_config,
    build_economics_preset_kwab150,
    run_health_economics,
    run_health_economics_for_config,
)
from engine.apy.results_bundle import build_results_bundle
from engine.apy.runner import run_scenario, run_scenario_with_do_nothing
from engine.apy.validation import collect_validation_issues


class PythonApyBackend:
    """Experimental pure-Python APY v9 backend adapter.

    This adapter intentionally does not import or call MATLAB. MATLAB remains
    the reference backend while Python parity validation is expanded.
    """

    def __init__(self, root: Path) -> None:
        self.root = Path(root)

    def status(self) -> JsonDict:
        return {
            "name": "python_apy",
            "started": True,
            "abm_path": "",
            "error": "",
            "experimental": True,
            "matlabRequired": False,
        }

    def default_config(self) -> JsonDict:
        return to_json_like(build_default_config())

    def validate_config(self, config: JsonDict) -> JsonDict:
        return to_json_like(collect_validation_issues(_matlab_empty_to_none(config)))

    def run_scenario(self, config: JsonDict) -> JsonDict:
        return to_json_like(run_scenario(_matlab_empty_to_none(config)))

    def build_results_bundle(
        self,
        results: JsonDict,
        validation_report: JsonDict | None = None,
        economics: JsonDict | None = None,
    ) -> JsonDict:
        bundle = build_results_bundle(results)
        if validation_report is not None:
            bundle["validation"] = {"report": validation_report}
        if economics is not None:
            bundle["economics"] = economics
        return to_json_like(bundle)

    def run_scenario_bundle(
        self,
        config: JsonDict,
        validation_report: JsonDict | None = None,
        progress_callback: Callable[[dict[str, Any]], None] | None = None,
    ) -> JsonDict:
        out = run_scenario_with_do_nothing(
            _matlab_empty_to_none(config),
            progress_callback=progress_callback,
        )
        bundle = out["bundle"]
        if validation_report is not None:
            bundle["validation"] = {"report": validation_report}
        return to_json_like(bundle)

    def save_scenario(
        self,
        config: JsonDict,
        path: str,
        economics_config: JsonDict | None = None,
    ) -> JsonDict:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "contractVersion": "ltbi_screening_scenario_v1",
            "backend": "python_apy",
            "scenarioLabel": config.get("scenarioLabel", ""),
            "scenario": to_json_like(config.get("scenario")),
            "config": to_json_like(config),
            "economics": to_json_like(economics_config),
        }
        target.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return {
            "filename": str(target),
            "saved": True,
            "backend": "python_apy",
        }

    def load_scenario(self, path: str) -> tuple[JsonDict, JsonDict, JsonDict]:
        source = Path(path)
        payload = json.loads(source.read_text(encoding="utf-8"))
        config = payload.get("config")
        if not isinstance(config, dict):
            raise ValueError("Scenario JSON does not contain a config object.")
        report = self.validate_config(config)
        load_info = {
            "filename": str(source),
            "contractVersion": payload.get("contractVersion", ""),
            "scenarioLabel": payload.get("scenarioLabel", ""),
            "backend": payload.get("backend", ""),
        }
        economics = payload.get("economics")
        return to_json_like(config), report, {**load_info, "economics": economics}

    def default_economics_config(self) -> JsonDict:
        return to_json_like(build_default_economics_config())

    def economics_preset_kwab150(self) -> JsonDict:
        return to_json_like(build_economics_preset_kwab150())

    def run_economics(self, results: JsonDict, economics_config: JsonDict) -> JsonDict:
        return to_json_like(
            run_health_economics(
                _matlab_empty_to_none(results),
                _matlab_empty_to_none(economics_config),
            )
        )

    def run_economics_for_config(
        self,
        config: JsonDict,
        economics_config: JsonDict,
    ) -> JsonDict:
        return to_json_like(
            run_health_economics_for_config(
                _matlab_empty_to_none(config),
                _matlab_empty_to_none(economics_config),
            )
        )


def _matlab_empty_to_none(value):
    if value == []:
        return None
    if isinstance(value, dict):
        return {
            key: _matlab_empty_to_none(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_matlab_empty_to_none(item) for item in value]
    return value
