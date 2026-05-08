from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from adapters.backend import JsonDict
from adapters.paths import abm_dir
from adapters.serialization import to_json_like, to_matlab_arg


class MatlabBackend:
    """Thin MATLAB Engine adapter around the frozen APY v9 backend."""

    def __init__(self, root: Path) -> None:
        self.root = Path(root)
        self.abm_path = abm_dir(self.root)
        self._engine: Any | None = None
        self._error = ""

    def status(self) -> JsonDict:
        return {
            "name": "matlab",
            "started": self._engine is not None,
            "abm_path": str(self.abm_path),
            "error": self._error,
        }

    def start(self) -> None:
        if self._engine is not None:
            return
        try:
            import matlab.engine  # type: ignore[import-not-found]

            self._engine = matlab.engine.start_matlab()
            self._engine.addpath(str(self.abm_path), nargout=0)
            self._error = ""
        except Exception as exc:  # pragma: no cover - depends on local MATLAB
            self._error = str(exc)
            raise

    def default_config(self) -> JsonDict:
        return self._call("build_default_config_v9")

    def default_economics_config(self) -> JsonDict:
        return self._call("build_default_economics_config_v9")

    def validate_config(self, config: JsonDict) -> JsonDict:
        return self._call_json("collect_validation_issues_json_v9", config)

    def save_scenario(
        self,
        config: JsonDict,
        path: str,
        economics_config: JsonDict | None = None,
    ) -> JsonDict:
        if economics_config is None:
            return self._call("save_scenario_v9", config, path)
        return self._call("save_scenario_v9", config, path, economics_config)

    def load_scenario(self, path: str) -> tuple[JsonDict, JsonDict, JsonDict]:
        self.start()
        assert self._engine is not None
        config, report, load_info = self._engine.load_scenario_v9(path, nargout=3)
        return (
            to_json_like(config),
            to_json_like(report),
            to_json_like(load_info),
        )

    def run_scenario(self, config: JsonDict) -> JsonDict:
        return self._call("run_scenario_v9", config)

    def build_results_bundle(
        self,
        results: JsonDict,
        validation_report: JsonDict | None = None,
        economics: JsonDict | None = None,
    ) -> JsonDict:
        args: list[Any] = [results]
        if validation_report is not None:
            args.extend(["validationReport", validation_report])
        if economics is not None:
            args.extend(["economics", economics])
        return self._call("build_results_bundle_v9", *args)

    def run_economics(self, results: JsonDict, economics_config: JsonDict) -> JsonDict:
        return self._call("run_health_economics_v9", results, economics_config)

    def run_scenario_bundle(
        self,
        config: JsonDict,
        validation_report: JsonDict | None = None,
    ) -> JsonDict:
        return self._call_json("run_scenario_bundle_json_v9", config)

    def _call(self, function_name: str, *args: Any) -> JsonDict:
        self.start()
        assert self._engine is not None
        function = getattr(self._engine, function_name)
        matlab_args = [to_matlab_arg(self._engine, arg) for arg in args]
        return to_json_like(function(*matlab_args))

    def _call_json(self, function_name: str, *args: Any) -> JsonDict:
        self.start()
        assert self._engine is not None
        function = getattr(self._engine, function_name)
        matlab_args = [to_matlab_arg(self._engine, arg) for arg in args]
        json_text = function(*matlab_args)
        return json.loads(str(json_text))
