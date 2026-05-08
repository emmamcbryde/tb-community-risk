from __future__ import annotations

from typing import Any, Protocol


JsonDict = dict[str, Any]


class Backend(Protocol):
    """Portable app-facing backend contract."""

    def status(self) -> JsonDict:
        ...

    def default_config(self) -> JsonDict:
        ...

    def default_economics_config(self) -> JsonDict:
        ...

    def validate_config(self, config: JsonDict) -> JsonDict:
        ...

    def save_scenario(
        self,
        config: JsonDict,
        path: str,
        economics_config: JsonDict | None = None,
    ) -> JsonDict:
        ...

    def load_scenario(self, path: str) -> tuple[JsonDict, JsonDict, JsonDict]:
        ...

    def run_scenario(self, config: JsonDict) -> JsonDict:
        ...

    def build_results_bundle(
        self,
        results: JsonDict,
        validation_report: JsonDict | None = None,
        economics: JsonDict | None = None,
    ) -> JsonDict:
        ...

    def run_economics(self, results: JsonDict, economics_config: JsonDict) -> JsonDict:
        ...

    def run_scenario_bundle(
        self,
        config: JsonDict,
        validation_report: JsonDict | None = None,
    ) -> JsonDict:
        ...
