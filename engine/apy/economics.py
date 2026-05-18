from __future__ import annotations

from collections.abc import Mapping
from typing import Any


_NOT_PORTED_MESSAGE = (
    "The full APY health-economics model is not yet ported to Python."
)


def calculate_economics(
    result_bundle: Mapping[str, Any],
    economics_config: Mapping[str, Any],
) -> dict[str, Any]:
    """Return the supported Python economics placeholder payload."""
    _validate_result_bundle(result_bundle)
    _validate_economics_config(economics_config)

    return {
        "available": False,
        "message": _NOT_PORTED_MESSAGE,
        "metadata": {
            "model": "apy_economics_python_placeholder",
            "contractVersion": "apy_economics_placeholder_v1",
        },
        "results": {},
    }


def _validate_result_bundle(result_bundle: Mapping[str, Any]) -> None:
    if not isinstance(result_bundle, Mapping):
        raise TypeError("result_bundle must be a mapping.")

    missing = [
        key for key in ("metadata", "results")
        if key not in result_bundle
    ]
    if missing:
        raise ValueError(
            "result_bundle is missing required top-level key(s): "
            + ", ".join(missing)
            + "."
        )

    for key in ("metadata", "results"):
        if not isinstance(result_bundle[key], Mapping):
            raise ValueError(f"result_bundle.{key} must be a mapping.")


def _validate_economics_config(economics_config: Mapping[str, Any]) -> None:
    if not isinstance(economics_config, Mapping):
        raise TypeError("economics_config must be a mapping.")
