from __future__ import annotations

import builtins
from contextlib import contextmanager
import importlib
import json
import sys
from pathlib import Path
from typing import Any, Callable, Iterator


REPO_ROOT = Path(__file__).resolve().parents[1]


def _add_repo_root_to_path() -> None:
    root = str(REPO_ROOT)
    if root not in sys.path:
        sys.path.insert(0, root)


def _matlab_modules() -> list[str]:
    return sorted(
        name for name in sys.modules if name == "matlab" or name.startswith("matlab.")
    )


def _assert_no_matlab_modules_loaded() -> None:
    loaded_matlab = _matlab_modules()
    if loaded_matlab:
        raise AssertionError(f"MATLAB modules loaded: {loaded_matlab}")


@contextmanager
def _guard_against_matlab_imports() -> Iterator[None]:
    for module_name in _matlab_modules():
        sys.modules.pop(module_name, None)

    original_import = builtins.__import__

    def fail_on_matlab_import(name, *args, **kwargs):
        if name == "matlab" or name.startswith("matlab."):
            raise AssertionError(f"attempted to import {name!r}")
        return original_import(name, *args, **kwargs)

    builtins.__import__ = fail_on_matlab_import
    try:
        yield
    finally:
        builtins.__import__ = original_import

    _assert_no_matlab_modules_loaded()


def _json_round_trip(name: str, payload: Any) -> None:
    json.dumps(payload, allow_nan=False, sort_keys=True)
    print(f"PASS {name} is JSON-serialisable")


def _run_check(name: str, check: Callable[[], None]) -> bool:
    try:
        check()
    except Exception as exc:
        print(f"FAIL {name}: {exc}")
        return False

    print(f"PASS {name}")
    return True


def _check_backend_import_without_matlab() -> None:
    with _guard_against_matlab_imports():
        backend_module = importlib.import_module("adapters.python_apy_backend")

        if not hasattr(backend_module, "PythonApyBackend"):
            raise AssertionError("PythonApyBackend was not exported")


def _check_capability_payloads_json() -> None:
    with _guard_against_matlab_imports():
        from adapters.python_apy_backend import PythonApyBackend
        from engine.apy.capabilities import get_apy_capabilities
        from engine.apy.reference_coverage import get_reference_coverage

        backend = PythonApyBackend(REPO_ROOT)

        _json_round_trip(
            "engine.apy.capabilities.get_apy_capabilities()",
            get_apy_capabilities(),
        )
        _json_round_trip(
            "engine.apy.reference_coverage.get_reference_coverage()",
            get_reference_coverage(REPO_ROOT),
        )
        _json_round_trip("PythonApyBackend.capabilities()", backend.capabilities())
        _json_round_trip(
            "PythonApyBackend.reference_coverage()",
            backend.reference_coverage(),
        )


def main() -> int:
    _add_repo_root_to_path()

    checks = (
        ("PythonApyBackend import does not load MATLAB", _check_backend_import_without_matlab),
        ("APY migration JSON payloads", _check_capability_payloads_json),
    )

    print("Running MATLAB-free Python APY migration checks...")
    passed = 0
    for name, check in checks:
        if _run_check(name, check):
            passed += 1

    failed = len(checks) - passed
    if failed:
        print(f"FAILED {failed}/{len(checks)} checks")
        return 1

    print(f"PASSED {passed}/{len(checks)} checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
