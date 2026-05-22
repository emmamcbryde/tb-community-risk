from __future__ import annotations

from pathlib import Path
from typing import Any

from engine.apy.capabilities import (
    COMPLETE,
    DOCUMENTATION_ONLY,
    MATLAB_REFERENCE,
    PARTIAL,
    REFERENCE_ONLY,
    UNSUPPORTED,
    capability_rows,
)


CONTRACT_VERSION = "apy_reference_coverage_v1"

_REFERENCE_ROOT = Path("validation") / "matlab_reference"
_REFERENCE_SUITE = _REFERENCE_ROOT / "scenario_suite_v1.json"
_TESTS_ROOT = Path("tests")

_PARITY_TEST_FILES = {"test_apy_parity.py"}
_DISTRIBUTIONAL_TEST_FILES = {"test_apy_distributional_validation.py"}

_CORE_REFERENCE_ELIGIBLE_PARTIAL_COMPONENTS = {
    "core_simulation",
    "scenario_runner",
}

_NEVER_COMPLETE_COMPONENTS = {
    "attributable_risk",
    "attributable_risk_analysis",
    "chart_rendering_parity",
    "do_nothing_dynamic_comparison",
    "economics",
    "export_display_helpers",
    "exports",
    "health_economics",
    "matlab_app_callbacks",
    "matlab_app_designer_helpers",
    "natural_history",
    "natural_history_addon_reporting",
    "targeting_compare",
    "targeting_strategy_comparison",
}

_PYTHON_ONLY_TEST_FILES_BY_COMPONENT = {
    "age_distribution": {"test_apy_data_loading.py"},
    "attributable_risk": {
        "test_apy_attributable_risk.py",
        "test_python_apy_backend_attributable_risk.py",
    },
    "attributable_risk_analysis": {
        "test_apy_attributable_risk.py",
        "test_python_apy_backend_attributable_risk.py",
    },
    "calibration": {"test_apy_calibration.py"},
    "chart_rendering_parity": {"test_apy_exports.py"},
    "cohort_primitives": {
        "test_apy_cohort.py",
        "test_apy_simulation_one_cohort.py",
    },
    "config_validation": {"test_apy_config_validation.py"},
    "core_simulation": {
        "test_apy_runner.py",
        "test_apy_simulation_one_cohort.py",
    },
    "data_loading": {"test_apy_data_loading.py"},
    "default_config": {"test_apy_config_validation.py"},
    "do_nothing_dynamic_comparison": {
        "test_apy_natural_history.py",
        "test_apy_results_bundle.py",
        "test_apy_runner.py",
    },
    "economics": {"test_apy_economics.py", "test_apy_economics_compare.py"},
    "export_display_helpers": {"test_apy_exports.py"},
    "exports": {"test_apy_exports.py"},
    "health_economics": {"test_apy_economics.py", "test_apy_economics_compare.py"},
    "natural_history": {"test_apy_natural_history.py"},
    "natural_history_addon_reporting": {"test_apy_natural_history.py"},
    "regimen_handling": {"test_apy_regimen.py"},
    "results_bundle": {"test_apy_results_bundle.py"},
    "scenario_runner": {"test_apy_runner.py"},
    "scenario_save_load": {"test_python_apy_backend.py"},
    "summary_rows": {"test_apy_summary.py"},
    "targeting_compare": {"test_apy_targeting.py"},
    "targeting_strategy_comparison": {"test_apy_targeting.py"},
}

_PARITY_COMPONENTS = {
    "core_simulation",
    "do_nothing_dynamic_comparison",
    "matlab_reference_fixtures",
    "scenario_runner",
}

_DISTRIBUTIONAL_COMPONENTS = {
    "core_simulation",
    "distributional_validation",
    "do_nothing_dynamic_comparison",
    "matlab_reference_fixtures",
    "scenario_runner",
}

_REFERENCE_COMPONENTS = {
    "distributional_validation",
    "matlab_reference_engine",
    "matlab_reference_fixtures",
    "matlab_user_options_compatibility",
}


def get_reference_coverage(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Return JSON-serialisable static APY reference/test coverage metadata.

    The helper intentionally does not import MATLAB or execute tests. It only
    checks committed reference-fixture files and statically visible Python test
    files, then combines those signals with the capability contract.
    """
    root = _resolve_repo_root(repo_root)
    reference_root = root / _REFERENCE_ROOT
    suite_path = root / _REFERENCE_SUITE
    tests_root = root / _TESTS_ROOT

    reference_fixture_dirs = _reference_fixture_dirs(reference_root)
    suite_exists = suite_path.is_file()
    test_files = _test_files(tests_root)

    components = [
        _coverage_row(
            capability=row,
            suite_exists=suite_exists,
            fixture_dirs_exist=bool(reference_fixture_dirs),
            test_files=test_files,
        )
        for row in capability_rows()
    ]

    return {
        "contractVersion": CONTRACT_VERSION,
        "matlabRequired": False,
        "repoRoot": str(root),
        "referenceSuite": str(_REFERENCE_SUITE),
        "referenceSuiteExists": suite_exists,
        "referenceFixtureCount": len(reference_fixture_dirs),
        "referenceFixtureDirs": [path.name for path in reference_fixture_dirs],
        "testsRoot": str(_TESTS_ROOT),
        "testsRootExists": tests_root.is_dir(),
        "components": components,
    }


def _resolve_repo_root(repo_root: str | Path | None) -> Path:
    if repo_root is None:
        return Path(__file__).resolve().parents[2]
    return Path(repo_root).expanduser().resolve()


def _reference_fixture_dirs(reference_root: Path) -> list[Path]:
    if not reference_root.is_dir():
        return []
    return sorted(
        path
        for path in reference_root.iterdir()
        if path.is_dir() and (path / "manifest.json").is_file()
    )


def _test_files(tests_root: Path) -> set[str]:
    if not tests_root.is_dir():
        return set()
    return {
        path.name
        for path in tests_root.rglob("test*.py")
        if path.is_file()
    }


def _coverage_row(
    *,
    capability: dict[str, Any],
    suite_exists: bool,
    fixture_dirs_exist: bool,
    test_files: set[str],
) -> dict[str, Any]:
    component = str(capability["component"])
    migration_status = str(capability["status"])
    matlab_reference_fixtures_exist = _matlab_reference_fixtures_exist(
        component=component,
        suite_exists=suite_exists,
        fixture_dirs_exist=fixture_dirs_exist,
    )
    parity_tests_exist = _tests_exist(
        test_files,
        _PARITY_TEST_FILES if component in _PARITY_COMPONENTS else set(),
    )
    distributional_tests_exist = _tests_exist(
        test_files,
        (
            _DISTRIBUTIONAL_TEST_FILES
            if component in _DISTRIBUTIONAL_COMPONENTS
            else set()
        ),
    )
    python_only_tests_exist = _tests_exist(
        test_files,
        _PYTHON_ONLY_TEST_FILES_BY_COMPONENT.get(component, set()),
    )

    coverage_status = _coverage_status(
        component=component,
        migration_status=migration_status,
        matlab_reference_fixtures_exist=matlab_reference_fixtures_exist,
        parity_tests_exist=parity_tests_exist,
        distributional_tests_exist=distributional_tests_exist,
        python_only_tests_exist=python_only_tests_exist,
    )

    return {
        "component": component,
        "label": capability["label"],
        "migrationStatus": migration_status,
        "pythonEntryPoint": capability["pythonEntryPoint"],
        "matlabReference": capability["matlabReference"],
        "matlabReferenceFixturesExist": matlab_reference_fixtures_exist,
        "parityTestsExist": parity_tests_exist,
        "distributionalTestsExist": distributional_tests_exist,
        "pythonOnlyTestsExist": python_only_tests_exist,
        "coverageStatus": coverage_status,
        "messages": _messages(
            migration_status=migration_status,
            coverage_status=coverage_status,
            matlab_reference_fixtures_exist=matlab_reference_fixtures_exist,
            parity_tests_exist=parity_tests_exist,
            distributional_tests_exist=distributional_tests_exist,
            python_only_tests_exist=python_only_tests_exist,
        ),
    }


def _matlab_reference_fixtures_exist(
    *,
    component: str,
    suite_exists: bool,
    fixture_dirs_exist: bool,
) -> bool:
    if component in _REFERENCE_COMPONENTS or component in _PARITY_COMPONENTS:
        return suite_exists and fixture_dirs_exist
    return fixture_dirs_exist


def _tests_exist(test_files: set[str], expected_files: set[str]) -> bool:
    return bool(expected_files) and bool(test_files.intersection(expected_files))


def _coverage_status(
    *,
    component: str,
    migration_status: str,
    matlab_reference_fixtures_exist: bool,
    parity_tests_exist: bool,
    distributional_tests_exist: bool,
    python_only_tests_exist: bool,
) -> str:
    if migration_status in {
        DOCUMENTATION_ONLY,
        MATLAB_REFERENCE,
        REFERENCE_ONLY,
        UNSUPPORTED,
    }:
        return migration_status

    if component in _NEVER_COMPLETE_COMPONENTS:
        return PARTIAL

    if migration_status == COMPLETE:
        return COMPLETE if python_only_tests_exist else PARTIAL

    if component in _CORE_REFERENCE_ELIGIBLE_PARTIAL_COMPONENTS:
        has_reference_test = parity_tests_exist or distributional_tests_exist
        if (
            matlab_reference_fixtures_exist
            and has_reference_test
            and python_only_tests_exist
        ):
            return COMPLETE

    return PARTIAL


def _messages(
    *,
    migration_status: str,
    coverage_status: str,
    matlab_reference_fixtures_exist: bool,
    parity_tests_exist: bool,
    distributional_tests_exist: bool,
    python_only_tests_exist: bool,
) -> list[str]:
    messages: list[str] = []

    if migration_status in {
        DOCUMENTATION_ONLY,
        MATLAB_REFERENCE,
        REFERENCE_ONLY,
        UNSUPPORTED,
    }:
        messages.append(
            f"Coverage status mirrors the capability status: {migration_status}."
        )

    if matlab_reference_fixtures_exist:
        messages.append("MATLAB reference suite and manifest fixtures were detected.")
    if parity_tests_exist:
        messages.append("Relevant parity test file was detected under tests/.")
    if distributional_tests_exist:
        messages.append(
            "Relevant distributional validation test file was detected under tests/."
        )
    if python_only_tests_exist:
        messages.append("Relevant Python-only test file was detected under tests/.")

    if coverage_status == PARTIAL:
        messages.append("Coverage is kept partial under conservative classification.")
    elif coverage_status == COMPLETE:
        messages.append("Static coverage signals support complete reference coverage.")

    return messages


__all__ = ["CONTRACT_VERSION", "get_reference_coverage"]
