from __future__ import annotations

from copy import deepcopy
from typing import Any


CONTRACT_VERSION = "apy_python_migration_capabilities_v1"

COMPLETE = "complete"
PARTIAL = "partial"
REFERENCE_ONLY = "reference_only"
DOCUMENTATION_ONLY = "documentation_only"
UNSUPPORTED = "unsupported"
MATLAB_REFERENCE = "matlab_reference"
UNKNOWN = "unknown"

_COMPONENTS: tuple[dict[str, Any], ...] = (
    {
        "component": "matlab_reference_engine",
        "label": "MATLAB APY v9 reference engine",
        "status": MATLAB_REFERENCE,
        "pythonEntryPoint": None,
        "matlabReference": "abm/tb_screening_mc_model_v9.m",
        "notes": "Authoritative APY v9 reference engine until parity criteria are updated.",
    },
    {
        "component": "default_config",
        "label": "Default APY v9 scenario config",
        "status": COMPLETE,
        "pythonEntryPoint": "engine.apy.config.build_default_config",
        "matlabReference": "build_default_config_v9.m",
        "notes": "Builds the APY v9 default config using portable repo paths.",
    },
    {
        "component": "config_validation",
        "label": "Scenario config validation",
        "status": COMPLETE,
        "pythonEntryPoint": "engine.apy.validation.collect_validation_issues",
        "matlabReference": "collect_validation_issues_v9.m",
        "notes": "Validates and reports model-blocking scenario configuration issues.",
    },
    {
        "component": "data_loading",
        "label": "APY input data loading",
        "status": COMPLETE,
        "pythonEntryPoint": "engine.apy.data",
        "matlabReference": "default_data.csv",
        "notes": "Loads APY source tables through portable Python paths.",
    },
    {
        "component": "age_distribution",
        "label": "Age distribution loading and sampling",
        "status": COMPLETE,
        "pythonEntryPoint": "engine.apy.age_distribution",
        "matlabReference": "default_age_distribution.csv",
        "notes": "Provides the Python age-distribution helpers used by APY runs.",
    },
    {
        "component": "regimen_handling",
        "label": "Regimen configuration and partial-course rules",
        "status": COMPLETE,
        "pythonEntryPoint": "engine.apy.regimen",
        "matlabReference": "tb_screening_mc_model_v9.m",
        "notes": "Ports APY regimen defaults, validation, and partial-efficacy rules.",
    },
    {
        "component": "calibration",
        "label": "APY calibration helpers",
        "status": COMPLETE,
        "pythonEntryPoint": "engine.apy.calibration",
        "matlabReference": "tb_screening_mc_model_v9.m",
        "notes": "Provides the tested calibration helpers used by the Python simulation.",
    },
    {
        "component": "cohort_primitives",
        "label": "Single-cohort simulation primitives",
        "status": COMPLETE,
        "pythonEntryPoint": "engine.apy.cohort",
        "matlabReference": "tb_screening_mc_model_v9.m",
        "notes": "Provides tested Python cohort primitives for APY simulation.",
    },
    {
        "component": "core_simulation",
        "label": "Core stochastic APY simulation",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.runner.run_scenario",
        "matlabReference": "tb_screening_mc_model_v9.m",
        "notes": (
            "Runs the Python APY port; MATLAB remains the reference engine while "
            "distributional parity is expanded."
        ),
    },
    {
        "component": "scenario_runner",
        "label": "Python APY scenario runner",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.runner.run_scenario",
        "matlabReference": "run_scenario_v9.m",
        "notes": "Runs Python scenarios and preserves the current backend-facing bundle shape.",
    },
    {
        "component": "scenario_save_load",
        "label": "Scenario save/load",
        "status": PARTIAL,
        "pythonEntryPoint": "adapters.python_apy_backend.PythonApyBackend",
        "matlabReference": "save_scenario_v9.m; load_scenario_v9.m",
        "notes": "Supports the current Python backend scenario JSON workflow without claiming full MATLAB UI parity.",
    },
    {
        "component": "summary_rows",
        "label": "APY summary rows",
        "status": COMPLETE,
        "pythonEntryPoint": "engine.apy.summary",
        "matlabReference": "summarise_results_v9.m",
        "notes": "Builds stable summary rows used by the Python result bundle.",
    },
    {
        "component": "results_bundle",
        "label": "Results bundle and headline summaries",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.results_bundle.build_results_bundle",
        "matlabReference": "build_results_bundle_v9.m",
        "notes": "Provides the Python result bundle consumed by current app flows.",
    },
    {
        "component": "exports",
        "label": "JSON, CSV, table, and chart export helpers",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.exports",
        "matlabReference": "export_results_v9.m",
        "notes": "Covers current Streamlit download/display helpers, not full MATLAB export parity.",
    },
    {
        "component": "export_display_helpers",
        "label": "Export and display helper parity",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.exports",
        "matlabReference": "export_results_v9.m",
        "notes": "Implemented for current JSON, CSV, and display-table needs without full MATLAB helper parity.",
    },
    {
        "component": "chart_rendering_parity",
        "label": "Chart rendering parity",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.exports.chart_numeric_series",
        "matlabReference": "plot_results_v9.m",
        "notes": "Provides chart-ready numeric series for current app flows; MATLAB chart rendering parity is not complete.",
    },
    {
        "component": "natural_history",
        "label": "Do-nothing and natural-history add-ons",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.natural_history",
        "matlabReference": "run_tb_screening_do_nothing_v9.m",
        "notes": "Supports do-nothing summaries and dynamic comparison inputs.",
    },
    {
        "component": "do_nothing_dynamic_comparison",
        "label": "Do-nothing dynamic comparison",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.runner.run_scenario_with_do_nothing",
        "matlabReference": "run_tb_screening_do_nothing_v9.m",
        "notes": "Builds Python do-nothing comparison payloads for current bundle flows; full MATLAB add-on parity is not claimed.",
    },
    {
        "component": "natural_history_addon_reporting",
        "label": "Natural-history add-on reporting",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.natural_history.build_natural_history_addon_report",
        "matlabReference": "run_tb_screening_do_nothing_v9.m",
        "notes": "Reports current Python natural-history add-ons without claiming full MATLAB add-on parity.",
    },
    {
        "component": "targeting_compare",
        "label": "Targeting result-bundle comparison",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.targeting.compare_targeting_result_bundles",
        "matlabReference": "run_tb_screening_compare_strategies_v9.m",
        "notes": "Compares already-run bundles; it does not run full MATLAB targeting sweeps.",
    },
    {
        "component": "targeting_strategy_comparison",
        "label": "Targeting and strategy comparison",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.targeting",
        "matlabReference": "run_tb_screening_compare_strategies_v9.m",
        "notes": "Covers implemented bundle comparison helpers, not full MATLAB targeting-profile execution.",
    },
    {
        "component": "economics",
        "label": "Health economics",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.economics.calculate_economics",
        "matlabReference": "run_health_economics_v9.m",
        "notes": (
            "Calculates the tested partial economics subset and reports unsupported "
            "economics components explicitly."
        ),
    },
    {
        "component": "health_economics",
        "label": "Health economics workflow",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.economics.calculate_economics",
        "matlabReference": "run_health_economics_v9.m",
        "notes": "Explicit health-economics scope row for the tested partial Python economics implementation.",
    },
    {
        "component": "attributable_risk",
        "label": "Reactivation attributable-risk add-on",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.attributable_risk.run_attributable_risk",
        "matlabReference": "run_tb_screening_reactivation_attributable_v9.m",
        "notes": "Python scaffold validates inputs and passes through precomputed rows; full MATLAB metric parity is not ported.",
    },
    {
        "component": "attributable_risk_analysis",
        "label": "Attributable-risk analysis",
        "status": PARTIAL,
        "pythonEntryPoint": "engine.apy.attributable_risk.run_attributable_risk",
        "matlabReference": "run_tb_screening_reactivation_attributable_v9.m",
        "notes": "Implemented scaffold and reporting contract exist, but reactivation attributable-risk calculations are not full parity.",
    },
    {
        "component": "matlab_reference_fixtures",
        "label": "MATLAB reference fixture suite",
        "status": REFERENCE_ONLY,
        "pythonEntryPoint": "engine.apy.reference_loader",
        "matlabReference": "validation/matlab_reference",
        "notes": "Reference fixtures support parity checks and are not ordinary model outputs.",
    },
    {
        "component": "distributional_validation",
        "label": "Distributional validation runner",
        "status": REFERENCE_ONLY,
        "pythonEntryPoint": "engine.apy.distributional_validation",
        "matlabReference": "export_matlab_reference_scenarios_v9.m",
        "notes": "Diagnostic validation against MATLAB reference fixtures; it is not a replacement model workflow.",
    },
    {
        "component": "matlab_user_options_compatibility",
        "label": "MATLAB user-options compatibility",
        "status": REFERENCE_ONLY,
        "pythonEntryPoint": None,
        "matlabReference": "run_tb_screening_user_options_v9.m",
        "notes": "Compatibility target for Python APY config fields and workflow expectations.",
    },
    {
        "component": "matlab_app_callbacks",
        "label": "MATLAB App Designer callbacks",
        "status": DOCUMENTATION_ONLY,
        "pythonEntryPoint": None,
        "matlabReference": "app_*_v9.m",
        "notes": "Callbacks are UI wrappers and are not part of the Python APY engine port.",
    },
    {
        "component": "matlab_app_designer_helpers",
        "label": "MATLAB App Designer helpers",
        "status": DOCUMENTATION_ONLY,
        "pythonEntryPoint": None,
        "matlabReference": "app_*_v9.m",
        "notes": "Documented as MATLAB UI helper behavior only; no Python engine parity is claimed.",
    },
)

_CAPABILITY_CONTRACT: dict[str, Any] = {
    "contractVersion": CONTRACT_VERSION,
    "backend": "python_apy",
    "modelVersion": "python_apy_v9_port",
    "referenceEngine": "abm/tb_screening_mc_model_v9.m",
    "matlabRequired": False,
    "overallStatus": PARTIAL,
    "principle": "Planning and sequencing support; not a basis for denying care.",
    "statusValues": [
        COMPLETE,
        PARTIAL,
        REFERENCE_ONLY,
        DOCUMENTATION_ONLY,
        UNSUPPORTED,
        MATLAB_REFERENCE,
    ],
    "components": list(_COMPONENTS),
}


def get_apy_capabilities() -> dict[str, Any]:
    """Return the JSON-serialisable APY migration capability contract."""
    return deepcopy(_CAPABILITY_CONTRACT)


def get_component_status(component: str) -> str:
    """Return the migration status for one component, or ``unknown``."""
    for row in _COMPONENTS:
        if row["component"] == component:
            return str(row["status"])
    return UNKNOWN


def capability_rows() -> list[dict[str, Any]]:
    """Return flat component capability rows for display or tests."""
    return deepcopy(list(_COMPONENTS))


def unsupported_or_partial_components() -> list[dict[str, Any]]:
    """Return rows whose current Python migration status is partial or unsupported."""
    return [
        row
        for row in capability_rows()
        if row["status"] in {PARTIAL, UNSUPPORTED}
    ]
