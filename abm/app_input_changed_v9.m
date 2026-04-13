function [appState, validationDisplay, resultsDisplay] = app_input_changed_v9(appState, uiState)
%APP_INPUT_CHANGED_V9 Update config/validation and mark prior results stale.

cfg = ui_state_to_config_v9(uiState, appState.Schema, appState.CurrentConfig);
report = collect_validation_issues_v9(cfg);

appState.CurrentConfig = cfg;
appState.CurrentUIState = uiState;
appState.LastValidationReport = report;
appState.ResultsAreStale = true;
appState.LastRunSucceeded = false;
appState.LastExports = [];
appState.LastDoNothing = [];
appState.LastAttributableRisk = [];
appState.LastTargetedVsRandom = [];
appState.LastCharts = [];

validationDisplay = validation_report_to_display_v9(report);
resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
end
