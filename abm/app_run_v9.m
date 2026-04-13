function [appState, resultsDisplay, runError] = app_run_v9(appState, uiState)
%APP_RUN_V9 Validate config, run the scenario, and rebuild the bundle.

if nargin < 2 || isempty(uiState)
    if isstruct(appState) && isfield(appState, 'CurrentUIState') && ~isempty(appState.CurrentUIState)
        uiState = appState.CurrentUIState;
    else
        schema = [];
        if isstruct(appState) && isfield(appState, 'Schema') && ~isempty(appState.Schema)
            schema = appState.Schema;
        end
        uiState = config_to_ui_state_v9(appState.CurrentConfig, schema);
    end
end

runError = '';
[appState, ~] = app_input_changed_v9(appState, uiState);

try
    cfg = validate_config_v9(appState.CurrentConfig);
    report = collect_validation_issues_v9(cfg);
    results = run_scenario_v9(cfg);

    appState.CurrentConfig = cfg;
    appState.CurrentUIState = uiState;
    appState.LastValidationReport = report;
    appState.LastResults = results;
    appState.LastRunConfig = cfg;
    appState.ResultsAreStale = false;
    appState.LastRunSucceeded = true;
    appState.LastExports = [];
    appState.LastDoNothing = [];
    appState.LastAttributableRisk = [];
    appState.LastTargetedVsRandom = [];
    appState.LastCharts = [];
    appState = app_refresh_bundle_v9(appState);
catch ME
    appState.ResultsAreStale = true;
    appState.LastRunSucceeded = false;
    runError = ME.message;
end

resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
end
