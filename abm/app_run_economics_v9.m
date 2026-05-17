function [appState, resultsDisplay, errMsg] = app_run_economics_v9(appState, economicsUiState)
%APP_RUN_ECONOMICS_V9 Run economics post-processing for fresh core results.

errMsg = '';
if ~isfield(appState, 'CurrentEconomicsConfig') || isempty(appState.CurrentEconomicsConfig)
    appState.CurrentEconomicsConfig = build_default_economics_config_v9();
end

if nargin >= 2 && ~isempty(economicsUiState)
    [appState, ~] = app_economics_input_changed_v9(appState, economicsUiState);
end

[tf, message] = app_has_fresh_results_v9(appState);
if ~tf
    errMsg = message;
    resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
    return;
end

validationReport = validate_economics_config_v9(appState.CurrentEconomicsConfig);
appState.LastEconomicsValidationReport = validationReport;
if ~validationReport.isValid
    errMsg = 'Economics run blocked: fix fatal economics validation errors first.';
    resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
    return;
end

try
    appState.LastEconomics = run_health_economics_v9(appState.LastResults, ...
        appState.CurrentEconomicsConfig);
    appState.LastEconomicsConfig = appState.CurrentEconomicsConfig;
    appState.EconomicsAreStale = false;
    appState = app_refresh_bundle_v9(appState);
catch ME
    appState.EconomicsAreStale = true;
    errMsg = ME.message;
end

resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
end
