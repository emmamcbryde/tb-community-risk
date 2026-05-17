function [appState, validationDisplay, resultsDisplay, loadInfo] = app_load_scenario_v9(appState, filename)
%APP_LOAD_SCENARIO_V9 Load a scenario JSON and mark prior results stale.

[cfg, report, loadInfo] = load_scenario_v9(filename);

appState.CurrentConfig = cfg;
appState.CurrentUIState = config_to_ui_state_v9(cfg, appState.Schema);
appState.CurrentEconomicsConfig = get_loaded_economics_config(loadInfo);
appState.LastEconomicsConfig = [];
appState.LastEconomics = [];
appState.EconomicsAreStale = true;
appState.LastEconomicsPreset = '';
appState.LastEconomicsValidationReport = validate_economics_config_v9(appState.CurrentEconomicsConfig);
appState.LastValidationReport = report;
appState.LastScenarioFile = filename;
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

function econConfig = get_loaded_economics_config(loadInfo)
econConfig = build_default_economics_config_v9();
if isstruct(loadInfo) && isfield(loadInfo, 'economicsConfig') && ~isempty(loadInfo.economicsConfig)
    econConfig = loadInfo.economicsConfig;
end
end
