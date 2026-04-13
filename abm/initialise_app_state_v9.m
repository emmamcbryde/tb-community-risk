function appState = initialise_app_state_v9()
%INITIALISE_APP_STATE_V9 Create the default APY v9 app state container.

schema = build_ui_schema_v9();
config = build_default_config_v9();
report = collect_validation_issues_v9(config);

appState = struct();
appState.Schema = schema;
appState.CurrentConfig = config;
appState.CurrentUIState = config_to_ui_state_v9(config, schema);
appState.LastValidationReport = report;
appState.LastResults = [];
appState.LastBundle = [];
appState.LastExports = [];
appState.LastDoNothing = [];
appState.LastAttributableRisk = [];
appState.LastTargetedVsRandom = [];
appState.LastCharts = [];
appState.LastRunConfig = [];
appState.ResultsAreStale = true;
appState.ValidationNeedsRefresh = false;
appState.LastScenarioFile = '';
appState.LastRunSucceeded = false;
end
