function [appState, validationReport] = app_economics_input_changed_v9(appState, economicsUiState)
%APP_ECONOMICS_INPUT_CHANGED_V9 Update economics config without touching core results.

if ~isfield(appState, 'CurrentEconomicsConfig') || isempty(appState.CurrentEconomicsConfig)
    appState.CurrentEconomicsConfig = build_default_economics_config_v9();
end

if nargin < 2 || isempty(economicsUiState)
    economicsUiState = economics_config_to_ui_state_v9(appState.CurrentEconomicsConfig);
end

appState.CurrentEconomicsConfig = ui_state_to_economics_config_v9(economicsUiState, ...
    appState.CurrentEconomicsConfig);
validationReport = validate_economics_config_v9(appState.CurrentEconomicsConfig);
appState.LastEconomicsValidationReport = validationReport;

if isfield(appState, 'LastEconomics') && ~isempty(appState.LastEconomics)
    appState.EconomicsAreStale = true;
end
end
