function [appState, validationDisplay, validationError] = app_validate_v9(appState, uiState)
%APP_VALIDATE_V9 Refresh validation and confirm fatal MATLAB validation.

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

[appState, validationDisplay] = app_input_changed_v9(appState, uiState);
validationError = '';

try
    appState.CurrentConfig = validate_config_v9(appState.CurrentConfig);
    appState.LastValidationReport = collect_validation_issues_v9(appState.CurrentConfig);
    validationDisplay = validation_report_to_display_v9(appState.LastValidationReport);
catch ME
    validationError = ME.message;
end
end
