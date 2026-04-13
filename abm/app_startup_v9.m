function [appState, validationDisplay, resultsDisplay] = app_startup_v9()
%APP_STARTUP_V9 Initialize controller state and display models.

appState = initialise_app_state_v9();
validationDisplay = validation_report_to_display_v9(appState.LastValidationReport);
resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
end
