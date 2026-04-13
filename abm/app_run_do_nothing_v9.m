function [appState, resultsDisplay, errMsg] = app_run_do_nothing_v9(appState)
%APP_RUN_DO_NOTHING_V9 Run the do-nothing add-on for the current results.

errMsg = '';
[tf, message] = app_has_fresh_results_v9(appState);
if ~tf
    errMsg = message;
    resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
    return;
end

try
    appState.LastDoNothing = run_tb_screening_do_nothing_v9(appState.LastResults);
    appState = app_refresh_bundle_v9(appState);
catch ME
    errMsg = ME.message;
end

resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
end
