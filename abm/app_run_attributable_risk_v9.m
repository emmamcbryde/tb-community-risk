function [appState, resultsDisplay, errMsg] = app_run_attributable_risk_v9(appState, Nattr)
%APP_RUN_ATTRIBUTABLE_RISK_V9 Run attributable-risk add-on for current results.

if nargin < 2 || isempty(Nattr)
    Nattr = 50000;
end

errMsg = '';
[tf, message] = app_has_fresh_results_v9(appState);
if ~tf
    errMsg = message;
    resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
    return;
end

try
    appState.LastAttributableRisk = run_tb_screening_reactivation_attributable_v9(appState.LastResults, Nattr);
    appState = app_refresh_bundle_v9(appState);
catch ME
    errMsg = ME.message;
end

resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
end
