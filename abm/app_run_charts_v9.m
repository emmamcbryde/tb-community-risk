function [appState, resultsDisplay, errMsg] = app_run_charts_v9(appState, options)
%APP_RUN_CHARTS_V9 Run explicit chart add-on for current results.

if nargin < 2 || isempty(options)
    options = struct();
end

errMsg = '';
[tf, message] = app_has_fresh_results_v9(appState);
if ~tf
    errMsg = message;
    resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
    return;
end

cfg = appState.LastRunConfig;
if isempty(cfg)
    cfg = appState.CurrentConfig;
end

regimen = get_option(options, 'regimen', cfg.regimen);
nReps = get_option(options, 'nReps', cfg.nReps);
coverageGrid = get_option(options, 'coverageGrid', [0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.80 1.00]);

try
    appState.LastCharts = run_tb_screening_igra_charts_v9(regimen, nReps, coverageGrid);
    appState = app_refresh_bundle_v9(appState);
catch ME
    errMsg = ME.message;
end

resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
end

function value = get_option(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
