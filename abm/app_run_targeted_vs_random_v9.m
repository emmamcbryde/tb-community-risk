function [appState, resultsDisplay, errMsg] = app_run_targeted_vs_random_v9(appState, options)
%APP_RUN_TARGETED_VS_RANDOM_V9 Run explicit comparison add-ons for current results.

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
testType = get_option(options, 'testType', cfg.testType);
nReps = get_option(options, 'nReps', cfg.nReps);
screeningStrategy = get_option(options, 'screeningStrategy', cfg.screeningStrategy);
screenCoverage = get_option(options, 'screenCoverage', cfg.screenCoverage);
coverageGrid = get_option(options, 'coverageGrid', [0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.80 1.00]);
Nprofile = get_option(options, 'Nprofile', 50000);
seed = get_option(options, 'seed', cfg.seed);

try
    out = struct();
    out.compare = run_tb_screening_compare_strategies_v9(screeningStrategy, screenCoverage, nReps);
    out.gradient = run_tb_screening_targeted_gradient_v9(testType, regimen, nReps);
    out.profile = run_tb_screening_targeting_profile_v9(regimen, coverageGrid, Nprofile, seed);
    out.optima = run_tb_screening_targeting_optima_v9(regimen, nReps, coverageGrid, Nprofile, seed);
    appState.LastTargetedVsRandom = out;
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
