function [appState, resultsDisplay, errMsg] = app_export_v9(appState, baseName)
%APP_EXPORT_V9 Export current fresh results and refresh the downloads bundle section.

if nargin < 2 || isempty(baseName)
    baseName = [];
end

errMsg = '';
[tf, message] = app_has_fresh_results_v9(appState);
if ~tf
    errMsg = message;
    resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
    return;
end

if isempty(baseName)
    baseName = default_export_base_name(appState.LastRunConfig);
end

try
    appState.LastExports = export_results_v9(appState.LastResults, baseName);
    appState = app_refresh_bundle_v9(appState);
catch ME
    errMsg = ME.message;
end

resultsDisplay = results_bundle_to_display_v9(appState.LastBundle);
end

function baseName = default_export_base_name(config)
baseName = 'tb_screening_v9';
if isstruct(config) && isfield(config, 'scenarioLabel') && ~isempty(config.scenarioLabel)
    baseName = regexprep(char(string(config.scenarioLabel)), '[^a-zA-Z0-9_-]', '_');
end
end
