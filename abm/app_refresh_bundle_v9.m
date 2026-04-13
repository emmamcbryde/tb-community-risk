function appState = app_refresh_bundle_v9(appState)
%APP_REFRESH_BUNDLE_V9 Rebuild the current results bundle from appState.

if ~isstruct(appState) || ~isfield(appState, 'LastResults') || isempty(appState.LastResults)
    appState.LastBundle = [];
    return;
end

args = {};
args = add_if_present(args, 'validationReport', get_if_present(appState, 'LastValidationReport'));
args = add_if_present(args, 'targetedVsRandom', get_if_present(appState, 'LastTargetedVsRandom'));
args = add_if_present(args, 'doNothing', get_if_present(appState, 'LastDoNothing'));
args = add_if_present(args, 'attributableRisk', get_if_present(appState, 'LastAttributableRisk'));
args = add_if_present(args, 'charts', get_if_present(appState, 'LastCharts'));
args = add_if_present(args, 'downloads', get_if_present(appState, 'LastExports'));

appState.LastBundle = build_results_bundle_v9(appState.LastResults, args{:});
end

function args = add_if_present(args, name, value)
if isempty(value)
    return;
end
args = [args, {name, value}]; %#ok<AGROW>
end

function value = get_if_present(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = [];
end
end
