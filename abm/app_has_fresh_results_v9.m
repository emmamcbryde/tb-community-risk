function [tf, message] = app_has_fresh_results_v9(appState)
%APP_HAS_FRESH_RESULTS_V9 Return whether fresh run results are available.

tf = false;
message = 'No run results are available.';

if ~isstruct(appState) || ~isfield(appState, 'LastResults') || isempty(appState.LastResults)
    return;
end
if ~isfield(appState, 'ResultsAreStale') || logical(appState.ResultsAreStale)
    message = 'Results are stale. Please rerun the scenario before using add-ons or export.';
    return;
end
if ~isfield(appState, 'LastRunSucceeded') || ~logical(appState.LastRunSucceeded)
    message = 'The last run did not complete successfully.';
    return;
end

tf = true;
message = '';
end
