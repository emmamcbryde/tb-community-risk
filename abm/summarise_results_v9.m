function summary = summarise_results_v9(results)
%SUMMARISE_RESULTS_V9 Build a compact app-facing summary from v9 results.

if ~isstruct(results) || ~isfield(results, 'summary') || ~istable(results.summary)
    error('Input must be a v9 results struct containing results.summary.');
end

keepMetrics = [ ...
    "nScreened", "nTestPositiveNonActive", "nFalsePositiveTreated", ...
    "nTotalCoursesStarted", "nTotalCoursesCompleted", ...
    "nCuredInfection", "nPreventedActiveTB", ...
    "NNS_cureInfection", "NNS_preventActiveTB", ...
    "NNT_started_cureInfection", "NNT_started_preventActiveTB" ];

summary = struct();
summary.strategy = results.strategy;
summary.calibration = results.calibration;
summary.keyMetrics = results.summary(ismember(results.summary.Metric, keepMetrics), :);
summary.summaryTable = results.summary;
end
