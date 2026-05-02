function [appState, saveInfo, errMsg] = app_save_scenario_v9(appState, filename)
%APP_SAVE_SCENARIO_V9 Save the current config as a scenario JSON file.

saveInfo = struct();
errMsg = '';

try
    report = collect_validation_issues_v9(appState.CurrentConfig);
    appState.LastValidationReport = report;
    if isfield(report, 'errors') && ~isempty(report.errors)
        errMsg = 'Save blocked: fix fatal validation errors before saving this scenario.';
        return;
    end

    saveInfo = save_scenario_v9(appState.CurrentConfig, filename);
    appState.LastScenarioFile = filename;
catch ME
    errMsg = ME.message;
end
end
