function [appState, saveInfo, errMsg] = app_save_scenario_v9(appState, filename)
%APP_SAVE_SCENARIO_V9 Save the current config as a scenario JSON file.

saveInfo = struct();
errMsg = '';

try
    saveInfo = save_scenario_v9(appState.CurrentConfig, filename);
    appState.LastScenarioFile = filename;
catch ME
    errMsg = ME.message;
end
end
