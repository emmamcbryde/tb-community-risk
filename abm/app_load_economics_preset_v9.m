function [appState, economicsUiState, validationReport, errMsg] = app_load_economics_preset_v9(appState, presetName)
%APP_LOAD_ECONOMICS_PRESET_V9 Load a named economics preset into app state.

if nargin < 2 || isempty(presetName)
    presetName = 'KWAB150';
end

errMsg = '';
try
    switch upper(strtrim(char(string(presetName))))
        case 'KWAB150'
            appState.CurrentEconomicsConfig = build_economics_preset_kwab150_v9();
            appState.LastEconomicsPreset = 'KWAB150';
        otherwise
            error('app_load_economics_preset_v9:UnknownPreset', ...
                'Unknown economics preset: %s', char(string(presetName)));
    end

    validationReport = validate_economics_config_v9(appState.CurrentEconomicsConfig);
    appState.LastEconomicsValidationReport = validationReport;
    if isfield(appState, 'LastEconomics') && ~isempty(appState.LastEconomics)
        appState.EconomicsAreStale = true;
    end
catch ME
    validationReport = get_if_present(appState, 'LastEconomicsValidationReport');
    errMsg = ME.message;
end

economicsUiState = economics_config_to_ui_state_v9(get_if_present_or_default(appState, ...
    'CurrentEconomicsConfig', build_default_economics_config_v9()));
end

function value = get_if_present(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = [];
end
end

function value = get_if_present_or_default(s, fieldName, defaultValue)
value = get_if_present(s, fieldName);
if isempty(value)
    value = defaultValue;
end
end
