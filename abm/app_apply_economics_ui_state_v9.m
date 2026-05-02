function app = app_apply_economics_ui_state_v9(app, uiState)
%APP_APPLY_ECONOMICS_UI_STATE_V9 Push economics UI state into controls.

if ~isstruct(uiState) || ~isfield(uiState, 'fieldValues') || ~isstruct(uiState.fieldValues)
    return;
end

names = fieldnames(uiState.fieldValues);
for i = 1:numel(names)
    name = names{i};
    if ~isprop(app, name)
        continue;
    end
    write_component_value(app.(name), uiState.fieldValues.(name));
end
end

function write_component_value(component, value)
if isempty(value)
    if isprop(component, 'AllowEmpty') && isprop(component, 'Value')
        component.AllowEmpty = 'on';
        component.Value = [];
    elseif isprop(component, 'Value')
        component.Value = '';
    end
    return;
end

if isprop(component, 'Value')
    if isprop(component, 'Items')
        [hasMatch, matchedValue] = match_dropdown_item(component.Items, value);
        if ~hasMatch
            return;
        end
        component.Value = matchedValue;
        return;
    end
    component.Value = value;
elseif isprop(component, 'Data')
    component.Data = value;
end
end

function [hasMatch, matchedValue] = match_dropdown_item(items, value)
hasMatch = false;
matchedValue = [];
valueText = string(value);
if ~isscalar(valueText)
    return;
end

itemTexts = string(items);
idx = find(itemTexts == valueText, 1, 'first');
if isempty(idx)
    idx = find(strcmpi(cellstr(itemTexts(:)), char(valueText)), 1, 'first');
end
if isempty(idx)
    return;
end

hasMatch = true;
if iscell(items)
    matchedValue = items{idx};
else
    matchedValue = items(idx);
end
end
