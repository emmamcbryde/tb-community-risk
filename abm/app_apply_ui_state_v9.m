function app = app_apply_ui_state_v9(app, uiState, schema)
%APP_APPLY_UI_STATE_V9 Push uiState values into App Designer controls.

if nargin < 3 || isempty(schema)
    schema = build_ui_schema_v9();
end
if ~isstruct(uiState) || ~isfield(uiState, 'fieldValues') || ~isstruct(uiState.fieldValues)
    return;
end

fields = [schema.simpleFields(:); schema.advancedFields(:)];
for i = 1:numel(fields)
    field = fields(i);
    key = schema_field_key_v9(field);
    if ~isfield(uiState.fieldValues, key)
        continue;
    end

    component = find_component_v9(app, field);
    if isempty(component)
        continue;
    end

    write_component_value_v9(component, uiState.fieldValues.(key));
end
end

function component = find_component_v9(app, field)
component = [];
candidateNames = { ...
    schema_component_name_v9(field), ...
    schema_field_key_v9(field)};

for i = 1:numel(candidateNames)
    name = candidateNames{i};
    if isprop(app, name)
        component = app.(name);
        return;
    end
end
end

function write_component_value_v9(component, value)
if isempty(value)
    if isprop(component, 'AllowEmpty') && isprop(component, 'Value')
        component.AllowEmpty = 'on';
        component.Value = [];
    end
    return;
end

if isprop(component, 'Value')
    if isprop(component, 'Items')
        [hasMatch, matchedValue] = match_dropdown_item_v9(component.Items, value);
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

function [hasMatch, matchedValue] = match_dropdown_item_v9(items, value)
hasMatch = false;
matchedValue = [];

valueText = scalar_text_v9(value);
if isempty(valueText)
    return;
end

itemTexts = string(items);
exactMatch = (itemTexts == valueText);
if any(exactMatch, 'all')
    hasMatch = true;
    matchedValue = get_item_value_v9(items, find(exactMatch, 1, 'first'));
    return;
end

caseInsensitiveMatch = strcmpi(cellstr(itemTexts(:)), char(valueText));
if any(caseInsensitiveMatch)
    hasMatch = true;
    matchedValue = get_item_value_v9(items, find(caseInsensitiveMatch, 1, 'first'));
end
end

function textValue = scalar_text_v9(value)
textValue = string.empty;

if iscell(value)
    if ~isscalar(value)
        return;
    end
    value = value{1};
end

try
    textValue = string(value);
catch
    textValue = string.empty;
    return;
end

if ~isscalar(textValue)
    textValue = string.empty;
end
end

function itemValue = get_item_value_v9(items, index)
if iscell(items)
    itemValue = items{index};
elseif ischar(items) && size(items, 1) > 1
    itemValue = deblank(items(index, :));
else
    itemValue = items(index);
end
end
