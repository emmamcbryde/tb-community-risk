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
if isprop(component, 'Value')
    component.Value = value;
elseif isprop(component, 'Data')
    component.Data = value;
end
end
