function uiState = app_collect_ui_state_v9(app, schema)
%APP_COLLECT_UI_STATE_V9 Collect App Designer control values into uiState.

if nargin < 2 || isempty(schema)
    schema = build_ui_schema_v9();
end

uiState = struct();
uiState.fieldValues = struct();

fields = [schema.simpleFields(:); schema.advancedFields(:)];
for i = 1:numel(fields)
    field = fields(i);
    component = find_component_v9(app, field);
    if isempty(component)
        continue;
    end

    value = read_component_value_v9(component);
    if isempty(value) && ~has_value_property_v9(component)
        continue;
    end

    key = schema_field_key_v9(field);
    uiState.fieldValues.(key) = value;
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

function tf = has_value_property_v9(component)
tf = isprop(component, 'Value') || isprop(component, 'Data');
end

function value = read_component_value_v9(component)
value = [];
if isprop(component, 'Value')
    value = component.Value;
elseif isprop(component, 'Data')
    value = component.Data;
end
end
