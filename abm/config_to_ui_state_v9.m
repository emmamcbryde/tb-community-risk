function uiState = config_to_ui_state_v9(config, schema)
%CONFIG_TO_UI_STATE_V9 Convert canonical config to schema-driven UI state.

if nargin < 2 || isempty(schema)
    schema = build_ui_schema_v9();
end

uiState = struct();
uiState.fieldValues = struct();

fields = [schema.simpleFields(:); schema.advancedFields(:)];
for i = 1:numel(fields)
    field = fields(i);
    key = schema_field_key_v9(field);
    uiState.fieldValues.(key) = get_nested_field_v9(config, field.sourceConfigField);
end
end

function value = get_nested_field_v9(s, path)
value = [];
if isempty(path)
    return;
end
parts = strsplit(path, '.');
value = s;
for i = 1:numel(parts)
    if ~isstruct(value) || ~isfield(value, parts{i})
        value = [];
        return;
    end
    value = value.(parts{i});
end
end

