function config = ui_state_to_config_v9(uiState, arg2, arg3)
%UI_STATE_TO_CONFIG_V9 Convert schema-driven UI state into canonical config.

schema = [];
baseConfig = [];

if nargin >= 2 && ~isempty(arg2)
    if is_schema_struct(arg2)
        schema = arg2;
    else
        baseConfig = arg2;
    end
end

if nargin >= 3 && ~isempty(arg3)
    baseConfig = arg3;
end

if isempty(schema)
    schema = build_ui_schema_v9();
end
if isempty(baseConfig)
    baseConfig = build_default_config_v9();
end

config = baseConfig;
if ~isstruct(uiState) || ~isfield(uiState, 'fieldValues') || ~isstruct(uiState.fieldValues)
    return;
end

fields = [schema.simpleFields(:); schema.advancedFields(:)];
for i = 1:numel(fields)
    field = fields(i);
    key = ui_state_key_v9(field);
    if isfield(uiState.fieldValues, key)
        config = set_nested_field_v9(config, field.sourceConfigField, uiState.fieldValues.(key));
    end
end
end

function s = set_nested_field_v9(s, path, value)
parts = strsplit(path, '.');
if numel(parts) == 1
    s.(parts{1}) = value;
    return;
end

head = parts{1};
tail = strjoin(parts(2:end), '.');
if ~isfield(s, head) || ~isstruct(s.(head))
    s.(head) = struct();
end
s.(head) = set_nested_field_v9(s.(head), tail, value);
end

function key = ui_state_key_v9(field)
source = field.sourceConfigField;
if isempty(source)
    source = field.internalName;
end
key = regexprep(source, '[^a-zA-Z0-9_]', '_');
if isempty(key)
    key = 'field';
end
if isstrprop(key(1), 'digit')
    key = ['f_' key];
end
end

function tf = is_schema_struct(value)
tf = isstruct(value) && isfield(value, 'simpleFields') && isfield(value, 'advancedFields');
end
