function key = schema_field_key_v9(fieldOrPath)
%SCHEMA_FIELD_KEY_V9 Return the stable uiState field key for a schema field.

source = '';
if isstruct(fieldOrPath)
    if isfield(fieldOrPath, 'sourceConfigField') && ~isempty(fieldOrPath.sourceConfigField)
        source = fieldOrPath.sourceConfigField;
    elseif isfield(fieldOrPath, 'internalName') && ~isempty(fieldOrPath.internalName)
        source = fieldOrPath.internalName;
    end
elseif ischar(fieldOrPath) || (isstring(fieldOrPath) && isscalar(fieldOrPath))
    source = char(string(fieldOrPath));
end

key = regexprep(source, '[^a-zA-Z0-9_]', '_');
if isempty(key)
    key = 'field';
end
if isstrprop(key(1), 'digit')
    key = ['f_' key];
end
end
