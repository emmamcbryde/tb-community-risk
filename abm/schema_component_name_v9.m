function componentName = schema_component_name_v9(fieldOrPath, inputType)
%SCHEMA_COMPONENT_NAME_V9 Return the canonical App Designer component name.

if nargin < 2 || isempty(inputType)
    inputType = '';
end

source = '';
if isstruct(fieldOrPath)
    if isfield(fieldOrPath, 'sourceConfigField') && ~isempty(fieldOrPath.sourceConfigField)
        source = fieldOrPath.sourceConfigField;
    elseif isfield(fieldOrPath, 'internalName') && ~isempty(fieldOrPath.internalName)
        source = fieldOrPath.internalName;
    end
    if isempty(inputType) && isfield(fieldOrPath, 'inputType')
        inputType = fieldOrPath.inputType;
    end
elseif ischar(fieldOrPath) || (isstring(fieldOrPath) && isscalar(fieldOrPath))
    source = char(string(fieldOrPath));
end

parts = regexp(source, '[a-zA-Z0-9]+', 'match');
if isempty(parts)
    baseName = 'Field';
else
    for i = 1:numel(parts)
        token = parts{i};
        parts{i} = [upper(token(1)) token(2:end)];
    end
    baseName = strjoin(parts, '');
end

switch lower(char(string(inputType)))
    case 'select'
        suffix = 'DropDown';
    case 'table'
        suffix = 'Table';
    case 'boolean'
        suffix = 'CheckBox';
    otherwise
        suffix = 'EditField';
end

if endsWith(baseName, suffix, 'IgnoreCase', true)
    componentName = baseName;
else
    componentName = [baseName suffix];
end
end
