function tf = configs_substantively_equal_v9(a, b)
%CONFIGS_SUBSTANTIVELY_EQUAL_V9 Compare configs ignoring UI metadata fields.

ignoreFields = {'scenarioLabel', 'useDefaults', 'outputDir'};

a = strip_ignored_fields(a, ignoreFields);
b = strip_ignored_fields(b, ignoreFields);

tf = isequaln(a, b);
end

function s = strip_ignored_fields(s, ignoreFields)
if ~isstruct(s)
    return;
end

for i = 1:numel(ignoreFields)
    fieldName = ignoreFields{i};
    if isfield(s, fieldName)
        s = rmfield(s, fieldName);
    end
end
end
