function jsonText = run_health_economics_json_v9(config, econConfig)
%RUN_HEALTH_ECONOMICS_JSON_V9 Run v9 economics and return JSON text.
%
% The Python Streamlit layer keeps only JSON-like state, so this wrapper keeps
% MATLAB structs/tables inside MATLAB and returns encoded text.

config = normalize_config_tables(config);
config = validate_config_v9(config);
results = run_scenario_v9(config);
econ = run_health_economics_v9(results, econConfig);

jsonText = char(jsonencode(json_safe(econ)));
end

function config = normalize_config_tables(config)
if ~isstruct(config)
    return;
end

if isfield(config, 'ageDistributionTableRows') && ~isempty(config.ageDistributionTableRows)
    config.ageDistributionTable = struct2table(config.ageDistributionTableRows, 'AsArray', true);
elseif ~isfield(config, 'ageDistributionTable') || isempty(config.ageDistributionTable)
    config.ageDistributionTable = table();
end
end

function out = json_safe(value)
if istable(value)
    out = table2struct(value, 'ToScalar', false);
    out = json_safe(out);
elseif isstruct(value)
    out = struct_safe(value);
elseif iscell(value)
    out = cell(size(value));
    for i = 1:numel(value)
        out{i} = json_safe(value{i});
    end
elseif isstring(value)
    out = cellstr(value);
    if isscalar(value)
        out = out{1};
    end
elseif ischar(value) || isnumeric(value) || islogical(value)
    out = value;
else
    out = value;
end
end

function out = struct_safe(value)
if isempty(value)
    out = {};
    return;
end

if numel(value) > 1
    out = cell(numel(value), 1);
    for i = 1:numel(value)
        out{i} = json_safe(value(i));
    end
    return;
end

out = struct();
fields = fieldnames(value);
for i = 1:numel(fields)
    name = fields{i};
    out.(name) = json_safe(value.(name));
end
end
