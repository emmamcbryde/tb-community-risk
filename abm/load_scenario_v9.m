function [config, report, loadInfo] = load_scenario_v9(filename)
%LOAD_SCENARIO_V9 Load APY v9 scenario JSON with structured validation output.

if nargin < 1 || isempty(filename)
    error('load_scenario_v9 requires a filename.');
end

raw = fileread(filename);
payload = jsondecode(raw);

config = build_default_config_v9();
loadInfo = struct();
loadInfo.filename = filename;
loadInfo.contractVersion = '';
loadInfo.success = false;
loadInfo.messages = {};

if isstruct(payload) && isfield(payload, 'contractVersion')
    loadInfo.contractVersion = char(string(payload.contractVersion));
end

if ~isstruct(payload) || ~isfield(payload, 'config') || ~isstruct(payload.config)
    report = invalid_load_report('config', 'Scenario file', 'invalid_payload', ...
        'Scenario JSON does not contain a valid config object.');
    return;
end

config = merge_loaded_struct(config, payload.config);
if isfield(payload.config, 'ageDistributionTableRows')
    config.ageDistributionTable = rows_to_table(payload.config.ageDistributionTableRows);
end

report = collect_validation_issues_v9(config);
loadInfo.success = report.isValid;
if ~report.isValid
    loadInfo.messages = arrayfun(@(x) x.message, report.errors, 'UniformOutput', false);
    return;
end

config = validate_config_v9(config);
if isfield(payload, 'scenarioLabel') && ~isempty(payload.scenarioLabel)
    config.scenarioLabel = char(string(payload.scenarioLabel));
end
end

function s = merge_loaded_struct(base, loaded)
s = base;
fields = fieldnames(loaded);
for i = 1:numel(fields)
    name = fields{i};
    if strcmp(name, 'ageDistributionTableRows')
        continue;
    end
    if isfield(s, name) && isstruct(s.(name)) && isstruct(loaded.(name))
        s.(name) = merge_loaded_struct(s.(name), loaded.(name));
    else
        s.(name) = loaded.(name);
    end
end
end

function T = rows_to_table(rows)
if isempty(rows)
    T = table();
elseif isstruct(rows)
    T = struct2table(rows, 'AsArray', true);
else
    T = table();
end
end

function report = invalid_load_report(fieldName, fieldLabel, code, message)
issue = struct('field', fieldName, 'severity', 'error', 'code', code, 'message', message, 'fieldLabel', fieldLabel);
report = struct();
report.isValid = false;
report.hasWarnings = false;
report.errors = issue;
report.warnings = repmat(issue, 0, 1);
report.infos = repmat(issue, 0, 1);
report.fieldIssues = struct('load', issue);
report.fatalFieldNames = {fieldName};
report.warningFieldNames = {};
end
