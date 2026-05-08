function jsonText = collect_validation_issues_json_v9(config)
%COLLECT_VALIDATION_ISSUES_JSON_V9 Return v9 validation report as JSON text.
%
% MATLAB Engine for Python can fail to marshal nested/non-scalar struct arrays.
% Returning JSON keeps the Streamlit adapter contract JSON-like and portable.

config = normalize_config_tables(config);
report = collect_validation_issues_v9(config);

out = struct();
out.isValid = logical(report.isValid);
out.hasWarnings = logical(report.hasWarnings);
out.errors = issue_array(report.errors);
out.warnings = issue_array(report.warnings);
out.infos = issue_array(report.infos);
out.fieldIssues = field_issue_struct(report.fieldIssues);
out.fatalFieldNames = cellstr(report.fatalFieldNames);
out.warningFieldNames = cellstr(report.warningFieldNames);

jsonText = char(jsonencode(out));
end

function issues = issue_array(issues)
if isempty(issues)
    issues = {};
    return;
end
issues = reshape(issues, [], 1);
issues = num2cell(issues);
end

function out = field_issue_struct(fieldIssues)
out = struct();
if isempty(fieldIssues) || ~isstruct(fieldIssues)
    return;
end

fields = fieldnames(fieldIssues);
for i = 1:numel(fields)
    name = fields{i};
    out.(name) = issue_array(fieldIssues.(name));
end
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
