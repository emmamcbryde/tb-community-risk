function report = validate_economics_config_v9(econConfig)
%VALIDATE_ECONOMICS_CONFIG_V9 Validate a v9 economics config struct.

report = init_report();

if nargin < 1 || isempty(econConfig) || ~isstruct(econConfig)
    report = add_issue(report, 'econConfig', 'error', 'Economics config must be a struct.');
    report.isValid = false;
    return;
end

if ~isfield(econConfig, 'metadata') || ~isstruct(econConfig.metadata)
    report = add_issue(report, 'metadata', 'error', 'metadata must be a struct.');
else
    report = validate_text_field(report, econConfig.metadata, 'currencyCode', 'metadata.currencyCode');
    report = validate_optional_numeric_scalar(report, econConfig.metadata, 'priceYear', 'metadata.priceYear');
    report = validate_text_field(report, econConfig.metadata, 'locationLabel', 'metadata.locationLabel');
    report = validate_text_field(report, econConfig.metadata, 'sourceNotes', 'metadata.sourceNotes');
    report = validate_text_field(report, econConfig.metadata, 'programCostBasis', 'metadata.programCostBasis');
end

if ~isfield(econConfig, 'costs') || ~isstruct(econConfig.costs)
    report = add_issue(report, 'costs', 'error', 'costs must be a struct.');
else
    costs = econConfig.costs;
    if ~isfield(costs, 'test') || ~isstruct(costs.test)
        report = add_issue(report, 'costs.test', 'error', 'costs.test must be a struct.');
    else
        report = validate_optional_cost(report, costs.test, 'IGRA', 'costs.test.IGRA');
        report = validate_optional_cost(report, costs.test, 'TST', 'costs.test.TST');
    end

    if ~isfield(costs, 'regimen') || ~isstruct(costs.regimen)
        report = add_issue(report, 'costs.regimen', 'error', 'costs.regimen must be a struct.');
    else
        regimenFields = {'x3HP', 'x4R', 'x3HR', 'x6H', 'x9H'};
        for i = 1:numel(regimenFields)
            report = validate_optional_cost(report, costs.regimen, regimenFields{i}, ...
                ['costs.regimen.' regimenFields{i}]);
        end
    end

    report = validate_optional_cost(report, costs, 'falsePositiveIncrementalPerPerson', ...
        'costs.falsePositiveIncrementalPerPerson');
    report = validate_optional_cost(report, costs, 'activeTBDiseasePerCase', ...
        'costs.activeTBDiseasePerCase');
    report = validate_optional_cost(report, costs, 'programSetupTotal', ...
        'costs.programSetupTotal');
    report = validate_optional_cost(report, costs, 'programRunningTotal', ...
        'costs.programRunningTotal');
end

report.isValid = isempty(report.errors);
report.hasWarnings = ~isempty(report.warnings);
end

function report = validate_text_field(report, parent, fieldName, fullName)
if ~isfield(parent, fieldName) || isempty(parent.(fieldName))
    return;
end
value = parent.(fieldName);
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    report = add_issue(report, fullName, 'error', [fullName ' must be text.']);
end
end

function report = validate_optional_numeric_scalar(report, parent, fieldName, fullName)
if ~isfield(parent, fieldName) || isempty(parent.(fieldName))
    return;
end
value = parent.(fieldName);
if ~(isnumeric(value) && isscalar(value) && isfinite(value))
    report = add_issue(report, fullName, 'error', [fullName ' must be a finite scalar.']);
end
end

function report = validate_optional_cost(report, parent, fieldName, fullName)
if ~isfield(parent, fieldName) || isempty(parent.(fieldName))
    return;
end
value = parent.(fieldName);
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
    report = add_issue(report, fullName, 'error', [fullName ' must be empty or a non-negative scalar.']);
end
end

function report = add_issue(report, fieldName, severity, message)
issue = struct('field', fieldName, 'severity', severity, 'message', message);
switch severity
    case 'error'
        report.errors(end+1) = issue;
    otherwise
        report.warnings(end+1) = issue;
end
end

function report = init_report()
emptyIssue = struct('field', '', 'severity', '', 'message', '');
report = struct();
report.isValid = true;
report.hasWarnings = false;
report.errors = repmat(emptyIssue, 0, 1);
report.warnings = repmat(emptyIssue, 0, 1);
end
