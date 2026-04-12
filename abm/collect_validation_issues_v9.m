function report = collect_validation_issues_v9(config)
%COLLECT_VALIDATION_ISSUES_V9 Collect structured validation issues for v9 config.

report = init_validation_report();

if ~isstruct(config)
    report = add_issue(report, 'config', 'Config', 'error', 'invalid_type', ...
        'Config must be a struct.');
    report.isValid = false;
    return;
end

requiredFields = { ...
    'csvFile', 'N', 'nReps', 'seed', 'screenWindow', 'followHorizon', ...
    'screenCoverage', 'screeningStrategy', 'ageDistributionFile', ...
    'ageDistributionTable', 'age85PlusMax', 'targetAgeOR', 'riskPrev', ...
    'diseaseOR', 'testType', 'regimen'};

for i = 1:numel(requiredFields)
    fieldName = requiredFields{i};
    if ~isfield(config, fieldName)
        report = add_issue(report, fieldName, fieldName, 'error', 'required_missing', ...
            sprintf('Required field is missing: %s', fieldName));
    end
end

report = validate_existing_file(report, config, 'csvFile', 'Main data file');
report = validate_positive_scalar_field(report, config, 'N', 'Cohort size');
report = validate_positive_scalar_field(report, config, 'nReps', 'Number of replicates');
report = validate_nonnegative_scalar_field(report, config, 'seed', 'Random seed');
report = validate_positive_scalar_field(report, config, 'screenWindow', 'Screen window');
report = validate_positive_scalar_field(report, config, 'followHorizon', 'Follow-up horizon');
report = validate_fraction_field(report, config, 'screenCoverage', 'Screening coverage');
report = validate_positive_scalar_field(report, config, 'age85PlusMax', 'Age 85+ upper bound');
report = validate_positive_scalar_field(report, config, 'targetAgeOR', 'Age OR target');

if isfield(config, 'followHorizon') && isfield(config, 'screenWindow') && ...
        isnumeric(config.followHorizon) && isscalar(config.followHorizon) && ...
        isnumeric(config.screenWindow) && isscalar(config.screenWindow) && ...
        isfinite(config.followHorizon) && isfinite(config.screenWindow) && ...
        config.followHorizon <= config.screenWindow
    report = add_issue(report, 'followHorizon', 'Follow-up horizon', 'error', 'invalid_range', ...
        'followHorizon must be greater than screenWindow.');
end

report = validate_choice_field(report, config, 'testType', 'Test choice', {'IGRA', 'TST'});
report = validate_choice_field(report, config, 'regimen', 'Treatment choice', {'3HP', '4R', '3HR', '6H', '9H'});
report = validate_choice_field(report, config, 'screeningStrategy', 'Screening strategy', {'random', 'ltbi', 'cure', 'prevent'});
report = validate_optional_choice_field(report, config, 'partialShortCourseMode', 'Partial-course efficacy rule', {'threshold80', 'linear', 'none'});

report = validate_fraction_field(report, config, 'pStartTPT', 'Treatment start probability', true);
report = validate_fraction_field(report, config, 'regimenPComplete', 'Completion probability', true);
report = validate_fraction_field(report, config, 'regimenADRstop', 'ADR stop probability', true);
report = validate_fraction_field(report, config, 'regimenEffFull', 'Full-course efficacy', true);
report = validate_fraction_field(report, config, 'partialDoseFractionADR', 'Partial dose fraction after ADR stop', true);
report = validate_fraction_field(report, config, 'partialDoseFractionOther', 'Partial dose fraction after other stop', true);

report = validate_age_distribution_consistency(report, config);
report = validate_risk_prevalence_struct(report, config);

if isfield(config, 'riskPrev') && isstruct(config.riskPrev)
    report = validate_optional_fraction_field(report, config.riskPrev, 'female', 'Female prevalence');
    report = validate_optional_fraction_field(report, config.riskPrev, 'BCG', 'BCG prevalence');
end

if isfield(config, 'diseaseOR') && isstruct(config.diseaseOR)
    diseaseFields = {'contact', 'MJ', 'renal', 'diabetes', 'smoking', 'cld', 'alcohol'};
    diseaseLabels = {'Close contact progression OR', 'Marijuana progression OR', 'Renal disease progression OR', ...
        'Diabetes progression OR', 'Smoking progression OR', 'Chronic lung disease progression OR', ...
        'Alcohol / drugs progression OR'};
    for i = 1:numel(diseaseFields)
        report = validate_optional_positive_scalar_field(report, config.diseaseOR, diseaseFields{i}, diseaseLabels{i}, ...
            ['diseaseOR.' diseaseFields{i}]);
    end
end

if isfield(config, 'nReps') && isnumeric(config.nReps) && isscalar(config.nReps) && isfinite(config.nReps) && config.nReps > 1000
    report = add_issue(report, 'nReps', 'Number of replicates', 'warning', 'browser_large_nreps', ...
        'nReps is large for browser use and may lead to a slow UI response.');
end

report.hasWarnings = ~isempty(report.warnings);
report.isValid = isempty(report.errors);
end

function report = validate_age_distribution_consistency(report, config)
if isfield(config, 'ageDistributionTable') && ~istable(config.ageDistributionTable)
    report = add_issue(report, 'ageDistributionTable', 'Age distribution table', 'error', 'invalid_type', ...
        'ageDistributionTable must be a table.');
end

if isfield(config, 'ageDistributionFile') && ~isempty(config.ageDistributionFile)
    report = validate_existing_file(report, config, 'ageDistributionFile', 'Age distribution file');
end

if isfield(config, 'ageDistributionFile') && isfield(config, 'ageDistributionTable') && ...
        ~isempty(config.ageDistributionFile) && istable(config.ageDistributionTable) && height(config.ageDistributionTable) > 0
    report = add_issue(report, 'ageDistributionTable', 'Age distribution input', 'error', 'conflict', ...
        'Provide either ageDistributionFile or ageDistributionTable, not both.');
end
end

function report = validate_risk_prevalence_struct(report, config)
if ~isfield(config, 'riskPrev') || ~isstruct(config.riskPrev)
    report = add_issue(report, 'riskPrev', 'Risk-factor prevalences', 'error', 'invalid_type', ...
        'riskPrev must be a struct.');
    return;
end

fields = {'contact', 'MJ', 'renal', 'diabetes', 'smoking', 'cld', 'alcohol'};
labels = {'Close contact prevalence', 'Marijuana use prevalence', 'Renal disease prevalence', ...
    'Diabetes prevalence', 'Smoking prevalence', 'Chronic lung disease prevalence', 'Alcohol / drugs prevalence'};
for i = 1:numel(fields)
    fieldName = fields{i};
    label = labels{i};
    dotted = ['riskPrev.' fieldName];
    if ~isfield(config.riskPrev, fieldName) || isempty(config.riskPrev.(fieldName))
        continue;
    end
    value = config.riskPrev.(fieldName);
    if ~(isnumeric(value) && isvector(value))
        report = add_issue(report, dotted, label, 'error', 'invalid_type', ...
            'Value must be numeric.');
        continue;
    end
    if ~(isscalar(value) || numel(value) == 3)
        report = add_issue(report, dotted, label, 'error', 'invalid_shape', ...
            'Value must be either a scalar or a 3-element vector.');
    end
    if any(~isfinite(value(:))) || any(value(:) < 0) || any(value(:) > 1)
        report = add_issue(report, dotted, label, 'error', 'invalid_range', ...
            'Value must be within [0,1].');
    end
end
end

function report = validate_existing_file(report, parent, fieldName, fieldLabel)
if ~isfield(parent, fieldName)
    return;
end
value = parent.(fieldName);
if isempty(value)
    report = add_issue(report, fieldName, fieldLabel, 'error', 'required_empty', ...
        sprintf('%s must not be empty.', fieldLabel));
    return;
end
pathValue = char(string(value));
if isempty(strtrim(pathValue))
    report = add_issue(report, fieldName, fieldLabel, 'error', 'required_empty', ...
        sprintf('%s must not be empty.', fieldLabel));
elseif ~isfile(pathValue)
    report = add_issue(report, fieldName, fieldLabel, 'error', 'file_missing', ...
        sprintf('%s not found: %s', fieldLabel, pathValue));
end
end

function report = validate_fraction_field(report, parent, fieldName, fieldLabel, allowEmpty)
if nargin < 5
    allowEmpty = false;
end
if ~isfield(parent, fieldName)
    return;
end
value = parent.(fieldName);
if isempty(value)
    if allowEmpty
        return;
    end
    report = add_issue(report, fieldName, fieldLabel, 'error', 'required_empty', ...
        sprintf('%s must not be empty.', fieldLabel));
    return;
end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0 && value <= 1)
    report = add_issue(report, fieldName, fieldLabel, 'error', 'invalid_range', ...
        sprintf('%s must be a scalar in [0,1].', fieldLabel));
end
end

function report = validate_optional_fraction_field(report, parent, fieldName, fieldLabel)
if ~isfield(parent, fieldName) || isempty(parent.(fieldName))
    return;
end
value = parent.(fieldName);
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0 && value <= 1)
    report = add_issue(report, ['riskPrev.' fieldName], fieldLabel, 'error', 'invalid_range', ...
        sprintf('%s must be a scalar in [0,1].', fieldLabel));
end
end

function report = validate_choice_field(report, config, fieldName, fieldLabel, allowedValues)
if ~isfield(config, fieldName)
    return;
end
value = config.(fieldName);
if isempty(value)
    report = add_issue(report, fieldName, fieldLabel, 'error', 'required_empty', ...
        sprintf('%s must not be empty.', fieldLabel));
    return;
end
value = char(string(value));
if ~any(strcmpi(value, allowedValues))
    report = add_issue(report, fieldName, fieldLabel, 'error', 'invalid_choice', ...
        sprintf('%s must be one of: %s.', fieldLabel, strjoin(allowedValues, ', ')));
end
end

function report = validate_optional_choice_field(report, config, fieldName, fieldLabel, allowedValues)
if ~isfield(config, fieldName) || isempty(config.(fieldName))
    return;
end
value = char(string(config.(fieldName)));
if ~any(strcmpi(value, allowedValues))
    report = add_issue(report, fieldName, fieldLabel, 'error', 'invalid_choice', ...
        sprintf('%s must be one of: %s.', fieldLabel, strjoin(allowedValues, ', ')));
end
end

function report = validate_positive_scalar_field(report, parent, fieldName, fieldLabel)
if ~isfield(parent, fieldName)
    return;
end
value = parent.(fieldName);
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
    report = add_issue(report, fieldName, fieldLabel, 'error', 'invalid_positive_scalar', ...
        sprintf('%s must be a positive finite scalar.', fieldLabel));
end
end

function report = validate_nonnegative_scalar_field(report, parent, fieldName, fieldLabel)
if ~isfield(parent, fieldName)
    return;
end
value = parent.(fieldName);
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
    report = add_issue(report, fieldName, fieldLabel, 'error', 'invalid_nonnegative_scalar', ...
        sprintf('%s must be a non-negative finite scalar.', fieldLabel));
end
end

function report = validate_optional_positive_scalar_field(report, parent, fieldName, fieldLabel, fullFieldName)
if nargin < 5
    fullFieldName = fieldName;
end
if ~isfield(parent, fieldName) || isempty(parent.(fieldName))
    return;
end
value = parent.(fieldName);
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
    report = add_issue(report, fullFieldName, fieldLabel, 'error', 'invalid_positive_scalar', ...
        sprintf('%s must be a positive finite scalar.', fieldLabel));
end
end

function report = add_issue(report, fieldName, fieldLabel, severity, code, message)
issue = struct( ...
    'field', fieldName, ...
    'severity', severity, ...
    'code', code, ...
    'message', message, ...
    'fieldLabel', fieldLabel);

switch severity
    case 'error'
        report.errors(end+1) = issue; %#ok<AGROW>
        report.fatalFieldNames = add_unique_name(report.fatalFieldNames, fieldName);
    case 'warning'
        report.warnings(end+1) = issue; %#ok<AGROW>
        report.warningFieldNames = add_unique_name(report.warningFieldNames, fieldName);
    otherwise
        report.infos(end+1) = issue; %#ok<AGROW>
end

fieldKey = sanitize_field_name(fieldName);
if isfield(report.fieldIssues, fieldKey)
    report.fieldIssues.(fieldKey)(end+1) = issue; %#ok<AGROW>
else
    report.fieldIssues.(fieldKey) = issue;
end
end

function names = add_unique_name(names, name)
if ~any(strcmp(names, name))
    names{end+1} = name;
end
end

function key = sanitize_field_name(fieldName)
key = regexprep(fieldName, '[^a-zA-Z0-9_]', '_');
if isempty(key)
    key = 'root';
end
if isstrprop(key(1), 'digit')
    key = ['f_' key];
end
end

function report = init_validation_report()
emptyIssue = struct('field', '', 'severity', '', 'code', '', 'message', '', 'fieldLabel', '');
report = struct();
report.isValid = true;
report.hasWarnings = false;
report.errors = repmat(emptyIssue, 0, 1);
report.warnings = repmat(emptyIssue, 0, 1);
report.infos = repmat(emptyIssue, 0, 1);
report.fieldIssues = struct();
report.fatalFieldNames = {};
report.warningFieldNames = {};
end
