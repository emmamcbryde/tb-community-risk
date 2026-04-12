function config = validate_config_v9(config)
%VALIDATE_CONFIG_V9 MATLAB-facing wrapper for v9 config validation.

if ~isstruct(config)
    error('Config must be a struct.');
end

config = apply_metadata_defaults(config);
config = normalize_config_text_fields(config);

report = collect_validation_issues_v9(config);
if ~report.isValid
    error('validate_config_v9:InvalidConfig', '%s', build_error_message(report));
end
end

function config = apply_metadata_defaults(config)
defaults = build_default_config_v9();
defaultFields = fieldnames(defaults);
for i = 1:numel(defaultFields)
    fieldName = defaultFields{i};
    if ~isfield(config, fieldName) || isempty(config.(fieldName))
        config.(fieldName) = defaults.(fieldName);
    end
end

if ~isfield(config, 'sourceDataFiles') || ~isstruct(config.sourceDataFiles) || isempty(fieldnames(config.sourceDataFiles))
    config.sourceDataFiles = defaults.sourceDataFiles;
end

if isfield(config, 'csvFile') && ~isempty(config.csvFile)
    [~, name, ext] = fileparts(char(string(config.csvFile)));
    config.sourceDataFiles.tbDataFile = [name ext];
elseif ~isfield(config.sourceDataFiles, 'tbDataFile') || isempty(config.sourceDataFiles.tbDataFile)
    config.sourceDataFiles.tbDataFile = defaults.sourceDataFiles.tbDataFile;
end

if isfield(config, 'ageDistributionFile') && ~isempty(config.ageDistributionFile)
    [~, name, ext] = fileparts(char(string(config.ageDistributionFile)));
    config.sourceDataFiles.ageDistributionFile = [name ext];
elseif ~isfield(config.sourceDataFiles, 'ageDistributionFile') || isempty(config.sourceDataFiles.ageDistributionFile)
    config.sourceDataFiles.ageDistributionFile = defaults.sourceDataFiles.ageDistributionFile;
end
end

function config = normalize_config_text_fields(config)
config.csvFile = char(string(config.csvFile));
config.configVersion = char(string(config.configVersion));
config.modelVersion = char(string(config.modelVersion));
config.scenarioLabel = char(string(config.scenarioLabel));
config.sourceDataFiles.tbDataFile = char(string(config.sourceDataFiles.tbDataFile));
config.sourceDataFiles.ageDistributionFile = char(string(config.sourceDataFiles.ageDistributionFile));

if ~isempty(config.ageDistributionFile)
    config.ageDistributionFile = char(string(config.ageDistributionFile));
end
if isfield(config, 'testType') && ~isempty(config.testType)
    config.testType = char(string(config.testType));
end
if isfield(config, 'regimen') && ~isempty(config.regimen)
    config.regimen = char(string(config.regimen));
end
if isfield(config, 'screeningStrategy') && ~isempty(config.screeningStrategy)
    config.screeningStrategy = char(string(config.screeningStrategy));
end
if isfield(config, 'partialShortCourseMode') && ~isempty(config.partialShortCourseMode)
    config.partialShortCourseMode = char(string(config.partialShortCourseMode));
end
end

function msg = build_error_message(report)
lines = cell(1, numel(report.errors));
for i = 1:numel(report.errors)
    issue = report.errors(i);
    lines{i} = sprintf('[%s] %s: %s', issue.code, issue.fieldLabel, issue.message);
end
msg = strjoin(lines, newline);
end
