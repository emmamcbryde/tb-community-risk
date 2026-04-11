function config = validate_config_v9(config)
%VALIDATE_CONFIG_V9 Lightweight validation for the v9 backend config.

if ~isstruct(config)
    error('Config must be a struct.');
end

config = apply_metadata_defaults(config);

requiredFields = { ...
    'csvFile', 'N', 'nReps', 'seed', 'screenWindow', 'followHorizon', ...
    'screenCoverage', 'screeningStrategy', 'ageDistributionFile', ...
    'ageDistributionTable', 'age85PlusMax', 'riskPrev', 'diseaseOR', ...
    'testType', 'regimen'};

for i = 1:numel(requiredFields)
    if ~isfield(config, requiredFields{i})
        error('Config is missing required field: %s', requiredFields{i});
    end
end

config.csvFile = char(string(config.csvFile));
if isempty(strtrim(config.csvFile))
    error('config.csvFile must be a non-empty path.');
end
if ~isfile(config.csvFile)
    error('config.csvFile not found: %s', config.csvFile);
end

validate_positive_scalar(config.N, 'config.N');
validate_positive_scalar(config.nReps, 'config.nReps');
validate_nonnegative_scalar(config.seed, 'config.seed');
validate_positive_scalar(config.screenWindow, 'config.screenWindow');
validate_positive_scalar(config.followHorizon, 'config.followHorizon');
validate_fraction_scalar(config.screenCoverage, 'config.screenCoverage');
validate_positive_scalar(config.age85PlusMax, 'config.age85PlusMax');

if config.followHorizon <= config.screenWindow
    error('config.followHorizon must be > config.screenWindow.');
end

if ~isempty(config.ageDistributionFile)
    config.ageDistributionFile = char(string(config.ageDistributionFile));
    if ~isfile(config.ageDistributionFile)
        error('config.ageDistributionFile not found: %s', config.ageDistributionFile);
    end
end

if ~istable(config.ageDistributionTable)
    error('config.ageDistributionTable must be a table.');
end

if ~isempty(config.ageDistributionFile) && height(config.ageDistributionTable) > 0
    error('Provide either config.ageDistributionFile or config.ageDistributionTable, not both.');
end

if ~isstruct(config.riskPrev)
    error('config.riskPrev must be a struct.');
end
if ~isstruct(config.diseaseOR)
    error('config.diseaseOR must be a struct.');
end

if ~ischar(config.configVersion) && ~(isstring(config.configVersion) && isscalar(config.configVersion))
    error('config.configVersion must be a character vector or scalar string.');
end
if ~ischar(config.modelVersion) && ~(isstring(config.modelVersion) && isscalar(config.modelVersion))
    error('config.modelVersion must be a character vector or scalar string.');
end
if ~ischar(config.scenarioLabel) && ~(isstring(config.scenarioLabel) && isscalar(config.scenarioLabel))
    error('config.scenarioLabel must be a character vector or scalar string.');
end
if ~(islogical(config.usesDefaults) && isscalar(config.usesDefaults))
    error('config.usesDefaults must be a logical scalar.');
end
if ~isstruct(config.sourceDataFiles)
    error('config.sourceDataFiles must be a struct.');
end
sourceRequired = {'tbDataFile', 'ageDistributionFile'};
for i = 1:numel(sourceRequired)
    if ~isfield(config.sourceDataFiles, sourceRequired{i})
        error('config.sourceDataFiles is missing required field: %s', sourceRequired{i});
    end
end

config.configVersion = char(string(config.configVersion));
config.modelVersion = char(string(config.modelVersion));
config.scenarioLabel = char(string(config.scenarioLabel));
config.sourceDataFiles.tbDataFile = char(string(config.sourceDataFiles.tbDataFile));
config.sourceDataFiles.ageDistributionFile = char(string(config.sourceDataFiles.ageDistributionFile));
end

function config = apply_metadata_defaults(config)
defaults = build_default_config_v9();
metaFields = {'configVersion', 'modelVersion', 'scenarioLabel', 'usesDefaults', 'sourceDataFiles'};
for i = 1:numel(metaFields)
    fieldName = metaFields{i};
    if ~isfield(config, fieldName) || isempty(config.(fieldName))
        config.(fieldName) = defaults.(fieldName);
    end
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

function validate_positive_scalar(x, name)
if ~(isnumeric(x) && isscalar(x) && isfinite(x) && x > 0)
    error('%s must be a positive finite scalar.', name);
end
end

function validate_nonnegative_scalar(x, name)
if ~(isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0)
    error('%s must be a non-negative finite scalar.', name);
end
end

function validate_fraction_scalar(x, name)
if ~(isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1)
    error('%s must be a scalar in [0, 1].', name);
end
end
