function export_matlab_reference_scenarios_v9()
%EXPORT_MATLAB_REFERENCE_SCENARIOS_V9 Export compact MATLAB v9 reference outputs.
%
% This helper does not change the core APY model. It runs fixed reference
% scenarios and writes portable JSON/CSV files for Python validation.

try
    evalin('base', 'clear functions');
catch ME
    warning('export_matlab_reference_scenarios_v9:clearFunctionsFailed', ...
        'Could not clear functions before export: %s', ME.message);
end
rehash;

fprintf('MATLAB engine resolution for tb_screening_mc_model_v9:\n');
fprintf('%s\n', strtrim(evalc('which tb_screening_mc_model_v9 -all')));

thisDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
referenceRoot = fullfile(repoRoot, 'validation', 'matlab_reference');
suiteFile = fullfile(referenceRoot, 'scenario_suite_v1.json');
scenarios = load_reference_suite(suiteFile);

for i = 1:numel(scenarios)
    export_one_reference_scenario(referenceRoot, scenarios(i));
end
end

function export_one_reference_scenario(referenceRoot, scenario)
scenarioName = char(scenario.scenario_id);
outDir = fullfile(referenceRoot, scenarioName);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

config = build_default_config_v9();
config.scenarioLabel = scenarioName;
config = apply_config_overrides(config, scenario.config_overrides);

fprintf('Running MATLAB APY reference scenario: %s\n', scenarioName);
results = run_scenario_v9(config);
dn = run_tb_screening_do_nothing_v9(results);
bundle = build_results_bundle_v9(results, 'doNothing', dn);

write_json_file(fullfile(outDir, 'scenario_config.json'), config);
write_json_file(fullfile(outDir, 'matlab_results_bundle.json'), bundle);
write_json_file(fullfile(outDir, 'matlab_dynamic_comparison.json'), ...
    bundle.technical.dynamicComparison);
write_json_file(fullfile(outDir, 'manifest.json'), reference_manifest(scenario));

summary = summarise_results_v9(results);
writetable(summary.summaryTable, fullfile(outDir, 'matlab_summary.csv'));

if isfield(results, 'raw') && istable(results.raw)
    maxRawRows = 10000;
    if height(results.raw) <= maxRawRows
        writetable(results.raw, fullfile(outDir, 'matlab_raw_replicates.csv'));
    else
        fprintf('Skipping raw replicate export: %d rows exceeds limit of %d.\n', ...
            height(results.raw), maxRawRows);
    end
end

fprintf('Reference outputs written to: %s\n', outDir);
end

function scenarios = load_reference_suite(suiteFile)
if ~exist(suiteFile, 'file')
    error('Reference scenario suite not found: %s', suiteFile);
end
payload = jsondecode(fileread(suiteFile));
if isstruct(payload)
    scenarios = payload;
elseif iscell(payload)
    scenarios = [payload{:}];
else
    error('Reference scenario suite must be a JSON object array.');
end
end

function config = apply_config_overrides(config, overrides)
if isempty(overrides)
    return;
end
fields = fieldnames(overrides);
for i = 1:numel(fields)
    config.(fields{i}) = overrides.(fields{i});
end
end

function manifest = reference_manifest(scenario)
scenarioName = char(scenario.scenario_id);
manifest = struct();
manifest.scenario_id = scenarioName;
manifest.description = scenario.description;
manifest.model = 'MATLAB APY v9 reference';
manifest.purpose = 'Compact reference fixture for Python APY port validation';
manifest.expected_files = { ...
    'scenario_config.json', ...
    'matlab_dynamic_comparison.json', ...
    'matlab_summary.csv'};
manifest.excluded_large_files = { ...
    'matlab_results_bundle.json', ...
    'matlab_raw_replicates.csv'};
manifest.validation_focus = { ...
    'summary metrics', ...
    'technical.dynamicComparison', ...
    'output contract compatibility'};
manifest.tolerance_notes = ...
    'Exact individual-level equality is not expected because MATLAB and NumPy random streams differ.';
end

function snapshot = parameter_snapshot(results)
pars = results.parameters;
snapshot = struct();
snapshot.source = 'results.parameters / results.strategy / results.calibration';
snapshot.parameters = struct( ...
    'ageNames', {pars.ageNames}, ...
    'popFrac', pars.popFrac, ...
    'baseInfByAge', pars.baseInfByAge, ...
    'mjPrevByAge', pars.mjPrevByAge, ...
    'contactPrevByAge', pars.contactPrevByAge, ...
    'renalPrevByAge', pars.renalPrevByAge, ...
    'diabetesPrevByAge', pars.diabetesPrevByAge, ...
    'smokingPrevByAge', pars.smokingPrevByAge, ...
    'cldPrevByAge', pars.cldPrevByAge, ...
    'alcoholPrevByAge', pars.alcoholPrevByAge, ...
    'infOR', pars.infOR, ...
    'disOR', pars.disOR, ...
    'totalFemalePrev', pars.totalFemalePrev, ...
    'totalContactPrev', pars.totalContactPrev, ...
    'totalCurrentSmokerPrev', pars.totalCurrentSmokerPrev, ...
    'totalMJPrev', pars.totalMJPrev, ...
    'totalRenalPrev', pars.totalRenalPrev, ...
    'totalDiabetesPrev', pars.totalDiabetesPrev, ...
    'totalCLDPrev', pars.totalCLDPrev, ...
    'totalBCGPrev', pars.totalBCGPrev, ...
    'comorbidityCountCats', pars.comorbidityCountCats, ...
    'comorbidityCountProb', pars.comorbidityCountProb, ...
    'exactAgeValues', pars.exactAgeValues, ...
    'exactAgeProb', pars.exactAgeProb, ...
    'ageDistributionSource', pars.ageDistributionSource);
snapshot.strategy = results.strategy;
snapshot.calibration = results.calibration;
end

function write_json_file(filename, value)
payload = json_safe(value);
jsonText = char(jsonencode(payload));
fid = fopen(filename, 'w');
if fid < 0
    error('Could not open JSON output file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonText);
delete(cleanup);
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
elseif isnumeric(value) || islogical(value) || ischar(value)
    out = value;
else
    out = value;
end
end

function out = struct_safe(value)
if isempty(value)
    out = struct();
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
    out.(fields{i}) = json_safe(value.(fields{i}));
end
end
