function results = run_scenario_v9(config)
%RUN_SCENARIO_V9 Validate a config struct and run the v9 reference engine.

config = validate_config_v9(config);
config.usesDefaults = is_default_execution_config(config);

args = {};
args = add_if_set(args, 'N', config.N);
args = add_if_set(args, 'nReps', config.nReps);
args = add_if_set(args, 'seed', config.seed);
args = add_if_set(args, 'screenWindow', config.screenWindow);
args = add_if_set(args, 'followHorizon', config.followHorizon);
args = add_if_set(args, 'screenCoverage', config.screenCoverage);
args = add_if_set(args, 'screeningStrategy', config.screeningStrategy);
args = add_if_set(args, 'ltbiPrevalence', config.ltbiPrevalence);
args = add_if_set(args, 'activeTBPrevalence', config.activeTBPrevalence);
args = add_if_set(args, 'targetAgeOR', config.targetAgeOR);
args = add_if_set(args, 'age85PlusMax', config.age85PlusMax);
args = add_if_set(args, 'ageDistributionSheet', config.ageDistributionSheet);
args = add_if_set(args, 'testType', config.testType);
args = add_if_set(args, 'testSensitivity', config.testSensitivity);
args = add_if_set(args, 'testSpecificity', config.testSpecificity);
args = add_if_set(args, 'tstSensitivity', config.tstSensitivity);
args = add_if_set(args, 'tstSpecificityBCG', config.tstSpecificityBCG);
args = add_if_set(args, 'tstSpecificityNoBCG', config.tstSpecificityNoBCG);
args = add_if_set(args, 'regimen', config.regimen);
args = add_if_set(args, 'pStartTPT', config.pStartTPT);
args = add_if_set(args, 'regimenPComplete', config.regimenPComplete);
args = add_if_set(args, 'regimenADRstop', config.regimenADRstop);
args = add_if_set(args, 'regimenEffFull', config.regimenEffFull);
args = add_if_set(args, 'partialShortCourseMode', config.partialShortCourseMode);
args = add_if_set(args, 'partialDoseFractionADR', config.partialDoseFractionADR);
args = add_if_set(args, 'partialDoseFractionOther', config.partialDoseFractionOther);
args = add_if_set(args, 'earlyLateRatio', config.earlyLateRatio);

if height(config.ageDistributionTable) > 0
    args = add_if_set(args, 'ageDistributionTable', config.ageDistributionTable);
elseif ~isempty(config.ageDistributionFile)
    args = add_if_set(args, 'ageDistributionFile', config.ageDistributionFile);
end

riskPrevOverrides = strip_empty_fields(config.riskPrev);
if ~isempty(fieldnames(riskPrevOverrides))
    args = add_if_set(args, 'riskPrevOverrides', riskPrevOverrides);
end

diseaseOROverrides = strip_empty_fields(config.diseaseOR);
if ~isempty(fieldnames(diseaseOROverrides))
    args = add_if_set(args, 'diseaseOROverrides', diseaseOROverrides);
end

results = tb_screening_mc_model_v9(config.csvFile, args{:});
results.interfaceConfig = config;
end

function tf = is_default_execution_config(config)
defaults = build_default_config_v9();
ignoreFields = {'scenarioLabel', 'usesDefaults'};
for i = 1:numel(ignoreFields)
    defaults = rmfield(defaults, ignoreFields{i});
    config = rmfield(config, ignoreFields{i});
end
tf = isequaln(config, defaults);
end

function args = add_if_set(args, name, value)
if isempty(value)
    return;
end
if isnumeric(value) && isscalar(value) && isnan(value)
    return;
end
if isstring(value) && strlength(value) == 0
    return;
end
if ischar(value) && isempty(strtrim(value))
    return;
end
if istable(value) && (isempty(value) || height(value) == 0)
    return;
end
args = [args, {name, value}]; %#ok<AGROW>
end

function s = strip_empty_fields(s)
if isempty(s)
    return;
end
fields = fieldnames(s);
remove = false(size(fields));
for i = 1:numel(fields)
    v = s.(fields{i});
    if isempty(v)
        remove(i) = true;
    elseif isnumeric(v) && all(isnan(v(:)))
        remove(i) = true;
    elseif ischar(v) && isempty(strtrim(v))
        remove(i) = true;
    elseif isstring(v) && all(strlength(v) == 0)
        remove(i) = true;
    end
end
if any(remove)
    s = rmfield(s, fields(remove));
end
end
