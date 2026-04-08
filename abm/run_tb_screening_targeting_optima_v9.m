function out = run_tb_screening_targeting_optima_v9(regimen, nReps, coverageGrid, Nprofile, seed)
% Summarise which categories dominate at the most efficient targeting cutoffs.
%
% Default: IGRA only, because this is intended as an explanatory add-on for
% the current example outputs.
%
% Outputs:
%   out.gradient         result table from run_tb_screening_targeted_gradient_v9
%   out.profile          result struct from run_tb_screening_targeting_profile_v9
%   out.optima           summary table of optimal cutoffs and dominant categories
%   out.optimaCsv        CSV file path
outDir = get_output_dir_v9();
if nargin < 1 || isempty(regimen)
    regimen = '3HP';
end
if nargin < 2 || isempty(nReps)
    nReps = 2000;
end
if nargin < 3 || isempty(coverageGrid)
    coverageGrid = [0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.80 1.00];
end
if nargin < 4 || isempty(Nprofile)
    Nprofile = 50000;
end
if nargin < 5 || isempty(seed)
    seed = 1;
end

thisFile = mfilename('fullpath');
[outDir, ~, ~] = fileparts(thisFile);

% Get efficiency frontier.
gradient = run_tb_screening_targeted_gradient_v9('IGRA', regimen, nReps);
if ~isempty(coverageGrid)
    gradient = gradient(ismember(round(gradient.screenCoverage,10), round(coverageGrid(:),10)), :);
end

% Get category composition at the same coverages.
profile = run_tb_screening_targeting_profile_v9(regimen, coverageGrid, Nprofile, seed);

% Build summary rows.
rows = {};
rows{end+1} = build_optimum_row('cure', 'pooledNNT_cure', gradient, profile); %#ok<AGROW>
rows{end+1} = build_optimum_row('prevent', 'pooledNNT_prevent', gradient, profile); %#ok<AGROW>
rows{end+1} = build_optimum_row('cure', 'pooledNNS_cure', gradient, profile); %#ok<AGROW>
rows{end+1} = build_optimum_row('prevent', 'pooledNNS_prevent', gradient, profile); %#ok<AGROW>

optima = vertcat(rows{:});
optima = sortrows(optima, {'objective','efficiencyMetric'});

optCsv = fullfile(outDir, sprintf('igra_%s_targeting_optima_v9.csv', lower(regimen)));
writetable(optima, optCsv);

out = struct();
out.gradient = gradient;
out.profile = profile;
out.optima = optima;
out.optimaCsv = optCsv;

fprintf('Wrote targeting-optima summary to: %s\n', optCsv);
disp(optima)
end

function tbl = build_optimum_row(strategyName, metricName, gradient, profile)
G = gradient(strcmpi(gradient.screeningStrategy, strategyName), :);
metric = G.(metricName);
valid = isfinite(metric) & ~isnan(metric);
if ~any(valid)
    tbl = empty_optimum_row(strategyName, metricName);
    return;
end

% Choose the minimum efficiency metric. Break ties by higher prevented/cured,
% then by lower coverage.
[bestMetric, ~] = min(metric(valid));
idxValid = find(valid);
idxBest = idxValid(metric(valid) == bestMetric);
if numel(idxBest) > 1
    if strcmpi(strategyName, 'cure')
        [~, j] = max(G.medianCuredInfection(idxBest));
    else
        [~, j] = max(G.medianPreventedActiveTB(idxBest));
    end
    idxBest = idxBest(j);
end
if numel(idxBest) > 1
    [~, j] = min(G.screenCoverage(idxBest));
    idxBest = idxBest(j);
end
best = G(idxBest, :);

% Random comparator at same coverage.
R = gradient(strcmpi(gradient.screeningStrategy, 'random') & abs(gradient.screenCoverage - best.screenCoverage) < 1e-9, :);
if isempty(R)
    R = empty_gradient_like(best);
end

% Cumulative and incremental profile rows at the same cutoff.
Pcum = profile.cumulative(strcmpi(profile.cumulative.screeningStrategy, strategyName) & ...
    abs(profile.cumulative.screenCoverageUpper - best.screenCoverage) < 1e-9, :);
Pinc = profile.incremental(strcmpi(profile.incremental.screeningStrategy, strategyName) & ...
    abs(profile.incremental.screenCoverageUpper - best.screenCoverage) < 1e-9, :);

if isempty(Pcum), Pcum = empty_profile_like(); end
if isempty(Pinc), Pinc = empty_profile_like(); end

row = table();
row.objective = string(strategyName);
row.efficiencyMetric = string(metricName);
row.optimalCoverage = best.screenCoverage;
row.targetedMetric = best.(metricName);
row.randomMetricSameCoverage = R.(metricName);
row.targetedMedianCuredInfection = best.medianCuredInfection;
row.randomMedianCuredInfectionSameCoverage = R.medianCuredInfection;
row.targetedMedianPreventedActiveTB = best.medianPreventedActiveTB;
row.randomMedianPreventedActiveTBSameCoverage = R.medianPreventedActiveTB;
row.targetedMedianStartedTPT = best.medianStartedTPT;
row.randomMedianStartedTPTSameCoverage = R.medianStartedTPT;
row.targetedDominantAgeGroup_Cumulative = string(Pcum.dominantAgeGroup);
row.targetedDominantRiskFactors_Cumulative = string(Pcum.dominantRiskFactors);
row.targetedDominantAgeGroup_Incremental = string(Pinc.dominantAgeGroup);
row.targetedDominantRiskFactors_Incremental = string(Pinc.dominantRiskFactors);
row.propAge15plus_Cumulative = Pcum.propAge15plus;
row.propContact_Cumulative = Pcum.propContact;
row.propMJ_Cumulative = Pcum.propMJ;
row.propRenal_Cumulative = Pcum.propRenal;
row.propDiabetes_Cumulative = Pcum.propDiabetes;
row.propSmoking_Cumulative = Pcum.propSmoking;
row.propAlcohol_Cumulative = Pcum.propAlcohol;
row.rrContact_Cumulative = Pcum.rrContact;
row.rrMJ_Cumulative = Pcum.rrMJ;
row.rrRenal_Cumulative = Pcum.rrRenal;
row.rrDiabetes_Cumulative = Pcum.rrDiabetes;
row.rrSmoking_Cumulative = Pcum.rrSmoking;
row.rrAlcohol_Cumulative = Pcum.rrAlcohol;
row.meanPInfection_Cumulative = Pcum.meanPInfection;
row.meanDiseaseMultiplier_Cumulative = Pcum.meanDiseaseMultiplier;
row.propContact_Incremental = Pinc.propContact;
row.propMJ_Incremental = Pinc.propMJ;
row.propRenal_Incremental = Pinc.propRenal;
row.propDiabetes_Incremental = Pinc.propDiabetes;
row.propSmoking_Incremental = Pinc.propSmoking;
row.propAlcohol_Incremental = Pinc.propAlcohol;
row.rrContact_Incremental = Pinc.rrContact;
row.rrMJ_Incremental = Pinc.rrMJ;
row.rrRenal_Incremental = Pinc.rrRenal;
row.rrDiabetes_Incremental = Pinc.rrDiabetes;
row.rrSmoking_Incremental = Pinc.rrSmoking;
row.rrAlcohol_Incremental = Pinc.rrAlcohol;

tbl = row;
end

function T = empty_optimum_row(strategyName, metricName)
T = table();
T.objective = string(strategyName);
T.efficiencyMetric = string(metricName);
T.optimalCoverage = NaN;
T.targetedMetric = NaN;
T.randomMetricSameCoverage = NaN;
T.targetedMedianCuredInfection = NaN;
T.randomMedianCuredInfectionSameCoverage = NaN;
T.targetedMedianPreventedActiveTB = NaN;
T.randomMedianPreventedActiveTBSameCoverage = NaN;
T.targetedMedianStartedTPT = NaN;
T.randomMedianStartedTPTSameCoverage = NaN;
T.targetedDominantAgeGroup_Cumulative = "";
T.targetedDominantRiskFactors_Cumulative = "";
T.targetedDominantAgeGroup_Incremental = "";
T.targetedDominantRiskFactors_Incremental = "";
vars = {'propAge15plus_Cumulative','propContact_Cumulative','propMJ_Cumulative','propRenal_Cumulative', ...
    'propDiabetes_Cumulative','propSmoking_Cumulative','propAlcohol_Cumulative', ...
    'rrContact_Cumulative','rrMJ_Cumulative','rrRenal_Cumulative','rrDiabetes_Cumulative', ...
    'rrSmoking_Cumulative','rrAlcohol_Cumulative','meanPInfection_Cumulative', ...
    'meanDiseaseMultiplier_Cumulative','propContact_Incremental','propMJ_Incremental', ...
    'propRenal_Incremental','propDiabetes_Incremental','propSmoking_Incremental', ...
    'propAlcohol_Incremental','rrContact_Incremental','rrMJ_Incremental', ...
    'rrRenal_Incremental','rrDiabetes_Incremental','rrSmoking_Incremental', ...
    'rrAlcohol_Incremental'};
for i = 1:numel(vars)
    T.(vars{i}) = NaN;
end
end

function R = empty_gradient_like(best)
R = best;
numVars = R.Properties.VariableNames(varfun(@isnumeric,R,'OutputFormat','uniform'));
for i = 1:numel(numVars)
    R.(numVars{i}) = NaN;
end
strVars = R.Properties.VariableNames(varfun(@iscellstr,R,'OutputFormat','uniform'));
for i = 1:numel(strVars)
    R.(strVars{i}) = {''};
end
if any(strcmp(R.Properties.VariableNames,'testType')), R.testType = ""; end
if any(strcmp(R.Properties.VariableNames,'regimen')), R.regimen = ""; end
if any(strcmp(R.Properties.VariableNames,'screeningStrategy')), R.screeningStrategy = ""; end
end

function P = empty_profile_like()
P = table();
P.dominantAgeGroup = "";
P.dominantRiskFactors = "";
vars = {'propAge15plus','propContact','propMJ','propRenal','propDiabetes','propSmoking','propAlcohol', ...
    'rrContact','rrMJ','rrRenal','rrDiabetes','rrSmoking','rrAlcohol', ...
    'meanPInfection','meanDiseaseMultiplier'};
for i = 1:numel(vars)
    P.(vars{i}) = NaN;
end
end
