function comparison = run_tb_screening_targeted_gradient_v9(testType, regimen, nReps)
% Compare graduated targeting strategies from highest risk to lowest risk.
%
% Default comparison: IGRA + 3HP across cumulative screening coverages.
% Strategies:
%   - random  : random selection
%   - ltbi    : rank by predicted LTBI probability
%   - cure    : rank by predicted LTBI still latent at random screen time
%   - prevent : rank by predicted preventable active TB risk within 20 years
%
% This v9 runner also carries through the BCG-attributable TST outputs.

thisFile = mfilename('fullpath');
[thisDir, ~, ~] = fileparts(thisFile);
csvFile = fullfile(thisDir, 'default_data.csv');

if nargin < 1 || isempty(testType)
    testType = 'IGRA';
end
if nargin < 2 || isempty(regimen)
    regimen = '3HP';
end
if nargin < 3 || isempty(nReps)
    nReps = 2000;
end

coverageGrid = [0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.80 1.00];
strategies = {'random','ltbi','cure','prevent'};

nRows = numel(coverageGrid) * numel(strategies);
rows = repmat(struct( ...
    'testType', '', ...
    'regimen', '', ...
    'screeningStrategy', '', ...
    'screenCoverage', NaN, ...
    'medianScreened', NaN, ...
    'medianStartedTPT', NaN, ...
    'medianFalsePositiveTests', NaN, ...
    'medianFalsePositiveTreated', NaN, ...
    'medianFalsePositiveTreatedBCG', NaN, ...
    'medianExcessCoursesDueToBCG', NaN, ...
    'medianCuredInfection', NaN, ...
    'medianPreventedActiveTB', NaN, ...
    'medianBCGAttribFalsePositiveRate', NaN, ...
    'pooledNNS_cure', NaN, ...
    'pooledNNS_prevent', NaN, ...
    'pooledNNT_cure', NaN, ...
    'pooledNNT_prevent', NaN), nRows, 1);

row = 0;
for s = 1:numel(strategies)
    for c = 1:numel(coverageGrid)
        row = row + 1;
        res = tb_screening_mc_model_v9(csvFile, ...
            'N', 1500, ...
            'nReps', nReps, ...
            'screenWindow', 2, ...
            'followHorizon', 20, ...
            'testType', testType, ...
            'regimen', regimen, ...
            'screeningStrategy', strategies{s}, ...
            'screenCoverage', coverageGrid(c), ...
            'seed', 1);

        rows(row).testType = testType;
        rows(row).regimen = regimen;
        rows(row).screeningStrategy = strategies{s};
        rows(row).screenCoverage = coverageGrid(c);
        rows(row).medianScreened = get_metric(res.summary, 'nScreened');
        rows(row).medianStartedTPT = get_metric(res.summary, 'nStartTPT');
        rows(row).medianFalsePositiveTests = get_metric(res.summary, 'nFalsePositiveTests');
        rows(row).medianFalsePositiveTreated = get_metric(res.summary, 'nFalsePositiveTreated');
        rows(row).medianFalsePositiveTreatedBCG = get_metric(res.summary, 'nFalsePositiveTreatedBCG');
        rows(row).medianExcessCoursesDueToBCG = get_metric(res.summary, 'nExcessCoursesDueToBCG');
        rows(row).medianCuredInfection = get_metric(res.summary, 'nCuredInfection');
        rows(row).medianPreventedActiveTB = get_metric(res.summary, 'nPreventedActiveTB');
        rows(row).medianBCGAttribFalsePositiveRate = get_metric(res.summary, 'bcgAttributableFalsePositiveRateObserved');
        rows(row).pooledNNS_cure = sum(res.raw.nScreened) / sum(res.raw.nCuredInfection);
        rows(row).pooledNNS_prevent = sum(res.raw.nScreened) / sum(res.raw.nPreventedActiveTB);
        rows(row).pooledNNT_cure = sum(res.raw.nStartTPT) / sum(res.raw.nCuredInfection);
        rows(row).pooledNNT_prevent = sum(res.raw.nStartTPT) / sum(res.raw.nPreventedActiveTB);
    end
end

comparison = struct2table(rows, 'AsArray', true);
comparison = sortrows(comparison, {'screeningStrategy','screenCoverage'});
disp(comparison)
end

function val = get_metric(summaryTbl, metricName)
idx = strcmp(summaryTbl.Metric, metricName);
if any(idx)
    val = summaryTbl.Median(idx);
else
    val = NaN;
end
end
