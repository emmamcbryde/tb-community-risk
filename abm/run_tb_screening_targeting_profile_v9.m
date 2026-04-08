function out = run_tb_screening_targeting_profile_v9(regimen, coverageGrid, Nprofile, seed)
% Profile which risk-factor categories are selected at each coverage cutoff
% for the v9 targeting rules that maximise LTBI cure or active-TB prevention.
%
% Default example:
%   - IGRA example context
%   - regimen = '3HP'
%   - coverageGrid = [0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.80 1.00]
%   - Nprofile = 50000
%
% Outputs:
%   out.overall              overall cohort composition
%   out.cumulative           cumulative composition at each coverage cutoff
%   out.incremental          incremental composition of newly added bands
%   out.topline              concise table with dominant age/risk factors
%   out.cumulativeCsv        file path written to CSV
%   out.incrementalCsv       file path written to CSV
%   out.toplineCsv           file path written to CSV
%
% Interpretation:
%   - cumulative: among everyone selected up to that cutoff, what is the mix?
%   - incremental: among the newly added people between two cutoffs, who are they?
%   - dominantRiskFactors is a simple label based on enriched factors among the
%     selected group or band; use the prevalence/enrichment columns for detail.
outDir = get_output_dir_v9();
%thisFile = mfilename('fullpath');
[outDir, ~, ~] = fileparts(thisFile);
csvFile = fullfile(outDir, 'default_data.csv');

if nargin < 1 || isempty(regimen)
    regimen = '3HP';
end
if nargin < 2 || isempty(coverageGrid)
    coverageGrid = [0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.80 1.00];
end
if nargin < 3 || isempty(Nprofile)
    Nprofile = 50000;
end
if nargin < 4 || isempty(seed)
    seed = 1;
end

base = tb_screening_mc_model_v9(csvFile, ...
    'N', Nprofile, ...
    'nReps', 1, ...
    'screenWindow', 2, ...
    'followHorizon', 20, ...
    'testType', 'IGRA', ...
    'regimen', regimen, ...
    'screeningStrategy', 'random', ...
    'screenCoverage', 1.0, ...
    'seed', seed);

T = base.exampleCohort;
N = height(T);
allIdx = true(N,1);
overall = summarise_selection(T, allIdx, struct(), 'overall', 0, 0, 1, NaN);

strategies = {'cure','prevent'};
nCum = numel(strategies) * numel(coverageGrid);
nBand = nCum;

cumRows = repmat(empty_profile_row(), nCum, 1);
bandRows = repmat(empty_profile_row(), nBand, 1);
rowCum = 0;
rowBand = 0;

overallRef = table2struct(overall);

for s = 1:numel(strategies)
    strategy = strategies{s};
    switch strategy
        case 'cure'
            score = T.cureTargetScore;
        case 'prevent'
            score = T.preventTargetScore;
        otherwise
            error('Unsupported strategy: %s', strategy);
    end

    [~, order] = sortrows([-score(:), (1:N)'], [1 2]);
    prevN = 0;
    prevCut = 0;

    for c = 1:numel(coverageGrid)
        cut = coverageGrid(c);
        nSel = max(min(round(cut * N), N), 0);
        selIdx = false(N,1);
        if nSel > 0
            selIdx(order(1:nSel)) = true;
            scoreThreshold = score(order(nSel));
        else
            scoreThreshold = NaN;
        end

        rowCum = rowCum + 1;
        cumRows(rowCum) = table2struct(summarise_selection(T, selIdx, overallRef, strategy, cut, 0, cut, scoreThreshold));

        bandIdx = false(N,1);
        if nSel > prevN
            bandIdx(order(prevN+1:nSel)) = true;
            bandThreshold = score(order(nSel));
        else
            bandThreshold = NaN;
        end
        rowBand = rowBand + 1;
        bandRows(rowBand) = table2struct(summarise_selection(T, bandIdx, overallRef, strategy, cut, prevCut, cut, bandThreshold));

        prevN = nSel;
        prevCut = cut;
    end
end

cumulative = struct2table(cumRows, 'AsArray', true);
incremental = struct2table(bandRows, 'AsArray', true);
cumulative = sortrows(cumulative, {'screeningStrategy','screenCoverageUpper'});
incremental = sortrows(incremental, {'screeningStrategy','screenCoverageUpper'});

keepCols = {'screeningStrategy','screenCoverageUpper','nSelected','scoreThreshold', ...
    'meanPInfection','meanDiseaseMultiplier','dominantAgeGroup','dominantRiskFactors', ...
    'propAge15plus','propContact','propMJ','propRenal','propDiabetes','propSmoking', ...
    'propCLD','propAlcohol','rrAge15plus','rrContact','rrMJ','rrRenal', ...
    'rrDiabetes','rrSmoking','rrCLD','rrAlcohol'};
topline = cumulative(:, keepCols);

cumCsv = fullfile(outDir, sprintf('igra_%s_targeting_profile_cumulative_v9.csv', lower(regimen)));
bandCsv = fullfile(outDir, sprintf('igra_%s_targeting_profile_incremental_v9.csv', lower(regimen)));
topCsv = fullfile(outDir, sprintf('igra_%s_targeting_profile_topline_v9.csv', lower(regimen)));

writetable(cumulative, cumCsv);
writetable(incremental, bandCsv);
writetable(topline, topCsv);

out = struct();
out.overall = overall;
out.cumulative = cumulative;
out.incremental = incremental;
out.topline = topline;
out.cumulativeCsv = cumCsv;
out.incrementalCsv = bandCsv;
out.toplineCsv = topCsv;

fprintf('Wrote cumulative targeting profile to: %s\n', cumCsv);
fprintf('Wrote incremental targeting profile to: %s\n', bandCsv);
fprintf('Wrote concise targeting profile to: %s\n', topCsv);
disp(topline)
end

function row = empty_profile_row()
row = struct( ...
    'screeningStrategy', '', ...
    'profileType', '', ...
    'screenCoverage', NaN, ...
    'screenCoverageLower', NaN, ...
    'screenCoverageUpper', NaN, ...
    'nSelected', NaN, ...
    'scoreThreshold', NaN, ...
    'meanPInfection', NaN, ...
    'meanDiseaseMultiplier', NaN, ...
    'propInfected', NaN, ...
    'propBCG', NaN, ...
    'propAge0_4', NaN, ...
    'propAge5_14', NaN, ...
    'propAge15plus', NaN, ...
    'propContact', NaN, ...
    'propMJ', NaN, ...
    'propRenal', NaN, ...
    'propDiabetes', NaN, ...
    'propSmoking', NaN, ...
    'propCLD', NaN, ...
    'propAlcohol', NaN, ...
    'rrAge0_4', NaN, ...
    'rrAge5_14', NaN, ...
    'rrAge15plus', NaN, ...
    'rrContact', NaN, ...
    'rrMJ', NaN, ...
    'rrRenal', NaN, ...
    'rrDiabetes', NaN, ...
    'rrSmoking', NaN, ...
    'rrCLD', NaN, ...
    'rrAlcohol', NaN, ...
    'dominantAgeGroup', '', ...
    'dominantRiskFactors', '');
end

function tbl = summarise_selection(T, idx, overallRef, strategy, screenCoverage, lowerCut, upperCut, scoreThreshold)
idx = logical(idx(:));
n = sum(idx);

if n == 0
    row = empty_profile_row();
    row.screeningStrategy = strategy;
    row.profileType = iff(lowerCut == 0, 'cumulative', 'incremental');
    row.screenCoverage = screenCoverage;
    row.screenCoverageLower = lowerCut;
    row.screenCoverageUpper = upperCut;
    row.nSelected = 0;
    row.scoreThreshold = scoreThreshold;
    tbl = struct2table(row, 'AsArray', true);
    return;
end

propAge = accumarray(T.ageGroup(idx), 1, [3 1]) / n;
propContact = mean(T.contact(idx));
propMJ = mean(T.MJ(idx));
propRenal = mean(T.renal(idx));
propDiabetes = mean(T.diabetes(idx));
propSmoking = mean(T.smoking(idx));
propCLD = mean(T.chronicLungDisease(idx));
propAlcohol = mean(T.alcoholDrugs(idx));
propBCG = mean(T.BCG(idx));
propInfected = mean(T.infected(idx));
meanPInfection = mean(T.pInfection(idx));
meanDiseaseMultiplier = mean(T.diseaseMultiplier(idx));

row = empty_profile_row();
row.screeningStrategy = strategy;
row.profileType = iff(lowerCut == 0, 'cumulative', 'incremental');
row.screenCoverage = screenCoverage;
row.screenCoverageLower = lowerCut;
row.screenCoverageUpper = upperCut;
row.nSelected = n;
row.scoreThreshold = scoreThreshold;
row.meanPInfection = meanPInfection;
row.meanDiseaseMultiplier = meanDiseaseMultiplier;
row.propInfected = propInfected;
row.propBCG = propBCG;
row.propAge0_4 = propAge(1);
row.propAge5_14 = propAge(2);
row.propAge15plus = propAge(3);
row.propContact = propContact;
row.propMJ = propMJ;
row.propRenal = propRenal;
row.propDiabetes = propDiabetes;
row.propSmoking = propSmoking;
row.propCLD = propCLD;
row.propAlcohol = propAlcohol;

if ~isempty(fieldnames(overallRef))
    row.rrAge0_4 = safe_ratio(row.propAge0_4, overallRef.propAge0_4);
    row.rrAge5_14 = safe_ratio(row.propAge5_14, overallRef.propAge5_14);
    row.rrAge15plus = safe_ratio(row.propAge15plus, overallRef.propAge15plus);
    row.rrContact = safe_ratio(row.propContact, overallRef.propContact);
    row.rrMJ = safe_ratio(row.propMJ, overallRef.propMJ);
    row.rrRenal = safe_ratio(row.propRenal, overallRef.propRenal);
    row.rrDiabetes = safe_ratio(row.propDiabetes, overallRef.propDiabetes);
    row.rrSmoking = safe_ratio(row.propSmoking, overallRef.propSmoking);
    row.rrCLD = safe_ratio(row.propCLD, overallRef.propCLD);
    row.rrAlcohol = safe_ratio(row.propAlcohol, overallRef.propAlcohol);

    ageLabels = {'0-4','5-14','15+'};
    ageRR = [row.rrAge0_4, row.rrAge5_14, row.rrAge15plus];
    [~, aIdx] = max(ageRR);
    row.dominantAgeGroup = ageLabels{aIdx};

    rfLabels = {'contact','MJ','renal','diabetes','smoking','CLD','alcohol'};
    rfProps = [row.propContact, row.propMJ, row.propRenal, row.propDiabetes, row.propSmoking, row.propCLD, row.propAlcohol];
    rfRR = [row.rrContact, row.rrMJ, row.rrRenal, row.rrDiabetes, row.rrSmoking, row.rrCLD, row.rrAlcohol];
    row.dominantRiskFactors = dominant_label(rfLabels, rfProps, rfRR);
else
    [~, aIdx] = max(propAge);
    ageLabels = {'0-4','5-14','15+'};
    row.dominantAgeGroup = ageLabels{aIdx};
    row.dominantRiskFactors = '';
end

tbl = struct2table(row, 'AsArray', true);
end

function out = dominant_label(labels, props, rr)
score = props .* rr;
valid = (props >= 0.03) & (rr > 1.02);
if ~any(valid)
    out = 'none';
    return;
end
score(~valid) = -Inf;
[~, ord] = sort(score, 'descend');
ord = ord(isfinite(score(ord)));
ord = ord(1:min(3, numel(ord)));
parts = cell(1, numel(ord));
for i = 1:numel(ord)
    parts{i} = sprintf('%s (%.2fx)', labels{ord(i)}, rr(ord(i)));
end
out = strjoin(parts, ', ');
end

function r = safe_ratio(a, b)
if isnan(a) || isnan(b) || b <= 0
    r = NaN;
else
    r = a / b;
end
end

function out = iff(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
