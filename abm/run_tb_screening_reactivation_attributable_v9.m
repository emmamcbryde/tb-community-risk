function out = run_tb_screening_reactivation_attributable_v9(resultsOrCsv, Nattr)
% Estimate progression-only attributable risk for TB reactivation factors.
%
% This add-on isolates the contribution of each progression risk factor by
% comparing the modelled 20-year active-TB risk under the current disease
% multiplier with a counterfactual where one factor at a time is removed
% from the progression multiplier, while keeping infection risk fixed.
%
% Inputs:
%   resultsOrCsv : either a v9 RESULTS struct or a CSV filename.
%   Nattr        : size of the large natural-history cohort used for stable
%                  attributable-risk estimates (default 50000).
%
% Outputs:
%   out.summary                natural-history calibration summary
%   out.attributableTable      factor-by-factor attributable-risk table
%   out.naturalHistoryResults  the large no-intervention run used
%   out.figureFiles            PNGs for the attributable-risk charts
%
% CSVs also written:
%   tb_reactivation_attributable_v9.csv
%   tb_natural_history_check_v9.csv

if nargin < 2 || isempty(Nattr)
    Nattr = 50000;
end

nat = build_large_natural_history_run(resultsOrCsv, Nattr);
T = nat.exampleCohort;
settings = nat.settings;
cal = nat.calibration;

lambdaEarly = cal.lambdaEarly;
lambdaLate  = cal.lambdaLate;
sw = settings.screenWindow;
fh = settings.followHorizon;

% 20-year active-TB risk conditional on infection, under do-nothing.
risk20_full = active_risk_20y(T.diseaseMultiplier, lambdaEarly, lambdaLate, sw, fh);

% Use pInfection as the weight to keep infection risk fixed while isolating
% the progression contribution of each factor.
wInf = T.pInfection;

factorVars = {'contact','MJ','renal','diabetes','smoking','chronicLungDisease','alcoholDrugs'};
factorLabels = {'Close contact','Marijuana use','Renal disease','Diabetes', ...
    'Smoking','Chronic lung disease','Harmful alcohol use'};
ord = nat.parameters.disOR;
orVals = [ord.contact, ord.MJ, ord.renal, ord.diabetes, ord.smoking, ord.cld, ord.alcohol];

nFactors = numel(factorVars);
rows = repmat(struct(), nFactors, 1);

totalExpectedCases20yPer1500 = 1500 * sum(wInf .* risk20_full) / height(T);

for j = 1:nFactors
    exposed = logical(T.(factorVars{j}));
    mult_cf = T.diseaseMultiplier;
    mult_cf(exposed) = mult_cf(exposed) ./ orVals(j);
    risk20_cf = active_risk_20y(mult_cf, lambdaEarly, lambdaLate, sw, fh);
    diff20 = max(risk20_full - risk20_cf, 0);

    wExp = wInf .* double(exposed);

    rows(j).Factor = string(factorLabels{j});
    rows(j).RiskFactorVariable = string(factorVars{j});
    rows(j).ORused = orVals(j);
    rows(j).PrevalencePopulation = safe_fraction(sum(exposed), height(T));
    rows(j).ExpectedPrevalenceAmongInfected = safe_fraction(sum(wExp), sum(wInf));
    rows(j).ExpectedRisk20yAmongInfected_Exposed = safe_fraction(sum(wExp .* risk20_full), sum(wExp));
    rows(j).ExpectedRisk20yAmongInfected_Exposed_NoFactor = safe_fraction(sum(wExp .* risk20_cf), sum(wExp));
    rows(j).AttributableRisk20y_AmongExposedInfected = safe_fraction(sum(wExp .* diff20), sum(wExp));
    rows(j).AttributableFraction20y_AmongExposedInfected = safe_fraction(sum(wExp .* diff20), sum(wExp .* risk20_full));
    rows(j).PopulationAttributableFraction20y_ReactivationOnly = safe_fraction(sum(wInf .* diff20), sum(wInf .* risk20_full));
    rows(j).ExpectedAttributableCases20y_Per1500 = 1500 * sum(wInf .* diff20) / height(T);
    rows(j).ExpectedActiveCases20y_Per1500_DoNothing = totalExpectedCases20yPer1500;
end

attrib = struct2table(rows, 'AsArray', true);
attrib = sortrows(attrib, 'ExpectedAttributableCases20y_Per1500', 'descend');

% -------------------- Charts: attributable risk by factor --------------------
[thisDir, ~, ~] = fileparts(mfilename('fullpath'));

f1 = figure('Name','TB reactivation attributable cases','Color','w');
barh(attrib.ExpectedAttributableCases20y_Per1500);
set(gca, ...
    'YDir','reverse', ...
    'YTick',1:height(attrib), ...
    'YTickLabel',cellstr(attrib.Factor));
xlabel('Expected attributable active TB cases per 1500 over 20 years');
title('Reactivation risk factors: attributable active TB burden');
grid on;
pngCases = fullfile(thisDir, 'tb_reactivation_attributable_cases_v9.png');
try
    exportgraphics(f1, pngCases, 'Resolution', 300);
catch
    saveas(f1, pngCases);
end

f2 = figure('Name','TB reactivation attributable fraction','Color','w');
barh(100 * attrib.PopulationAttributableFraction20y_ReactivationOnly);
set(gca, ...
    'YDir','reverse', ...
    'YTick',1:height(attrib), ...
    'YTickLabel',cellstr(attrib.Factor));
xlabel('Population attributable fraction of 20-year active TB risk (%)');
title('Reactivation risk factors: attributable fraction');
grid on;
pngPaf = fullfile(thisDir, 'tb_reactivation_attributable_fraction_v9.png');

try
    exportgraphics(f2, pngPaf, 'Resolution', 300);
catch
    saveas(f2, pngPaf);
end

% Small natural-history summary to accompany the attributable table.
summaryRows = struct();
summaryRows.CohortSizeUsed = height(T);
summaryRows.TargetInfPrev = cal.targetInfPrev;
summaryRows.ModelExpectedInfPrev = cal.expectedInfPrev;
summaryRows.TargetAgeORge25 = cal.targetAgeOR;
summaryRows.ModelExpectedAgeORge25 = cal.expectedAgeOR;
summaryRows.TargetActive2y = cal.targetActive2y;
summaryRows.ModelExpectedActive2y = cal.expectedActive2y;
summaryRows.ExpectedActiveCases20y_Per1500_DoNothing = totalExpectedCases20yPer1500;
summary = struct2table(summaryRows, 'AsArray', true);

out = struct();
out.summary = summary;
out.attributableTable = attrib;
out.naturalHistoryResults = nat;
out.figureFiles = struct( ...
    'attributableCases', pngCases, ...
    'populationAttributableFraction', pngPaf);

try
    writetable(attrib, 'tb_reactivation_attributable_v9.csv');
    writetable(summary, 'tb_natural_history_check_v9.csv');
catch
end
end

function nat = build_large_natural_history_run(resultsOrCsv, Nattr)
thisFile = mfilename('fullpath');
[thisDir, ~, ~] = fileparts(thisFile);
defaultCsv = fullfile(thisDir, 'default_data.csv');

if nargin < 1 || isempty(resultsOrCsv)
    nat = tb_screening_mc_model_v9(defaultCsv, ...
        'N', Nattr, 'nReps', 1, 'screenCoverage', 0, 'screeningStrategy', 'random', ...
        'testType', 'IGRA', 'regimen', '3HP');
    return;
end

if isstruct(resultsOrCsv)
    base = resultsOrCsv;

    if isfield(base, 'settings') && isfield(base.settings, 'csvFile') && ~isempty(base.settings.csvFile)
        csvFile = char(base.settings.csvFile);
    else
        csvFile = defaultCsv;
    end

    args = {'N', Nattr, 'nReps', 1, 'screenCoverage', 0, 'screeningStrategy', 'random'};

    args = add_setting_arg(args, base.settings, 'seed');
    args = add_setting_arg(args, base.settings, 'screenWindow');
    args = add_setting_arg(args, base.settings, 'followHorizon');
    args = add_setting_arg(args, base.settings, 'earlyLateRatio');
    args = add_setting_arg(args, base.settings, 'targetInfPrev');
    args = add_setting_arg(args, base.settings, 'targetAgeOR');
    args = add_setting_arg(args, base.settings, 'targetActive2y');

    % Preserve exact-age generation and any user overrides that affect the
    % natural-history cohort composition.
    args = add_setting_arg(args, base.settings, 'ageDistributionFile');
    args = add_setting_arg(args, base.settings, 'ageDistributionTable');
    args = add_setting_arg(args, base.settings, 'ageDistributionSheet');
    args = add_setting_arg(args, base.settings, 'age85PlusMax');
    args = add_setting_arg(args, base.settings, 'riskPrevOverrides');
    args = add_setting_arg(args, base.settings, 'alcoholPrevByAge');

    % Carry through test / regimen choice for consistency, although these do
    % not change the no-screening natural history when screenCoverage = 0.
    args = add_setting_arg(args, base.settings, 'testType');
    args = add_setting_arg(args, base.settings, 'testSensitivity');
    args = add_setting_arg(args, base.settings, 'testSpecificity');
    args = add_setting_arg(args, base.settings, 'tstSensitivity');
    args = add_setting_arg(args, base.settings, 'tstSpecificityBCG');
    args = add_setting_arg(args, base.settings, 'tstSpecificityNoBCG');
    args = add_setting_arg(args, base.settings, 'regimen');
    args = add_setting_arg(args, base.settings, 'pStartTPT');
    args = add_setting_arg(args, base.settings, 'regimenPComplete');
    args = add_setting_arg(args, base.settings, 'regimenADRstop');
    args = add_setting_arg(args, base.settings, 'regimenEffFull');
    args = add_setting_arg(args, base.settings, 'partialShortCourseMode');
    args = add_setting_arg(args, base.settings, 'partialDoseFractionADR');
    args = add_setting_arg(args, base.settings, 'partialDoseFractionOther');

    nat = tb_screening_mc_model_v9(csvFile, args{:});
    return;
end

csvFile = char(resultsOrCsv);
nat = tb_screening_mc_model_v9(csvFile, ...
    'N', Nattr, 'nReps', 1, 'screenCoverage', 0, 'screeningStrategy', 'random', ...
    'testType', 'IGRA', 'regimen', '3HP');
end

function args = add_setting_arg(args, settings, name)
if ~isstruct(settings) || ~isfield(settings, name)
    return;
end
value = settings.(name);

if isempty(value)
    return;
end
if isnumeric(value) && isscalar(value) && isnan(value)
    return;
end
if ischar(value) && isempty(strtrim(value))
    return;
end
if isstring(value) && all(strlength(value) == 0)
    return;
end
if istable(value) && (isempty(value) || height(value) == 0)
    return;
end

args = [args, {name, value}]; %#ok<AGROW>
end

function risk20 = active_risk_20y(multDisease, lambdaEarly, lambdaLate, screenWindow, followHorizon)
risk20 = 1 - exp(-(lambdaEarly .* multDisease .* screenWindow) ...
    - (lambdaLate .* multDisease .* max(followHorizon - screenWindow, 0)));
end

function y = safe_fraction(num, den)
if den == 0
    y = NaN;
else
    y = num ./ den;
end
end
