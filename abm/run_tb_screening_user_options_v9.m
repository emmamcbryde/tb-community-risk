function results = run_tb_screening_user_options_v9()
%RUN_TB_SCREENING_USER_OPTIONS_V9 User-facing runner for tb_screening_mc_model_v9.
% Edit the USER INPUTS section below. Leave any field empty ([]) to use the
% default values already built into tb_screening_mc_model_v9.
%
% Risk-factor prevalence inputs can be supplied as:
%   - []            : use values from mock_data.csv
%   - scalar        : apply the same prevalence to all 3 broad age groups
%   - [a b c]       : use prevalences for [0-4, 5-14, >=15]
%
% Overall-only prevalence inputs:
%   - female, BCG   : scalar only
%
% Age distribution can be supplied either as:
%   - a file path (e.g. default_age_distribution.csv), or
%   - a table with one age-band column and one proportion column.
%
% The model still calibrates to the requested LTBI prevalence and active TB
% prevalence after applying your chosen inputs.

thisFile = mfilename('fullpath');
[thisDir, ~, ~] = fileparts(thisFile);
csvFile = fullfile(thisDir, 'default_data.csv');

%% ----------------------------- USER INPUTS ------------------------------
user = struct();

% Cohort / simulation size
user.N = 1500;
user.nReps = 2000;
user.seed = 1;

% Time horizon and screening setup
user.screenWindow = 2;
user.followHorizon = 20;
user.screenCoverage = 0.30;
user.screeningStrategy = 'prevent';   % 'random', 'ltbi', 'cure', 'prevent'

% Calibration targets
% Leave empty to use the Barry defaults embedded in the model.
user.ltbiPrevalence = [];             % e.g. 0.0753
user.activeTBPrevalence = [];         % e.g. 0.0130

% Age distribution
user.ageDistributionFile = '';        % '' = auto-detect default_age_distribution.csv next to default_data.csv
user.ageDistributionTable = table();  % optional manual table override
user.age85PlusMax = 89;

% Risk-factor prevalences
% [] = use mock_data.csv.
% scalar = same prevalence in all 3 age groups.
% [a b c] = [0-4, 5-14, >=15].
user.riskPrev = struct();
user.riskPrev.MJ = [];
user.riskPrev.contact = [];
user.riskPrev.renal = [];
user.riskPrev.diabetes = [];
user.riskPrev.smoking = [];
user.riskPrev.cld = [];
user.riskPrev.alcohol = [];
user.riskPrev.female = [];            % scalar overall only
user.riskPrev.BCG = [];               % scalar overall only

% Optional overrides for progression / reactivation odds ratios
% Leave empty ([]) to use the built-in defaults. The current default renal
% disease progression OR is set to 3.6 to reflect renal impairment rather
% than end-stage kidney disease / dialysis.
user.diseaseOR = struct();
user.diseaseOR.MJ = [];
user.diseaseOR.contact = [];
user.diseaseOR.renal = 3.6;
user.diseaseOR.diabetes = [];
user.diseaseOR.smoking = [];
user.diseaseOR.cld = [];
user.diseaseOR.alcohol = [];

% Test choice
user.testType = 'IGRA';               % 'IGRA' or 'TST'
user.testSensitivity = [];            % IGRA only; [] = default
user.testSpecificity = [];            % IGRA only; [] = default
user.tstSensitivity = [];             % TST only;  [] = default
user.tstSpecificityBCG = [];          % TST only;  [] = default
user.tstSpecificityNoBCG = [];        % TST only;  [] = default

% Treatment choice and cascade of care
user.regimen = '3HP';                 % '3HP', '4R', '3HR', '6H', '9H'
user.pStartTPT = [];                  % probability of starting after test-positive / eligible
user.regimenPComplete = [];           % marginal completion probability among starters
user.regimenADRstop = [];             % treatment-limiting ADR probability among starters
user.regimenEffFull = [];             % full-course protective efficacy
user.partialShortCourseMode = [];     % [] = default, or 'threshold80', 'linear', 'none'
user.partialDoseFractionADR = [];     % fraction of course taken if ADR stop occurs
user.partialDoseFractionOther = [];   % fraction of course taken if other stop occurs

% Optional natural-history tuning
user.earlyLateRatio = [];             % [] = default
%% -----------------------------------------------------------------------

clear functions
rehash

args = {};
args = add_if_set(args, 'N', user.N);
args = add_if_set(args, 'nReps', user.nReps);
args = add_if_set(args, 'seed', user.seed);
args = add_if_set(args, 'screenWindow', user.screenWindow);
args = add_if_set(args, 'followHorizon', user.followHorizon);
args = add_if_set(args, 'screenCoverage', user.screenCoverage);
args = add_if_set(args, 'screeningStrategy', user.screeningStrategy);
args = add_if_set(args, 'ltbiPrevalence', user.ltbiPrevalence);
args = add_if_set(args, 'activeTBPrevalence', user.activeTBPrevalence);
args = add_if_set(args, 'age85PlusMax', user.age85PlusMax);
args = add_if_set(args, 'testType', user.testType);
args = add_if_set(args, 'testSensitivity', user.testSensitivity);
args = add_if_set(args, 'testSpecificity', user.testSpecificity);
args = add_if_set(args, 'tstSensitivity', user.tstSensitivity);
args = add_if_set(args, 'tstSpecificityBCG', user.tstSpecificityBCG);
args = add_if_set(args, 'tstSpecificityNoBCG', user.tstSpecificityNoBCG);
args = add_if_set(args, 'regimen', user.regimen);
args = add_if_set(args, 'pStartTPT', user.pStartTPT);
args = add_if_set(args, 'regimenPComplete', user.regimenPComplete);
args = add_if_set(args, 'regimenADRstop', user.regimenADRstop);
args = add_if_set(args, 'regimenEffFull', user.regimenEffFull);
args = add_if_set(args, 'partialShortCourseMode', user.partialShortCourseMode);
args = add_if_set(args, 'partialDoseFractionADR', user.partialDoseFractionADR);
args = add_if_set(args, 'partialDoseFractionOther', user.partialDoseFractionOther);
args = add_if_set(args, 'earlyLateRatio', user.earlyLateRatio);

if istable(user.ageDistributionTable) && ~isempty(user.ageDistributionTable) && height(user.ageDistributionTable) > 0
    args = add_if_set(args, 'ageDistributionTable', user.ageDistributionTable);
elseif ~isempty(user.ageDistributionFile)
    args = add_if_set(args, 'ageDistributionFile', user.ageDistributionFile);
end

riskPrevOverrides = strip_empty_fields(user.riskPrev);
if ~isempty(fieldnames(riskPrevOverrides))
    args = add_if_set(args, 'riskPrevOverrides', riskPrevOverrides);
end

diseaseOROverrides = strip_empty_fields(user.diseaseOR);
if ~isempty(fieldnames(diseaseOROverrides))
    args = add_if_set(args, 'diseaseOROverrides', diseaseOROverrides);
end

results = tb_screening_mc_model_v9(csvFile, args{:});

% Concise on-screen summary
fprintf('\nChosen strategy\n');
disp(results.strategy)

fprintf('\nCalibration check\n');
fprintf('  Target LTBI prevalence:   %.4f\n', results.calibration.targetInfPrev);
fprintf('  Expected LTBI prevalence: %.4f\n', results.calibration.expectedInfPrev);
fprintf('  Target active TB prev:    %.4f\n', results.calibration.targetActive2y);
fprintf('  Expected active TB prev:  %.4f\n', results.calibration.expectedActive2y);

keepMetrics = [ ...
    "nScreened", "nTestPositiveNonActive", "nFalsePositiveTreated", ...
    "nTotalCoursesStarted", "nTotalCoursesCompleted", ...
    "nCuredInfection", "nPreventedActiveTB", ...
    "NNS_cureInfection", "NNS_preventActiveTB", ...
    "NNT_started_cureInfection", "NNT_started_preventActiveTB" ];

fprintf('\nKey summary metrics\n');
disp(results.summary(ismember(results.summary.Metric, keepMetrics), :));
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
