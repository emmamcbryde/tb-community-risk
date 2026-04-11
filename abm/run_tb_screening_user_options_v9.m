function results = run_tb_screening_user_options_v9()
%RUN_TB_SCREENING_USER_OPTIONS_V9 User-facing runner for tb_screening_mc_model_v9.
% Edit the USER INPUTS section below. Leave any field empty ([]) to use the
% canonical defaults returned by build_default_config_v9.
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

user = build_default_config_v9();

%% ----------------------------- USER INPUTS ------------------------------
% Defaults come from build_default_config_v9.

% Time horizon and screening setup
user.screeningStrategy = 'prevent';   % 'random', 'ltbi', 'cure', 'prevent'

% Calibration targets
% Leave empty to use the Barry defaults embedded in the model.
user.ltbiPrevalence = [];             % e.g. 0.0753
user.activeTBPrevalence = [];         % e.g. 0.0130

% Age distribution
user.ageDistributionFile = '';        % '' = auto-detect default_age_distribution.csv next to default_data.csv
user.ageDistributionTable = table();  % optional manual table override

% Risk-factor prevalences
% [] = use mock_data.csv.
% scalar = same prevalence in all 3 age groups.
% [a b c] = [0-4, 5-14, >=15].
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

results = run_scenario_v9(user);
summary = summarise_results_v9(results);

% Concise on-screen summary
fprintf('\nChosen strategy\n');
disp(summary.strategy)

fprintf('\nCalibration check\n');
fprintf('  Target LTBI prevalence:   %.4f\n', summary.calibration.targetInfPrev);
fprintf('  Expected LTBI prevalence: %.4f\n', summary.calibration.expectedInfPrev);
fprintf('  Target active TB prev:    %.4f\n', summary.calibration.targetActive2y);
fprintf('  Expected active TB prev:  %.4f\n', summary.calibration.expectedActive2y);

fprintf('\nKey summary metrics\n');
disp(summary.keyMetrics);
end
