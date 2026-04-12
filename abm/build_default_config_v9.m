function config = build_default_config_v9()
%BUILD_DEFAULT_CONFIG_V9 Return the canonical default v9 scenario config.

thisFile = mfilename('fullpath');
[thisDir, ~, ~] = fileparts(thisFile);

config = struct();

% Stable scenario metadata
config.configVersion = 'apy_v9_config_v1';
config.modelVersion = 'v9';
config.scenarioLabel = 'APY default scenario';
config.usesDefaults = true;

% Input/output paths
config.csvFile = fullfile(thisDir, 'default_data.csv');
config.outputDir = get_output_dir_v9();
config.sourceDataFiles = struct( ...
    'tbDataFile', 'default_data.csv', ...
    'ageDistributionFile', 'default_age_distribution.csv');

% Cohort / simulation size
config.N = 1500;
config.nReps = 2000;
config.seed = 1;

% Time horizon and screening setup
config.screenWindow = 2;
config.followHorizon = 20;
config.screenCoverage = 0.30;
config.screeningStrategy = 'prevent';

% Calibration targets
config.ltbiPrevalence = [];
config.activeTBPrevalence = [];
config.targetAgeOR = 7.54;

% Age distribution
config.ageDistributionFile = '';
config.ageDistributionTable = table();
config.age85PlusMax = 89;
config.ageDistributionSheet = 1;

% Risk-factor prevalences
config.riskPrev = struct();
config.riskPrev.MJ = [];
config.riskPrev.contact = [];
config.riskPrev.renal = [];
config.riskPrev.diabetes = [];
config.riskPrev.smoking = [];
config.riskPrev.cld = [];
config.riskPrev.alcohol = [];
config.riskPrev.female = [];
config.riskPrev.BCG = [];

% Optional overrides for progression / reactivation odds ratios
config.diseaseOR = struct();
config.diseaseOR.MJ = [];
config.diseaseOR.contact = [];
config.diseaseOR.renal = 3.6;
config.diseaseOR.diabetes = [];
config.diseaseOR.smoking = [];
config.diseaseOR.cld = [];
config.diseaseOR.alcohol = [];

% Test choice
config.testType = 'IGRA';
config.testSensitivity = [];
config.testSpecificity = [];
config.tstSensitivity = [];
config.tstSpecificityBCG = [];
config.tstSpecificityNoBCG = [];

% Treatment choice and cascade of care
config.regimen = '3HP';
config.pStartTPT = [];
config.regimenPComplete = [];
config.regimenADRstop = [];
config.regimenEffFull = [];
config.partialShortCourseMode = [];
config.partialDoseFractionADR = [];
config.partialDoseFractionOther = [];

% Optional natural-history tuning
config.earlyLateRatio = [];
end
