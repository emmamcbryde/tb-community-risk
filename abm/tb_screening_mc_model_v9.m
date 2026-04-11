function results = tb_screening_mc_model_v9(csvFile, varargin)
%TB_SCREENING_MC_MODEL_V9 Monte Carlo microsimulation for TB screening/TPT.
% v8 extends v7 with exact-age generation from an external age-distribution file
% and a complementary log-log LTBI model of the form:
%   p(LTBI) = 1 - exp(-exp(eta))
% where eta = log(lambda) + gamma*log(age+0.5) + betaX.
% The exact-age LTBI model is calibrated to reproduce:
%   - overall TB infection prevalence (Barry et al.)
%   - the crude age OR for >=25 vs <25 years (Barry et al.)
% while the active-TB natural-history hazard is separately calibrated to match
% the observed TB disease prevalence over the screening window.
%
% v7 extends v6 with:
%   - v4/v4b test-by-regimen strategy comparisons
%   - v4b BCG-attributable TST false-positive / extra-course outputs
%   - v5 risk-ranked targeted screening strategies
%
% Supported targeting modes:
%   - random  : random selection up to screenCoverage
%   - ltbi    : rank by predicted latent TB infection probability
%   - cure    : rank by predicted latent infection still present at screening
%   - prevent : rank by predicted preventable active TB risk within follow-up
% Supported tests:
%   - IGRA
%   - TST (with BCG-dependent specificity)
%
% Supported LTBI regimens:
%   - 3HP : 3 months once-weekly isoniazid + rifapentine (12 doses)
%   - 4R  : 4 months daily rifampin
%   - 3HR : 3 months daily isoniazid + rifampin
%   - 6H  : 6 months daily isoniazid
%   - 9H  : 9 months daily isoniazid
%
% Each regimen carries its own default assumptions for:
%   - completion probability among starters (marginal)
%   - treatment-limiting ADR probability among starters
%   - full-course protective efficacy
%   - partial-course efficacy rule
%
% The model remains a person-level closed-cohort microsimulation.
% Infection probability is assigned from age + selected risk factors.
% Time to active TB is generated from a piecewise exponential hazard.
% Screening is spread uniformly over the first 'screenWindow' years.
%
% Notes on interpretation:
%   - 'nCuredInfection' is a model-equivalent "protected/cured infection" state.
%     It represents people whose future TB progression is removed by treatment in
%     the simulation. For full courses this approximates sterilising protection;
%     for partial courses this is an efficacy assumption, not direct microbiologic proof.
%   - 'nPreventedActiveTB' is the number of untreated active TB cases that would
%     have occurred within the follow-up horizon but are prevented by treatment.

p = inputParser;
p.addRequired('csvFile', @(x) ischar(x) || isstring(x));
p.addParameter('N', 1500, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('nReps', 5000, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('screenWindow', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('followHorizon', 20, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('screenCoverage', 624/770, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('screeningStrategy', 'random', @(x) ischar(x) || isstring(x));

% Test strategy
p.addParameter('testType', 'IGRA', @(x) ischar(x) || isstring(x));
p.addParameter('testSensitivity', 0.95, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('testSpecificity', 0.98, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('tstSensitivity', 0.80, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('tstSpecificityBCG', 0.55, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('tstSpecificityNoBCG', 0.97, @(x) isnumeric(x) && isscalar(x));

% Treatment strategy / cascade
p.addParameter('regimen', '3HP', @(x) ischar(x) || isstring(x));
p.addParameter('pStartTPT', 0.85, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('regimenLibrary', struct(), @(x) isstruct(x));
p.addParameter('regimenPComplete', NaN, @(x) isnumeric(x) && isscalar(x));
p.addParameter('regimenADRstop', NaN, @(x) isnumeric(x) && isscalar(x));
p.addParameter('regimenEffFull', NaN, @(x) isnumeric(x) && isscalar(x));
p.addParameter('partialShortCourseMode', 'threshold80', @(x) ischar(x) || isstring(x));
p.addParameter('partialDoseFractionADR', 0.30, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('partialDoseFractionOther', 0.60, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);

% Legacy override parameters kept for backward compatibility
p.addParameter('pCompleteTPT', NaN, @(x) isnumeric(x) && isscalar(x));
p.addParameter('pSterilise', NaN, @(x) isnumeric(x) && isscalar(x));
p.addParameter('tptDuration', NaN, @(x) isnumeric(x) && isscalar(x));

% Natural history / calibration
p.addParameter('earlyLateRatio', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('targetInfPrev', 47/624, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('targetAgeOR', 7.54, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('targetActive2y', 10/770, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('ltbiPrevalence', NaN, @(x) isnumeric(x) && isscalar(x));
p.addParameter('activeTBPrevalence', NaN, @(x) isnumeric(x) && isscalar(x));
p.addParameter('ageDistributionFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('ageDistributionTable', table(), @(x) istable(x) || isempty(x));
p.addParameter('ageDistributionSheet', 1, @(x) (isnumeric(x) && isscalar(x)) || ischar(x) || isstring(x));
p.addParameter('riskPrevOverrides', struct(), @(x) isstruct(x));
p.addParameter('diseaseOROverrides', struct('renal', 3.6), @(x) isstruct(x));
p.addParameter('age85PlusMax', 89, @(x) isnumeric(x) && isscalar(x) && x >= 85);
p.addParameter('alcoholPrevByAge', [NaN NaN NaN], @(x) isnumeric(x) && numel(x) == 3);
p.addParameter('seed', 1, @(x) isnumeric(x) && isscalar(x));
p.parse(csvFile, varargin{:});
opts = p.Results;

if ~isnan(opts.ltbiPrevalence)
    opts.targetInfPrev = opts.ltbiPrevalence;
end
if ~isnan(opts.activeTBPrevalence)
    opts.targetActive2y = opts.activeTBPrevalence;
end

if opts.followHorizon <= opts.screenWindow
    error('followHorizon must be > screenWindow.');
end

opts.testType = upper(strtrim(string(opts.testType)));
if ~ismember(opts.testType, ["IGRA", "TST"])
    error('testType must be either "IGRA" or "TST".');
end

opts.screeningStrategy = lower(strtrim(string(opts.screeningStrategy)));
if ~ismember(opts.screeningStrategy, ["random", "ltbi", "cure", "prevent"])
    error('screeningStrategy must be one of "random", "ltbi", "cure", or "prevent".');
end
if isnan(opts.tstSpecificityNoBCG)
    opts.tstSpecificityNoBCG = opts.testSpecificity;
end

opts.regimen = upper(strtrim(string(opts.regimen)));
if ~ismember(opts.regimen, ["3HP", "4R", "3HR", "6H", "9H"])
    error('regimen must be one of "3HP", "4R", "3HR", "6H", or "9H".');
end
opts.partialShortCourseMode = lower(strtrim(string(opts.partialShortCourseMode)));
if ~ismember(opts.partialShortCourseMode, ["threshold80", "linear", "none"])
    error('partialShortCourseMode must be one of "threshold80", "linear", or "none".');
end

rng(opts.seed);

pars = load_tb_csv(csvFile);
ageFile = resolve_age_distribution_file(opts.ageDistributionFile, csvFile);
opts.ageDistributionFile = ageFile;
pars = load_age_distribution(pars, ageFile, opts.ageDistributionTable, opts.ageDistributionSheet, opts.age85PlusMax);
defaultAlcohol = [0; 0; 0.4];
if all(isnan(opts.alcoholPrevByAge))
    if ~isfield(pars, 'alcoholPrevByAge') || all(isnan(pars.alcoholPrevByAge))
        pars.alcoholPrevByAge = defaultAlcohol;
    else
        pars.alcoholPrevByAge = pars.alcoholPrevByAge(:);
        missAlcohol = isnan(pars.alcoholPrevByAge);
        pars.alcoholPrevByAge(missAlcohol) = defaultAlcohol(missAlcohol);
    end
else
    pars.alcoholPrevByAge = opts.alcoholPrevByAge(:);
end
pars = apply_risk_prevalence_overrides(pars, opts.riskPrevOverrides);
pars = apply_disease_or_overrides(pars, opts.diseaseOROverrides);
pars.totalFemalePrev = safe_default(pars.totalFemalePrev, 0.56);
pars.totalBCGPrev = safe_default(pars.totalBCGPrev, 0.68);
pars = refresh_total_prevalence_summaries(pars);

library = default_regimen_library(opts.partialShortCourseMode);
library = merge_regimen_library(library, opts.regimenLibrary);
reg = get_regimen_from_library(opts.regimen, library);
reg = apply_regimen_overrides(reg, opts);
validate_regimen(reg);

[ageInfLogLambda, ageInfGamma] = calibrate_age_infection_model(pars, opts.targetInfPrev, opts.targetAgeOR);
lambdaEarly = calibrate_early_hazard(pars, ageInfLogLambda, ageInfGamma, opts.targetActive2y, opts.screenWindow, opts.earlyLateRatio);
lambdaLate = lambdaEarly / opts.earlyLateRatio;

[firstRaw, exampleCohort] = simulate_one_cohort(opts.N, pars, reg, ageInfLogLambda, ageInfGamma, lambdaEarly, lambdaLate, opts);
raw = repmat(firstRaw, opts.nReps, 1);
raw(1) = firstRaw;
for r = 2:opts.nReps
    raw(r) = simulate_one_cohort(opts.N, pars, reg, ageInfLogLambda, ageInfGamma, lambdaEarly, lambdaLate, opts);
end

results = struct();
results.settings = opts;
results.parameters = pars;
results.regimenLibrary = library;
results.strategy = build_strategy_summary(opts, reg);
results.calibration = struct( ...
    'ageInfLogLambda', ageInfLogLambda, ...
    'ageInfGamma', ageInfGamma, ...
    'expectedInfPrev', expected_infection_prevalence_exact(ageInfLogLambda, ageInfGamma, pars), ...
    'expectedAgeOR', expected_age_or_exact(ageInfLogLambda, ageInfGamma, pars), ...
    'lambdaEarly', lambdaEarly, ...
    'lambdaLate', lambdaLate, ...
    'expectedActive2y', expected_active_within_window(lambdaEarly, pars, ageInfLogLambda, ageInfGamma, opts.screenWindow), ...
    'targetInfPrev', opts.targetInfPrev, ...
    'targetAgeOR', opts.targetAgeOR, ...
    'targetActive2y', opts.targetActive2y);
results.raw = struct2table(raw);
results.summary = summarise_results(results.raw);
results.exampleCohort = exampleCohort;
end

function pars = load_tb_csv(csvFile)
C = readcell(csvFile, 'Delimiter', ',');

pars = struct();
pars.ageNames = {'0-4','5-14','>=15'};
pars.popFrac = get_row_values(C, 1, 'population_fraction', 2:4);
pars.baseInfByAge = get_row_values(C, 1, 'infection_risk_by_age', 2:4);

pars.mjPrevByAge = get_row_values(C, 1, 'MJ_use', 2:4);
pars.contactPrevByAge = get_row_values(C, 1, 'contact', 2:4);
pars.renalPrevByAge = get_row_values(C, 1, 'renal', 2:4);
pars.diabetesPrevByAge = get_row_values(C, 1, 'diabetes', 2:4);
pars.smokingPrevByAge = get_row_values(C, 1, 'smoking', 2:4);
pars.cldPrevByAge = get_row_values(C, 1, 'chronic lung disease', 2:4);
pars.alcoholPrevByAge = get_row_values(C, 1, 'alcohol/drugs', 2:4);

pars.infOR = struct();
pars.infOR.MJ = safe_default(get_row_values(C, 1, 'MJ_use', 6, true), 1);
pars.infOR.contact = safe_default(get_row_values(C, 1, 'contact', 6, true), 1);
pars.infOR.renal = safe_default(get_row_values(C, 1, 'renal', 6, true), 1);

pars.disOR = struct();
pars.disOR.MJ = safe_default(get_row_values(C, 1, 'MJ_use', 5, true), 1);
pars.disOR.contact = safe_default(get_row_values(C, 1, 'contact', 5, true), 1);
pars.disOR.renal = safe_default(get_row_values(C, 1, 'renal', 5, true), 1);
pars.disOR.diabetes = safe_default(get_row_values(C, 1, 'diabetes', 5, true), 1);
pars.disOR.smoking = safe_default(get_row_values(C, 1, 'smoking', 5, true), 1);
pars.disOR.cld = safe_default(get_row_values(C, 1, 'chronic lung disease', 5, true), 1);
pars.disOR.alcohol = safe_default(get_row_values(C, 1, 'alcohol/drugs', 5, true), 1);

pars.totalFemalePrev = get_row_values(C, 9, 'female sex', 10, true);
pars.totalContactPrev = get_row_values(C, 9, 'close contact', 10, true);
pars.totalCurrentSmokerPrev = get_row_values(C, 9, 'current smoker', 10, true);
pars.totalMJPrev = get_row_values(C, 9, 'water-pipe marijuana', 10, true);
pars.totalRenalPrev = get_row_values(C, 9, 'renal disease', 10, true);
pars.totalDiabetesPrev = get_row_values(C, 9, 'diabetes', 10, true);
pars.totalCLDPrev = get_row_values(C, 9, 'chronic lung disease', 10, true);
pars.totalBCGPrev = get_row_values(C, 9, 'bcg vaccinated', 10, true);

pars.comorbidityCountCats = [0 1 2 3];
pars.comorbidityCountProb = [ ...
    get_row_values(C, 8, 'proportion comorbidity', 9), ...
    get_row_values(C, 8, 'proportion comorbidity', 10), ...
    get_row_values(C, 8, 'proportion comorbidity', 11), ...
    get_row_values(C, 8, 'proportion comorbidity', 12)];

pars.popFrac = pars.popFrac(:);
pars.baseInfByAge = pars.baseInfByAge(:);
pars.mjPrevByAge = pars.mjPrevByAge(:);
pars.contactPrevByAge = pars.contactPrevByAge(:);
pars.renalPrevByAge = pars.renalPrevByAge(:);
pars.diabetesPrevByAge = pars.diabetesPrevByAge(:);
pars.smokingPrevByAge = pars.smokingPrevByAge(:);
pars.cldPrevByAge = pars.cldPrevByAge(:);
pars.alcoholPrevByAge = pars.alcoholPrevByAge(:);
end

function pars = apply_risk_prevalence_overrides(pars, overrides)
if isempty(overrides)
    return;
end
fields = fieldnames(overrides);
if isempty(fields)
    return;
end
for i = 1:numel(fields)
    rawName = fields{i};
    val = overrides.(rawName);
    if isempty(val)
        continue;
    end
    if isnumeric(val) && all(isnan(val(:)))
        continue;
    end
    name = lower(regexprep(strtrim(rawName), '[^a-z0-9]', ''));
    switch name
        case {'mj','marijuana','waterpipemarijuana','waterpipemarijuanause'}
            pars.mjPrevByAge = normalize_age_prevalence_override(val, rawName);
        case {'contact','closecontact','tbcontact','closecontactwithtbcase'}
            pars.contactPrevByAge = normalize_age_prevalence_override(val, rawName);
        case {'renal','renaldisease'}
            pars.renalPrevByAge = normalize_age_prevalence_override(val, rawName);
        case {'diabetes'}
            pars.diabetesPrevByAge = normalize_age_prevalence_override(val, rawName);
        case {'smoking','currentsmoker','eversmoker'}
            pars.smokingPrevByAge = normalize_age_prevalence_override(val, rawName);
        case {'cld','chroniclungdisease','chroniclung'}
            pars.cldPrevByAge = normalize_age_prevalence_override(val, rawName);
        case {'alcohol','alcoholdrugs','harmfulalcohol'}
            pars.alcoholPrevByAge = normalize_age_prevalence_override(val, rawName);
        case {'female','femalesex','sexfemale'}
            pars.totalFemalePrev = normalize_overall_prevalence_override(val, rawName);
        case {'bcg','priorbcg','bcgvaccinated'}
            pars.totalBCGPrev = normalize_overall_prevalence_override(val, rawName);
        otherwise
            warning('tb_screening_mc_model_v9:UnknownRiskOverride', ...
                'Unknown riskPrevOverrides field "%s". Ignoring it.', rawName);
    end
end
end


function pars = apply_disease_or_overrides(pars, overrides)
if isempty(overrides)
    return;
end
fields = fieldnames(overrides);
if isempty(fields)
    return;
end
for i = 1:numel(fields)
    rawName = fields{i};
    val = overrides.(rawName);
    if isempty(val)
        continue;
    end
    if isnumeric(val) && all(isnan(val(:)))
        continue;
    end
    name = lower(regexprep(strtrim(rawName), '[^a-z0-9]', ''));
    switch name
        case {'mj','marijuana','waterpipemarijuana','waterpipemarijuanause'}
            pars.disOR.MJ = normalize_positive_or_override(val, rawName);
        case {'contact','closecontact','tbcontact','closecontactwithtbcase'}
            pars.disOR.contact = normalize_positive_or_override(val, rawName);
        case {'renal','renaldisease'}
            pars.disOR.renal = normalize_positive_or_override(val, rawName);
        case {'diabetes'}
            pars.disOR.diabetes = normalize_positive_or_override(val, rawName);
        case {'smoking','currentsmoker','eversmoker'}
            pars.disOR.smoking = normalize_positive_or_override(val, rawName);
        case {'cld','chroniclungdisease','chroniclung'}
            pars.disOR.cld = normalize_positive_or_override(val, rawName);
        case {'alcohol','alcoholdrugs','harmfulalcohol'}
            pars.disOR.alcohol = normalize_positive_or_override(val, rawName);
        otherwise
            warning('tb_screening_mc_model_v9:UnknownDiseaseOROverride', ...
                'Unknown diseaseOROverrides field "%s". Ignoring it.', rawName);
    end
end
end

function val = normalize_positive_or_override(v, label)
v = double(v(:));
if ~isscalar(v)
    error('Disease OR override for %s must be a positive scalar.', label);
end
val = v;
if isnan(val) || val <= 0
    error('Disease OR override for %s must be a positive scalar.', label);
end
end

function vec = normalize_age_prevalence_override(val, label)
val = double(val(:));
if isscalar(val)
    vec = repmat(val, 3, 1);
elseif numel(val) == 3
    vec = val(:);
else
    error('Risk prevalence override for %s must be a scalar or a 3-element vector [0-4, 5-14, >=15].', label);
end
if any(isnan(vec)) || any(vec < 0) || any(vec > 1)
    error('Risk prevalence override for %s must be probabilities between 0 and 1.', label);
end
end

function val = normalize_overall_prevalence_override(v, label)
v = double(v(:));
if ~isscalar(v)
    error('Overall prevalence override for %s must be a scalar probability between 0 and 1.', label);
end
val = v;
if isnan(val) || val < 0 || val > 1
    error('Overall prevalence override for %s must be a scalar probability between 0 and 1.', label);
end
end

function pars = refresh_total_prevalence_summaries(pars)
if isfield(pars, 'exactAgeValues') && isfield(pars, 'exactAgeProb') ...
        && ~isempty(pars.exactAgeValues) && ~isempty(pars.exactAgeProb)
    grp = broad_age_group_from_years(pars.exactAgeValues);
    w = pars.exactAgeProb(:);
    pars.popFrac = [sum(w(grp == 1)); sum(w(grp == 2)); sum(w(grp == 3))];
else
    grp = [ones(5,1); 2*ones(10,1); 3*ones(75,1)];
    w = build_fallback_exact_age_prob((0:89)', pars.popFrac);
end
pars.totalMJPrev = sum(age_lookup(pars.mjPrevByAge, grp) .* w);
pars.totalContactPrev = sum(age_lookup(pars.contactPrevByAge, grp) .* w);
pars.totalCurrentSmokerPrev = sum(age_lookup(pars.smokingPrevByAge, grp) .* w);
pars.totalRenalPrev = sum(age_lookup(pars.renalPrevByAge, grp) .* w);
pars.totalDiabetesPrev = sum(age_lookup(pars.diabetesPrevByAge, grp) .* w);
pars.totalCLDPrev = sum(age_lookup(pars.cldPrevByAge, grp) .* w);
if ~isfield(pars, 'totalFemalePrev') || isnan(pars.totalFemalePrev)
    pars.totalFemalePrev = 0.56;
end
if ~isfield(pars, 'totalBCGPrev') || isnan(pars.totalBCGPrev)
    pars.totalBCGPrev = 0.68;
end
end

function library = default_regimen_library(shortMode)
library = struct();

% 3HP base-case choices are deliberately conservative and intended for
% program planning rather than as fixed truths. They are anchored to CDC's
% current recommendation of 3HP as a preferred short-course regimen and to
% trial/programme evidence showing high completion but somewhat higher ADR
% discontinuation than 4R. Full-course efficacy is kept aligned with other
% short-course regimens because the best available evidence suggests similar
% prevention of TB disease versus 9H and no difference versus 4R.
library.r3HP = make_regimen('3HP', 3, 12, 0.80, 0.05, 0.85, shortMode, ...
    [0 (10/12) (11/12) 1.00], [0 0.00 0.85 0.85], (11/12));
library.r4R = make_regimen('4R', 4, 120, 0.78, 0.02, 0.85, shortMode, ...
    [0 0.80 1.00], [0 0.00 0.85], 0.80);
library.r3HR = make_regimen('3HR', 3, 90, 0.80, 0.03, 0.85, shortMode, ...
    [0 0.80 1.00], [0 0.00 0.85], 0.80);
library.r6H = make_regimen('6H', 6, 180, 0.67, 0.03, 0.65, 'knots', ...
    [0 0.50 1.00], [0 0.30 0.65], NaN);
library.r9H = make_regimen('9H', 9, 270, 0.60, 0.04, 0.85, 'knots', ...
    [0 (90/270) (180/270) 1.00], [0 0.30 0.65 0.85], NaN);
end

function reg = make_regimen(label, months, targetDoses, pComplete, pADRstop, effFull, partialMode, partialDoseKnots, partialEffKnots, partialThreshold)
reg = struct();
reg.label = char(label);
reg.months = months;
reg.targetDoses = targetDoses;
reg.pComplete = pComplete;
reg.pADRstop = pADRstop;
reg.effFull = effFull;
reg.partialMode = char(partialMode);
reg.partialDoseKnots = partialDoseKnots(:)';
reg.partialEffKnots = partialEffKnots(:)';
reg.partialThreshold = partialThreshold;
end

function library = merge_regimen_library(library, customLibrary)
if isempty(fieldnames(customLibrary))
    return;
end
customFields = fieldnames(customLibrary);
for i = 1:numel(customFields)
    f = customFields{i};
    library.(f) = customLibrary.(f);
end
end

function reg = get_regimen_from_library(label, library)
field = regimen_field_name(label);
if ~isfield(library, field)
    error('Regimen %s not found in regimenLibrary.', char(label));
end
reg = library.(field);
end

function reg = apply_regimen_overrides(reg, opts)
if ~isnan(opts.regimenPComplete)
    reg.pComplete = opts.regimenPComplete;
end
if ~isnan(opts.regimenADRstop)
    reg.pADRstop = opts.regimenADRstop;
end
if ~isnan(opts.regimenEffFull)
    reg.effFull = opts.regimenEffFull;
end

% Backward-compatible overrides
if ~isnan(opts.pCompleteTPT)
    reg.pComplete = opts.pCompleteTPT;
end
if ~isnan(opts.pSterilise)
    reg.effFull = opts.pSterilise;
end
if ~isnan(opts.tptDuration)
    reg.months = 12 * opts.tptDuration;
    reg.targetDoses = round(30 * reg.months);
end
end

function validate_regimen(reg)
fieldsNeeded = {'label','months','targetDoses','pComplete','pADRstop','effFull','partialMode','partialDoseKnots','partialEffKnots','partialThreshold'};
for i = 1:numel(fieldsNeeded)
    if ~isfield(reg, fieldsNeeded{i})
        error('Regimen definition is missing field "%s".', fieldsNeeded{i});
    end
end
if reg.months <= 0 || reg.targetDoses <= 0
    error('Regimen months and targetDoses must be > 0.');
end
if reg.pComplete < 0 || reg.pComplete > 1
    error('Regimen pComplete must be in [0,1].');
end
if reg.pADRstop < 0 || reg.pADRstop > 1
    error('Regimen pADRstop must be in [0,1].');
end
if reg.effFull < 0 || reg.effFull > 1
    error('Regimen effFull must be in [0,1].');
end
if reg.pComplete > (1 - reg.pADRstop + 1e-12)
    error(['Regimen pComplete exceeds the maximum possible completion after accounting for ' ...
        'treatment-limiting ADR. Reduce pComplete or pADRstop.']);
end
if isempty(reg.partialDoseKnots) || isempty(reg.partialEffKnots)
    error('Regimen partialDoseKnots and partialEffKnots cannot be empty.');
end
if numel(reg.partialDoseKnots) ~= numel(reg.partialEffKnots)
    error('Regimen partialDoseKnots and partialEffKnots must have the same length.');
end
end

function field = regimen_field_name(label)
label = upper(strtrim(char(label)));
switch label
    case '3HP'
        field = 'r3HP';
    case '4R'
        field = 'r4R';
    case '3HR'
        field = 'r3HR';
    case '6H'
        field = 'r6H';
    case '9H'
        field = 'r9H';
    otherwise
        error('Unknown regimen label: %s', label);
end
end

function strategy = build_strategy_summary(opts, reg)
strategy = struct();
strategy.testType = char(opts.testType);
strategy.screeningStrategy = char(opts.screeningStrategy);
strategy.regimen = reg.label;
strategy.regimenMonths = reg.months;
strategy.targetDoses = reg.targetDoses;
strategy.pStartTPT = opts.pStartTPT;
strategy.pComplete = reg.pComplete;
strategy.pADRstop = reg.pADRstop;
strategy.effFull = reg.effFull;
strategy.partialMode = reg.partialMode;
strategy.partialDoseFractionADR = opts.partialDoseFractionADR;
strategy.partialDoseFractionOther = opts.partialDoseFractionOther;
strategy.testSensitivityIGRA = opts.testSensitivity;
strategy.testSpecificityIGRA = opts.testSpecificity;
strategy.tstSensitivity = opts.tstSensitivity;
strategy.tstSpecificityBCG = opts.tstSpecificityBCG;
strategy.tstSpecificityNoBCG = opts.tstSpecificityNoBCG;
strategy.bcgSpecificityPenaltyTST = max(0, opts.tstSpecificityNoBCG - opts.tstSpecificityBCG);
strategy.partialEfficacyAt50pct = regimen_partial_efficacy(reg, 0.50);
strategy.partialEfficacyAt80pct = regimen_partial_efficacy(reg, 0.80);
strategy.partialEfficacyAt100pct = reg.effFull;
end

function value = get_row_values(C, searchCol, label, outCols, allowMissing)
if nargin < 5
    allowMissing = false;
end
colText = lower(string(C(:, searchCol)));
row = find(contains(colText, lower(label)), 1, 'first');
if isempty(row)
    if allowMissing
        if numel(outCols) == 1
            value = NaN;
        else
            value = NaN(numel(outCols),1);
        end
        return;
    else
        error('Could not find label "%s" in column %d.', label, searchCol);
    end
end

if numel(outCols) == 1
    value = parse_number(C{row, outCols});
else
    value = zeros(numel(outCols),1);
    for k = 1:numel(outCols)
        value(k) = parse_number(C{row, outCols(k)});
    end
end
end

function x = parse_number(v)
if isempty(v)
    x = NaN;
    return;
end
if isnumeric(v)
    if isempty(v) || all(isnan(v(:)))
        x = NaN;
    else
        x = v(1);
    end
    return;
end
if iscell(v)
    if isempty(v)
        x = NaN;
    else
        x = parse_number(v{1});
    end
    return;
end
try
    tfMissing = all(ismissing(v), 'all');
catch
    tfMissing = false;
end
if tfMissing
    x = NaN;
    return;
end
s = strtrim(string(v));
if numel(s) > 1
    s = s(1);
end
if ismissing(s) || strlength(s) == 0
    x = NaN;
    return;
end
s = replace(s, ',', '');
if endsWith(s, "%")
    s2 = extractBefore(s, strlength(s));
    x = str2double(s2) / 100;
else
    x = str2double(s);
end
if isnan(x)
    tok = regexp(char(s), '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match', 'once');
    if isempty(tok)
        x = NaN;
    else
        x = str2double(tok);
        if endsWith(s, "%")
            x = x / 100;
        end
    end
end
end

function out = safe_default(x, fallback)
if isnan(x)
    out = fallback;
else
    out = x;
end
end

function [logLambda, gamma] = calibrate_age_infection_model(pars, targetInfPrev, targetAgeOR)
obj = @(logGamma) log(expected_age_or_exact(solve_loglambda_for_gamma(exp(logGamma), pars, targetInfPrev), exp(logGamma), pars) / targetAgeOR);
grid = linspace(-3, 3, 49);
vals = nan(size(grid));
for i = 1:numel(grid)
    vals(i) = obj(grid(i));
end

bracketFound = false;
for i = 1:(numel(grid) - 1)
    if isfinite(vals(i)) && isfinite(vals(i+1)) && sign(vals(i)) ~= sign(vals(i+1))
        logGamma = fzero(obj, [grid(i), grid(i+1)]);
        bracketFound = true;
        break;
    end
end
if ~bracketFound
    logGamma = fminsearch(@(x) obj(x).^2, log(1.2), optimset('Display','off'));
end
gamma = exp(logGamma);
logLambda = solve_loglambda_for_gamma(gamma, pars, targetInfPrev);
end

function logLambda = solve_loglambda_for_gamma(gamma, pars, targetInfPrev)
obj = @(ll) expected_infection_prevalence_exact(ll, gamma, pars) - targetInfPrev;
lo = -30;
hi = 10;
fLo = obj(lo);
fHi = obj(hi);
while fLo > 0
    lo = lo - 5;
    fLo = obj(lo);
end
while fHi < 0
    hi = hi + 5;
    fHi = obj(hi);
end
logLambda = fzero(obj, [lo, hi]);
end

function prev = expected_infection_prevalence_exact(logLambda, gamma, pars)
prev = 0;
for iAge = 1:numel(pars.exactAgeValues)
    age = pars.exactAgeValues(iAge);
    pAge = pars.exactAgeProb(iAge);
    ageGroup = broad_age_group_from_years(age);
    pMJ = pars.mjPrevByAge(ageGroup);
    pContact = pars.contactPrevByAge(ageGroup);
    pRenal = pars.renalPrevByAge(ageGroup);
    baseCumHaz = age_cumulative_infection_hazard(age, logLambda, gamma);
    for mj = 0:1
        prMJ = bern_prob(mj, pMJ);
        for ct = 0:1
            prCT = bern_prob(ct, pContact);
            for rn = 0:1
                prRN = bern_prob(rn, pRenal);
                eta = log(baseCumHaz) ...
                    + log(pars.infOR.MJ) * mj ...
                    + log(pars.infOR.contact) * ct ...
                    + log(pars.infOR.renal) * rn;
                pInf = 1 - exp(-exp(eta));
                prev = prev + pAge * prMJ * prCT * prRN * pInf;
            end
        end
    end
end
end

function ageOR = expected_age_or_exact(logLambda, gamma, pars)
[pOlder, pYounger] = expected_infection_by_age_threshold(logLambda, gamma, pars, 25);
ageOR = odds_from_prob(pOlder) / max(odds_from_prob(pYounger), eps);
end

function [pOlder, pYounger] = expected_infection_by_age_threshold(logLambda, gamma, pars, ageCut)
olderMask = pars.exactAgeValues >= ageCut;
youngerMask = pars.exactAgeValues < ageCut;
pOlder = expected_infection_in_age_subset(logLambda, gamma, pars, olderMask);
pYounger = expected_infection_in_age_subset(logLambda, gamma, pars, youngerMask);
end

function prev = expected_infection_in_age_subset(logLambda, gamma, pars, ageMask)
ageVals = pars.exactAgeValues(ageMask);
ageProb = pars.exactAgeProb(ageMask);
ageProb = ageProb ./ sum(ageProb);
prev = 0;
for iAge = 1:numel(ageVals)
    age = ageVals(iAge);
    pAge = ageProb(iAge);
    ageGroup = broad_age_group_from_years(age);
    pMJ = pars.mjPrevByAge(ageGroup);
    pContact = pars.contactPrevByAge(ageGroup);
    pRenal = pars.renalPrevByAge(ageGroup);
    baseCumHaz = age_cumulative_infection_hazard(age, logLambda, gamma);
    for mj = 0:1
        prMJ = bern_prob(mj, pMJ);
        for ct = 0:1
            prCT = bern_prob(ct, pContact);
            for rn = 0:1
                prRN = bern_prob(rn, pRenal);
                eta = log(baseCumHaz) ...
                    + log(pars.infOR.MJ) * mj ...
                    + log(pars.infOR.contact) * ct ...
                    + log(pars.infOR.renal) * rn;
                pInf = 1 - exp(-exp(eta));
                prev = prev + pAge * prMJ * prCT * prRN * pInf;
            end
        end
    end
end
end

function lambdaEarly = calibrate_early_hazard(pars, logLambda, gamma, targetActive2y, screenWindow, earlyLateRatio)
obj = @(lam) expected_active_within_window(lam, pars, logLambda, gamma, screenWindow) - targetActive2y;
lambdaEarly = fzero(obj, [1e-8, 10]);
if earlyLateRatio < 1
    error('earlyLateRatio must be >= 1');
end
end

function prev = expected_active_within_window(lambdaEarly, pars, logLambda, gamma, screenWindow)
prev = 0;
for iAge = 1:numel(pars.exactAgeValues)
    age = pars.exactAgeValues(iAge);
    pAge = pars.exactAgeProb(iAge);
    ageGroup = broad_age_group_from_years(age);
    pMJ = pars.mjPrevByAge(ageGroup);
    pContact = pars.contactPrevByAge(ageGroup);
    pRenal = pars.renalPrevByAge(ageGroup);
    pDiab = pars.diabetesPrevByAge(ageGroup);
    pSmoke = pars.smokingPrevByAge(ageGroup);
    pCLD = pars.cldPrevByAge(ageGroup);
    pAlcohol = pars.alcoholPrevByAge(ageGroup);
    baseCumHaz = age_cumulative_infection_hazard(age, logLambda, gamma);
    for mj = 0:1
        prMJ = bern_prob(mj, pMJ);
        for ct = 0:1
            prCT = bern_prob(ct, pContact);
            for rn = 0:1
                prRN = bern_prob(rn, pRenal);
                eta = log(baseCumHaz) ...
                    + log(pars.infOR.MJ) * mj ...
                    + log(pars.infOR.contact) * ct ...
                    + log(pars.infOR.renal) * rn;
                pInf = 1 - exp(-exp(eta));
                for db = 0:1
                    prDB = bern_prob(db, pDiab);
                    for sm = 0:1
                        prSM = bern_prob(sm, pSmoke);
                        for cl = 0:1
                            prCL = bern_prob(cl, pCLD);
                            for al = 0:1
                                prAL = bern_prob(al, pAlcohol);
                                mult = (pars.disOR.MJ ^ mj) ...
                                     * (pars.disOR.contact ^ ct) ...
                                     * (pars.disOR.renal ^ rn) ...
                                     * (pars.disOR.diabetes ^ db) ...
                                     * (pars.disOR.smoking ^ sm) ...
                                     * (pars.disOR.cld ^ cl) ...
                                     * (pars.disOR.alcohol ^ al);
                                pAct = 1 - exp(-screenWindow * lambdaEarly * mult);
                                prev = prev + pAge * prMJ * prCT * prRN * prDB * prSM * prCL * prAL * pInf * pAct;
                            end
                        end
                    end
                end
            end
        end
    end
end
end

function [out, cohort] = simulate_one_cohort(N, pars, reg, logLambda, gamma, lambdaEarly, lambdaLate, opts)
ageYears = discrete_draw_values(pars.exactAgeValues, pars.exactAgeProb, N);
ageGroup = broad_age_group_from_years(ageYears);

female = rand(N,1) < pars.totalFemalePrev;
BCG = rand(N,1) < pars.totalBCGPrev;

pMJ = age_lookup(pars.mjPrevByAge, ageGroup);
pContact = age_lookup(pars.contactPrevByAge, ageGroup);
pRenal = age_lookup(pars.renalPrevByAge, ageGroup);
pDiabetes = age_lookup(pars.diabetesPrevByAge, ageGroup);
pSmoking = age_lookup(pars.smokingPrevByAge, ageGroup);
pCLD = age_lookup(pars.cldPrevByAge, ageGroup);
pAlcohol = age_lookup(pars.alcoholPrevByAge, ageGroup);

MJ = rand(N,1) < pMJ;
contact = rand(N,1) < pContact;
renal = rand(N,1) < pRenal;
diabetes = rand(N,1) < pDiabetes;
smoking = rand(N,1) < pSmoking;
cld = rand(N,1) < pCLD;
alcohol = rand(N,1) < pAlcohol;

baseCumHaz = age_cumulative_infection_hazard(ageYears, logLambda, gamma);
linPredInf = log(pars.infOR.MJ) .* MJ ...
    + log(pars.infOR.contact) .* contact ...
    + log(pars.infOR.renal) .* renal;
etaInf = log(baseCumHaz) + linPredInf;
pInf = 1 - exp(-exp(etaInf));
infected = rand(N,1) < pInf;

multDisease = (pars.disOR.MJ .^ MJ) ...
    .* (pars.disOR.contact .^ contact) ...
    .* (pars.disOR.renal .^ renal) ...
    .* (pars.disOR.diabetes .^ diabetes) ...
    .* (pars.disOR.smoking .^ smoking) ...
    .* (pars.disOR.cld .^ cld) ...
    .* (pars.disOR.alcohol .^ alcohol);

avgLatentAtScreenGivenInf = average_survival_to_random_screen(multDisease, lambdaEarly, opts.screenWindow);
preventableActiveRiskGivenInf = preventable_active_risk(multDisease, lambdaEarly, lambdaLate, opts.screenWindow, opts.followHorizon);
ltbiRiskScore = pInf;
cureTargetScore = pInf .* avgLatentAtScreenGivenInf;
preventTargetScore = pInf .* preventableActiveRiskGivenInf;

assert_nx1(N, ageGroup, 'ageGroup');
assert_nx1(N, pMJ, 'pMJ');
assert_nx1(N, MJ, 'MJ');
assert_nx1(N, pInf, 'pInf');
assert_nx1(N, infected, 'infected');
assert_nx1(N, multDisease, 'multDisease');
assert_nx1(N, ltbiRiskScore, 'ltbiRiskScore');
assert_nx1(N, cureTargetScore, 'cureTargetScore');
assert_nx1(N, preventTargetScore, 'preventTargetScore');

tActive = inf(N,1);
idxInf = find(infected);
if ~isempty(idxInf)
    rate1 = lambdaEarly .* multDisease(idxInf);
    t1 = exprnd_local(1 ./ rate1);
    early = t1 <= opts.screenWindow;
    tActive(idxInf(early)) = t1(early);
    idxLate = idxInf(~early);
    if ~isempty(idxLate)
        rate2 = lambdaLate .* multDisease(idxLate);
        t2 = exprnd_local(1 ./ rate2);
        tActive(idxLate) = opts.screenWindow + t2;
    end
end

[screened, screenPriorityScore, screenPriorityRank] = select_screened_people( ...
    ltbiRiskScore, cureTargetScore, preventTargetScore, opts.screenCoverage, opts.screeningStrategy);
screenTime = inf(N,1);
screenTime(screened) = rand(sum(screened),1) * opts.screenWindow;

assert_nx1(N, tActive, 'tActive');
assert_nx1(N, screened, 'screened');
assert_nx1(N, screenTime, 'screenTime');

activeAtScreen = screened & infected & (tActive <= screenTime);
latentAtScreen = screened & infected & (tActive > screenTime);

[testSensitivityUsed, testSpecificityUsed] = get_test_performance(BCG, opts);
testSpecificityNoBCGCounterfactual = get_counterfactual_noBCG_specificity(BCG, opts);
falsePositiveRate = 1 - testSpecificityUsed;
falsePositiveRateNoBCGCounterfactual = 1 - testSpecificityNoBCGCounterfactual;

testPositive = false(N,1);
falsePositiveCounterfactualNoBCG = false(N,1);
infScreened = screened & infected;
uninfScreened = screened & ~infected;
testPositive(infScreened) = rand(sum(infScreened),1) < testSensitivityUsed(infScreened);
if any(uninfScreened)
    uUninf = rand(sum(uninfScreened),1);
    testPositive(uninfScreened) = uUninf < falsePositiveRate(uninfScreened);
    falsePositiveCounterfactualNoBCG(uninfScreened) = uUninf < falsePositiveRateNoBCGCounterfactual(uninfScreened);
end

falsePositiveTest = screened & ~infected & testPositive;
falsePositiveTestBCG = falsePositiveTest & BCG;
falsePositiveTestNoBCG = falsePositiveTest & ~BCG;
falsePositiveDueToBCG = falsePositiveTest & ~falsePositiveCounterfactualNoBCG;
eligibleTPT = screened & testPositive & ~activeAtScreen;
startedTPT = eligibleTPT & (rand(N,1) < opts.pStartTPT);

adrStop = false(N,1);
completedTPT = false(N,1);
stoppedOther = false(N,1);

doseFraction = zeros(N,1);
courseStopTime = inf(N,1);
partialEffAssigned = zeros(N,1);
fullEffAssigned = zeros(N,1);

treatIdx = find(startedTPT);
if ~isempty(treatIdx)
    adrDraw = rand(numel(treatIdx),1) < reg.pADRstop;
    adrStop(treatIdx(adrDraw)) = true;

    remainIdx = treatIdx(~adrDraw);
    if (1 - reg.pADRstop) > 0
        pCompleteGivenNoADR = reg.pComplete / (1 - reg.pADRstop);
    else
        pCompleteGivenNoADR = 0;
    end
    pCompleteGivenNoADR = min(max(pCompleteGivenNoADR, 0), 1);

    if ~isempty(remainIdx)
        compDraw = rand(numel(remainIdx),1) < pCompleteGivenNoADR;
        completedTPT(remainIdx(compDraw)) = true;
        stoppedOther(remainIdx(~compDraw)) = true;
    end

    doseFraction(completedTPT) = 1;
    doseFraction(adrStop) = opts.partialDoseFractionADR;
    doseFraction(stoppedOther) = opts.partialDoseFractionOther;

    courseStopTime(startedTPT) = screenTime(startedTPT) + (reg.months / 12) .* doseFraction(startedTPT);
    fullEffAssigned(completedTPT) = reg.effFull;
    incomplete = startedTPT & ~completedTPT;
    if any(incomplete)
        partialEffAssigned(incomplete) = regimen_partial_efficacy(reg, doseFraction(incomplete));
    end
end

falsePositiveTreated = startedTPT & ~infected;
falsePositiveCompleted = completedTPT & ~infected;
falsePositiveTreatedBCG = falsePositiveTreated & BCG;
falsePositiveTreatedNoBCG = falsePositiveTreated & ~BCG;
falsePositiveCompletedBCG = falsePositiveCompleted & BCG;
falsePositiveCompletedNoBCG = falsePositiveCompleted & ~BCG;
excessCourseStartedDueToBCG = falsePositiveTreated & falsePositiveDueToBCG;
excessCourseCompletedDueToBCG = falsePositiveCompleted & falsePositiveDueToBCG;

protectedFull = completedTPT & infected & (tActive > courseStopTime) ...
    & (rand(N,1) < fullEffAssigned);
protectedPartial = (startedTPT & ~completedTPT) & infected & (tActive > courseStopTime) ...
    & (rand(N,1) < partialEffAssigned);
protectedAny = protectedFull | protectedPartial;

curedInfection = protectedAny;
curedInfectionFull = protectedFull;
curedInfectionPartial = protectedPartial;
preventedActiveTB = protectedAny & (tActive <= opts.followHorizon);
preventedActiveTBFull = protectedFull & (tActive <= opts.followHorizon);
preventedActiveTBPartial = protectedPartial & (tActive <= opts.followHorizon);

activeBy2y = infected & (tActive <= opts.screenWindow);
activeBy20y = infected & (tActive <= opts.followHorizon);

nScreened = sum(screened);
nScreenedBCG = sum(screened & BCG);
nScreenedNoBCG = sum(screened & ~BCG);
nUninfectedScreenedBCG = sum(uninfScreened & BCG);
nUninfectedScreenedNoBCG = sum(uninfScreened & ~BCG);
nStarted = sum(startedTPT);
nCompleted = sum(completedTPT);
nADRstop = sum(adrStop);
nStoppedOther = sum(stoppedOther);
nCured = sum(curedInfection);
nCuredFull = sum(curedInfectionFull);
nCuredPartial = sum(curedInfectionPartial);
nPrevented = sum(preventedActiveTB);
nPreventedFull = sum(preventedActiveTBFull);
nPreventedPartial = sum(preventedActiveTBPartial);
nFalsePosTests = sum(falsePositiveTest);
nFalsePosTestsBCG = sum(falsePositiveTestBCG);
nFalsePosTestsNoBCG = sum(falsePositiveTestNoBCG);
nFalsePosTreated = sum(falsePositiveTreated);
nFalsePosTreatedBCG = sum(falsePositiveTreatedBCG);
nFalsePosTreatedNoBCG = sum(falsePositiveTreatedNoBCG);
nFalsePosCompleted = sum(falsePositiveCompleted);
nFalsePosCompletedBCG = sum(falsePositiveCompletedBCG);
nFalsePosCompletedNoBCG = sum(falsePositiveCompletedNoBCG);
nExcessFalsePosDueToBCG = sum(falsePositiveDueToBCG);
nExcessCoursesStartedDueToBCG = sum(excessCourseStartedDueToBCG);
nExcessCoursesCompletedDueToBCG = sum(excessCourseCompletedDueToBCG);

out = struct();
out.testType = char(opts.testType);
out.screeningStrategy = char(opts.screeningStrategy);
out.regimen = reg.label;
out.nScreened = nScreened;
out.nInfected = sum(infected);
out.nLatentAtScreen = sum(latentAtScreen);
out.nActiveAtScreen = sum(activeAtScreen);
out.nTestPositive = sum(testPositive & screened);
out.nTestPositiveNonActive = sum(testPositive & screened & ~activeAtScreen);
out.nIGRApos = out.nTestPositiveNonActive; % backward compatibility name
out.nFalsePositiveTests = nFalsePosTests;
out.nFalsePositiveTestsBCG = nFalsePosTestsBCG;
out.nFalsePositiveTestsNoBCG = nFalsePosTestsNoBCG;
out.nFalsePositiveTreated = nFalsePosTreated;
out.nFalsePositiveTreatedBCG = nFalsePosTreatedBCG;
out.nFalsePositiveTreatedNoBCG = nFalsePosTreatedNoBCG;
out.nFalsePositiveCompleted = nFalsePosCompleted;
out.nFalsePositiveCompletedBCG = nFalsePosCompletedBCG;
out.nFalsePositiveCompletedNoBCG = nFalsePosCompletedNoBCG;
out.nExcessFalsePositiveTestsDueToBCG = nExcessFalsePosDueToBCG;
out.nExcessCoursesStartedDueToBCG = nExcessCoursesStartedDueToBCG;
out.nExcessCoursesCompletedDueToBCG = nExcessCoursesCompletedDueToBCG;
out.nExcessCoursesDueToBCG = nExcessCoursesStartedDueToBCG;
out.nScreenedBCG = nScreenedBCG;
out.nScreenedNoBCG = nScreenedNoBCG;
out.nUninfectedScreenedBCG = nUninfectedScreenedBCG;
out.nUninfectedScreenedNoBCG = nUninfectedScreenedNoBCG;
out.nStartTPT = nStarted;
out.nCompleteTPT = nCompleted;
out.nADRstop = nADRstop;
out.nStoppedOther = nStoppedOther;
out.nPartialCourses = nADRstop + nStoppedOther;
out.nTotalCoursesStarted = nStarted;
out.nTotalCoursesCompleted = nCompleted;
out.nCuredInfection = nCured;
out.nCuredInfectionFull = nCuredFull;
out.nCuredInfectionPartial = nCuredPartial;
out.nPreventedActiveTB = nPrevented;
out.nPreventedActiveTBFull = nPreventedFull;
out.nPreventedActiveTBPartial = nPreventedPartial;
out.nActiveBy2y = sum(activeBy2y);
out.nActiveBy20y = sum(activeBy20y);
out.completionRateObserved = safe_fraction(nCompleted, nStarted);
out.adrRateObserved = safe_fraction(nADRstop, nStarted);
out.partialCourseRateObserved = safe_fraction(nADRstop + nStoppedOther, nStarted);
out.propCoursesFalsePositive = safe_fraction(nFalsePosTreated, nStarted);
out.falsePositiveRateObservedBCG = safe_fraction(nFalsePosTestsBCG, nUninfectedScreenedBCG);
out.falsePositiveRateObservedNoBCG = safe_fraction(nFalsePosTestsNoBCG, nUninfectedScreenedNoBCG);
out.bcgAttributableFalsePositiveRateObserved = safe_fraction(nExcessFalsePosDueToBCG, nUninfectedScreenedBCG);
out.propFalsePositiveTestsDueToBCGAmongBCG = safe_fraction(nExcessFalsePosDueToBCG, nFalsePosTestsBCG);
out.propCoursesDueToBCGAmongFalsePositiveCoursesBCG = safe_fraction(nExcessCoursesStartedDueToBCG, nFalsePosTreatedBCG);
out.NNS_cureInfection = safe_divide(nScreened, nCured);
out.NNS_preventActiveTB = safe_divide(nScreened, nPrevented);
out.NNS_falsePositiveTreated = safe_divide(nScreened, nFalsePosTreated);
out.NNT_started_cureInfection = safe_divide(nStarted, nCured);
out.NNT_started_preventActiveTB = safe_divide(nStarted, nPrevented);

cohort = table((1:N)', ageGroup, ageYears, female, BCG, MJ, contact, renal, ...
    diabetes, smoking, cld, alcohol, pInf, infected, multDisease, ...
    ltbiRiskScore, cureTargetScore, preventTargetScore, ...
    screenPriorityScore, screenPriorityRank, tActive, ...
    screened, screenTime, activeAtScreen, latentAtScreen, ...
    testSensitivityUsed, testSpecificityUsed, testSpecificityNoBCGCounterfactual, ...
    testPositive, falsePositiveTest, falsePositiveCounterfactualNoBCG, falsePositiveDueToBCG, ...
    eligibleTPT, startedTPT, adrStop, stoppedOther, completedTPT, ...
    doseFraction, courseStopTime, fullEffAssigned, partialEffAssigned, ...
    falsePositiveTreated, falsePositiveTreatedBCG, falsePositiveCompleted, falsePositiveCompletedBCG, ...
    excessCourseStartedDueToBCG, excessCourseCompletedDueToBCG, ...
    curedInfection, curedInfectionFull, curedInfectionPartial, ...
    preventedActiveTB, preventedActiveTBFull, preventedActiveTBPartial, ...
    'VariableNames', {'id','ageGroup','ageYears','female','BCG', ...
    'MJ','contact','renal','diabetes','smoking','chronicLungDisease', ...
    'alcoholDrugs','pInfection','infected','diseaseMultiplier', ...
    'ltbiRiskScore','cureTargetScore','preventTargetScore', ...
    'screenPriorityScore','screenPriorityRank','tActiveUntreated', ...
    'screened','screenTime','activeAtScreen', ...
    'latentAtScreen','testSensitivityUsed','testSpecificityUsed', ...
    'testSpecificityNoBCGCounterfactual','testPositive', ...
    'falsePositiveTest','falsePositiveCounterfactualNoBCG','falsePositiveDueToBCG', ...
    'eligibleTPT','startedTPT','adrStop','stoppedOther', ...
    'completedTPT','doseFractionTaken','treatmentStopTime','fullEffAssigned', ...
    'partialEffAssigned','falsePositiveTreated','falsePositiveTreatedBCG', ...
    'falsePositiveCompleted','falsePositiveCompletedBCG', ...
    'excessCourseStartedDueToBCG','excessCourseCompletedDueToBCG', ...
    'curedInfection','curedInfectionFull','curedInfectionPartial', ...
    'preventedActiveTB','preventedActiveTBFull','preventedActiveTBPartial'});
end

function avgSurv = average_survival_to_random_screen(multDisease, lambdaEarly, screenWindow)
rate = lambdaEarly .* multDisease;
avgSurv = ones(size(rate));
idx = rate > 0;
avgSurv(idx) = (1 - exp(-rate(idx) .* screenWindow)) ./ max(rate(idx) .* screenWindow, eps);
end

function preventable = preventable_active_risk(multDisease, lambdaEarly, lambdaLate, screenWindow, followHorizon)
avgSurvToScreen = average_survival_to_random_screen(multDisease, lambdaEarly, screenWindow);
survToHorizon = exp(-(lambdaEarly .* multDisease .* screenWindow) ...
    - (lambdaLate .* multDisease .* max(followHorizon - screenWindow, 0)));
preventable = max(avgSurvToScreen - survToHorizon, 0);
end

function [screened, priorityScore, priorityRank] = select_screened_people(ltbiRiskScore, cureTargetScore, preventTargetScore, screenCoverage, strategy)
N = numel(ltbiRiskScore);
switch lower(strtrim(char(strategy)))
    case 'random'
        priorityScore = rand(N,1);
    case 'ltbi'
        priorityScore = ltbiRiskScore;
    case 'cure'
        priorityScore = cureTargetScore;
    case 'prevent'
        priorityScore = preventTargetScore;
    otherwise
        error('Unknown screeningStrategy: %s', char(strategy));
end

nToScreen = min(max(round(screenCoverage * N), 0), N);
priorityScore = priorityScore(:);
tieBreaker = rand(N,1) * 1e-12;
[~, order] = sort(priorityScore + tieBreaker, 'descend');
priorityRank = zeros(N,1);
priorityRank(order) = 1:N;
screened = false(N,1);
if nToScreen > 0
    screened(order(1:nToScreen)) = true;
end
end

function [sens, spec] = get_test_performance(BCG, opts)
BCG = logical(BCG(:));
sens = zeros(size(BCG));
spec = zeros(size(BCG));
if opts.testType == "IGRA"
    sens(:) = opts.testSensitivity;
    spec(:) = opts.testSpecificity;
else
    sens(:) = opts.tstSensitivity;
    spec(~BCG) = opts.tstSpecificityNoBCG;
    spec(BCG) = opts.tstSpecificityBCG;
end
end

function spec = get_counterfactual_noBCG_specificity(BCG, opts)
BCG = logical(BCG(:)); %#ok<NASGU>
spec = zeros(size(BCG));
if opts.testType == "IGRA"
    spec(:) = opts.testSpecificity;
else
    spec(:) = opts.tstSpecificityNoBCG;
end
end

function eff = regimen_partial_efficacy(reg, doseFrac)
frac = max(0, min(1, doseFrac));
mode = lower(strtrim(string(reg.partialMode)));
if isscalar(frac)
    frac = double(frac);
end
switch mode
    case "none"
        eff = zeros(size(frac));
    case "linear"
        eff = reg.effFull .* frac;
    case "threshold80"
        thr = reg.partialThreshold;
        if isnan(thr)
            thr = 0.80;
        end
        eff = reg.effFull .* double(frac >= thr);
    case "knots"
        eff = interp1(reg.partialDoseKnots, reg.partialEffKnots, frac, 'linear', 'extrap');
        eff = max(0, min(reg.effFull, eff));
    otherwise
        error('Unknown regimen partialMode: %s', reg.partialMode);
end
end

function v = age_lookup(vec, ageGroup)
vec = vec(:);
v = vec(ageGroup);
v = v(:);
end

function assert_nx1(N, x, name)
sz = size(x);
if ~isequal(sz, [N 1])
    error('Variable %s has size [%d %d], expected [%d 1].', name, sz(1), sz(2), N);
end
end

function idx = discrete_draw(prob, N)
prob = prob(:);
prob = prob ./ sum(prob);
edges = [0; cumsum(prob)];
u = rand(N,1);
idx = zeros(N,1);
for k = 1:numel(prob)
    idx(u > edges(k) & u <= edges(k+1)) = k;
end
idx(idx == 0) = numel(prob);
end

function vals = discrete_draw_values(valueVec, probVec, N)
idx = discrete_draw(probVec, N);
valueVec = valueVec(:);
vals = valueVec(idx);
vals = vals(:);
end

function pars = load_age_distribution(pars, ageFile, ageTableOverride, sheetName, age85PlusMax)
if nargin >= 3 && ~isempty(ageTableOverride)
    T = ageTableOverride;
    pars.ageDistributionSource = 'provided_table';
elseif isempty(ageFile) || ~isfile(ageFile)
    pars.exactAgeValues = (0:89)';
    pars.exactAgeProb = build_fallback_exact_age_prob(pars.exactAgeValues, pars.popFrac);
    pars.ageDistributionSource = '';
    pars.ageDistributionTable = table();
    return;
else
    [~, ~, ext] = fileparts(ageFile);
    if strcmpi(ext, '.csv')
        T = readtable(ageFile);
    else
        T = readtable(ageFile, 'Sheet', sheetName);
    end
    pars.ageDistributionSource = ageFile;
end
varNames = string(T.Properties.VariableNames);
lowerNames = lower(varNames);

ageCol = find(contains(lowerNames, 'age'), 1, 'first');
if isempty(ageCol)
    ageCol = 1;
end
weightCol = find(contains(lowerNames, 'smoothed'), 1, 'first');
if isempty(weightCol)
    weightCol = find(contains(lowerNames, 'proportion'), 1, 'first');
end
if isempty(weightCol)
    error('Could not identify a proportion column in the age distribution file: %s', ageFile);
end

ageColName = T.Properties.VariableNames{ageCol};
weightColName = T.Properties.VariableNames{weightCol};
nRows = height(T);
exactAges = [];
exactProb = [];
for i = 1:nRows
    bandVal = T{i, ageColName};
    if iscell(bandVal)
        bandVal = bandVal{1};
    end
    weightVal = T{i, weightColName};
    if iscell(weightVal)
        weightVal = weightVal{1};
    end
    weightVal = parse_number(weightVal);
    if isnan(weightVal)
        continue;
    end
    [lo, hi] = parse_age_band(bandVal, age85PlusMax);
    ages = (lo:hi)';
    exactAges = [exactAges; ages]; %#ok<AGROW>
    exactProb = [exactProb; repmat(weightVal / numel(ages), numel(ages), 1)]; %#ok<AGROW>
end

exactProb = exactProb ./ sum(exactProb);
pars.exactAgeValues = exactAges(:);
pars.exactAgeProb = exactProb(:);
pars.ageDistributionTable = T;
end

function ageFile = resolve_age_distribution_file(ageDistributionFile, csvFile)
ageDistributionFile = strtrim(string(ageDistributionFile));
if strlength(ageDistributionFile) > 0
    ageFile = char(ageDistributionFile);
    return;
end
csvDir = fileparts(char(csvFile));
candidateCsv = fullfile(csvDir, 'default_age_distribution.csv');
candidateXlsx = fullfile(csvDir, 'default_age_distribution.xlsx');
legacyXlsx = fullfile(csvDir, 'age_groups_SA.xlsx');
if isfile(candidateCsv)
    ageFile = candidateCsv;
elseif isfile(candidateXlsx)
    ageFile = candidateXlsx;
elseif isfile(legacyXlsx)
    ageFile = legacyXlsx;
else
    ageFile = '';
end
end

function prob = build_fallback_exact_age_prob(exactAges, popFrac3)
prob = zeros(size(exactAges));
prob(exactAges <= 4) = popFrac3(1) / 5;
prob(exactAges >= 5 & exactAges <= 14) = popFrac3(2) / 10;
prob(exactAges >= 15) = popFrac3(3) / sum(exactAges >= 15);
prob = prob ./ sum(prob);
end

function [lo, hi] = parse_age_band(ageBandValue, age85PlusMax)
s = strtrim(string(ageBandValue));
s = replace(s, '–', '-');
s = replace(s, '—', '-');
s = replace(s, '>=', '');
s = replace(s, '≥', '');
if contains(s, '+')
    lo = str2double(erase(s, '+'));
    hi = age85PlusMax;
elseif contains(s, '-')
    toks = split(s, '-');
    lo = str2double(toks(1));
    hi = str2double(toks(2));
else
    lo = str2double(s);
    hi = lo;
end
if isnan(lo) || isnan(hi)
    error('Could not parse age band value: %s', char(s));
end
end

function ageGroup = broad_age_group_from_years(ageYears)
ageYears = ageYears(:);
ageGroup = zeros(size(ageYears));
ageGroup(ageYears <= 4) = 1;
ageGroup(ageYears >= 5 & ageYears <= 14) = 2;
ageGroup(ageYears >= 15) = 3;
end

function H = age_cumulative_infection_hazard(ageYears, logLambda, gamma)
H = exp(logLambda) .* ((double(ageYears(:)) + 0.5) .^ gamma);
H = max(H, eps);
end

function t = exprnd_local(mu)
t = -mu .* log(rand(size(mu)));
end

function p = bern_prob(x, px)
if x == 1
    p = px;
else
    p = 1 - px;
end
end

function odds = odds_from_prob(p)
odds = p ./ max(1 - p, eps);
end

function x = safe_divide(a, b)
if b <= 0
    x = Inf;
else
    x = a / b;
end
end

function x = safe_fraction(a, b)
if b <= 0
    x = NaN;
else
    x = a / b;
end
end

function summary = summarise_results(raw)
vars = raw.Properties.VariableNames;
keep = false(numel(vars), 1);
for i = 1:numel(vars)
    keep(i) = isnumeric(raw.(vars{i})) || islogical(raw.(vars{i}));
end
vars = vars(keep);
Metric = strings(numel(vars), 1);
Median = zeros(numel(vars), 1);
Low95 = zeros(numel(vars), 1);
High95 = zeros(numel(vars), 1);
for i = 1:numel(vars)
    x = raw.(vars{i});
    x = double(x(~isnan(x)));
    Metric(i) = string(vars{i});
    if isempty(x)
        Median(i) = NaN;
        Low95(i) = NaN;
        High95(i) = NaN;
    else
        x = sort(x(:));
        Median(i) = median(x);
        Low95(i) = empirical_quantile(x, 0.025);
        High95(i) = empirical_quantile(x, 0.975);
    end
end
summary = table(Metric, Median, Low95, High95);
end

function q = empirical_quantile(x, p)
if isempty(x)
    q = NaN;
    return;
end
n = numel(x);
pos = 1 + (n - 1) * p;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    q = x(lo);
else
    q = x(lo) + (pos - lo) * (x(hi) - x(lo));
end
end
