function econ = run_health_economics_v9(results, econConfig)
%RUN_HEALTH_ECONOMICS_V9 Calculate first-pass economics from v9 results.

if nargin < 2 || isempty(econConfig)
    econConfig = build_default_economics_config_v9();
end

validationReport = validate_economics_config_v9(econConfig);
if ~validationReport.isValid
    error('run_health_economics_v9:InvalidEconomicsConfig', ...
        'Economics config contains fatal validation errors.');
end

validate_results(results);

econ = init_economics_result(econConfig);
econ.status.validationReport = validationReport;

testType = get_config_text(results, 'testType');
regimen = get_config_text(results, 'regimen');
econ.strategy.testType = testType;
econ.strategy.regimen = regimen;

econ.quantities.nScreened = summary_metric(results, 'nScreened');
econ.quantities.nTotalCoursesStarted = summary_metric(results, 'nTotalCoursesStarted');
econ.quantities.nFalsePositiveTreated = summary_metric(results, 'nFalsePositiveTreated');
econ.quantities.nCuredInfection = summary_metric(results, 'nCuredInfection');
econ.quantities.nPreventedActiveTB = summary_metric(results, 'nPreventedActiveTB');

[baselineCases, interventionCases] = baseline_intervention_cases(results);
econ.quantities.baselineActiveTBCases = baselineCases;
econ.quantities.interventionActiveTBCases = interventionCases;

[econ.unitCosts.testPerPerson, econ.status] = selected_test_cost(econConfig, testType, econ.status);
[econ.unitCosts.treatmentPerStarted, econ.status] = selected_regimen_cost(econConfig, regimen, econ.status);
econ.unitCosts.falsePositiveIncrementalPerPerson = get_optional_cost(econConfig.costs, ...
    'falsePositiveIncrementalPerPerson');
econ.unitCosts.activeTBDiseasePerCase = get_optional_cost(econConfig.costs, 'activeTBDiseasePerCase');

econ.costs.testingCost = multiply_if_available(econ.quantities.nScreened, econ.unitCosts.testPerPerson);
econ.status = note_missing_if_empty(econ.status, econ.costs.testingCost, ...
    'testingCost', 'Testing cost not calculated because nScreened or selected test cost is missing.');

econ.costs.treatmentCost = multiply_if_available(econ.quantities.nTotalCoursesStarted, ...
    econ.unitCosts.treatmentPerStarted);
econ.status = note_missing_if_empty(econ.status, econ.costs.treatmentCost, ...
    'treatmentCost', 'Treatment cost not calculated because starts or selected regimen cost is missing.');

econ.costs.falsePositiveIncrementalCost = multiply_if_available(econ.quantities.nFalsePositiveTreated, ...
    econ.unitCosts.falsePositiveIncrementalPerPerson);
if isempty(econ.unitCosts.falsePositiveIncrementalPerPerson)
    econ.status.messages{end+1, 1} = 'No false-positive incremental unit cost supplied; no extra false-positive cost is included.';
end

econ.costs.programSetupCost = get_optional_cost(econConfig.costs, 'programSetupTotal');
econ.costs.programRunningCost = get_optional_cost(econConfig.costs, 'programRunningTotal');

econ.costs.tbDiseaseCostsAverted = multiply_if_available(econ.quantities.nPreventedActiveTB, ...
    econ.unitCosts.activeTBDiseasePerCase);
econ.status = note_missing_if_empty(econ.status, econ.costs.tbDiseaseCostsAverted, ...
    'tbDiseaseCostsAverted', 'TB disease costs averted not calculated because prevented cases or disease cost is missing.');

econ.costs.baselineTBDiseaseCost = multiply_if_available(econ.quantities.baselineActiveTBCases, ...
    econ.unitCosts.activeTBDiseasePerCase);
econ.costs.interventionTBDiseaseCost = multiply_if_available(econ.quantities.interventionActiveTBCases, ...
    econ.unitCosts.activeTBDiseasePerCase);
if isempty(econ.costs.baselineTBDiseaseCost) || isempty(econ.costs.interventionTBDiseaseCost)
    econ.status.partialCalculations{end+1, 1} = 'Baseline/intervention TB disease costs were not calculated because case counts are unavailable.';
end

econ.costs.totalProgramCost = sum_required_costs({ ...
    econ.costs.testingCost, ...
    econ.costs.treatmentCost, ...
    zero_if_empty(econ.costs.falsePositiveIncrementalCost), ...
    econ.costs.programSetupCost, ...
    econ.costs.programRunningCost});
econ.status = note_missing_if_empty(econ.status, econ.costs.totalProgramCost, ...
    'totalProgramCost', 'Total program cost not calculated because one or more required cost components is missing.');

econ.costs.netCostVsBaseline = subtract_if_available(econ.costs.totalProgramCost, ...
    econ.costs.tbDiseaseCostsAverted);
econ.status = note_missing_if_empty(econ.status, econ.costs.netCostVsBaseline, ...
    'netCostVsBaseline', 'Net cost versus baseline not calculated because total program cost or TB costs averted is missing.');

[econ.costEffectiveness.costPerInfectionCured, econ.status] = divide_if_available( ...
    econ.costs.netCostVsBaseline, econ.quantities.nCuredInfection, ...
    'costPerInfectionCured', econ.status);
[econ.costEffectiveness.costPerTBCasePrevented, econ.status] = divide_if_available( ...
    econ.costs.netCostVsBaseline, econ.quantities.nPreventedActiveTB, ...
    'costPerTBCasePrevented', econ.status);

econ.summaryTable = build_summary_table(econ);
econ.summaryRows = table2struct(econ.summaryTable, 'ToScalar', false);
econ.status.isComplete = isempty(econ.status.missingInputs) && isempty(econ.status.notCalculated);
end

function validate_results(results)
if ~isstruct(results) || ~isfield(results, 'summary') || ~istable(results.summary)
    error('run_health_economics_v9:InvalidResults', ...
        'Input must be a v9 results struct containing results.summary.');
end
if ~ismember('Metric', results.summary.Properties.VariableNames) || ...
        ~ismember('Median', results.summary.Properties.VariableNames)
    error('run_health_economics_v9:InvalidResultsSummary', ...
        'results.summary must contain Metric and Median variables.');
end
end

function econ = init_economics_result(econConfig)
econ = struct();
econ.available = true;
econ.source = 'run_health_economics_v9';
econ.contractVersion = 'apy_v9_economics_results_v1';
econ.metadata = econConfig.metadata;
econ.inputs = econConfig;
econ.strategy = struct('testType', '', 'regimen', '');
econ.quantities = empty_quantity_struct();
econ.unitCosts = empty_unit_cost_struct();
econ.costs = empty_cost_struct();
econ.costEffectiveness = struct('costPerInfectionCured', [], 'costPerTBCasePrevented', []);
econ.summaryTable = table();
econ.summaryRows = struct.empty(0, 1);
econ.status = struct();
econ.status.isComplete = false;
econ.status.missingInputs = {};
econ.status.partialCalculations = {};
econ.status.notCalculated = {};
econ.status.messages = {};
econ.status.validationReport = [];
end

function s = empty_quantity_struct()
s = struct();
s.nScreened = [];
s.nTotalCoursesStarted = [];
s.nFalsePositiveTreated = [];
s.nCuredInfection = [];
s.nPreventedActiveTB = [];
s.baselineActiveTBCases = [];
s.interventionActiveTBCases = [];
end

function s = empty_unit_cost_struct()
s = struct();
s.testPerPerson = [];
s.treatmentPerStarted = [];
s.falsePositiveIncrementalPerPerson = [];
s.activeTBDiseasePerCase = [];
end

function s = empty_cost_struct()
s = struct();
s.testingCost = [];
s.treatmentCost = [];
s.falsePositiveIncrementalCost = [];
s.programSetupCost = [];
s.programRunningCost = [];
s.baselineTBDiseaseCost = [];
s.interventionTBDiseaseCost = [];
s.tbDiseaseCostsAverted = [];
s.totalProgramCost = [];
s.netCostVsBaseline = [];
end

function value = summary_metric(results, metricName)
value = [];
metricText = string(results.summary.Metric);
idx = find(metricText == string(metricName), 1, 'first');
if isempty(idx)
    return;
end
value = results.summary.Median(idx);
end

function textValue = get_config_text(results, fieldName)
textValue = '';
if isfield(results, 'interfaceConfig') && isstruct(results.interfaceConfig) ...
        && isfield(results.interfaceConfig, fieldName) && ~isempty(results.interfaceConfig.(fieldName))
    textValue = char(string(results.interfaceConfig.(fieldName)));
end
end

function [value, status] = selected_test_cost(econConfig, testType, status)
value = [];
if isempty(testType)
    status.missingInputs{end+1, 1} = 'results.interfaceConfig.testType';
    return;
end
if isfield(econConfig.costs.test, testType)
    value = econConfig.costs.test.(testType);
end
if isempty(value)
    status.missingInputs{end+1, 1} = ['costs.test.' testType];
end
end

function [value, status] = selected_regimen_cost(econConfig, regimen, status)
value = [];
if isempty(regimen)
    status.missingInputs{end+1, 1} = 'results.interfaceConfig.regimen';
    return;
end
try
    fieldName = regimen_cost_field_v9(regimen);
catch ME
    status.missingInputs{end+1, 1} = ME.message;
    return;
end
if isfield(econConfig.costs.regimen, fieldName)
    value = econConfig.costs.regimen.(fieldName);
end
if isempty(value)
    status.missingInputs{end+1, 1} = ['costs.regimen.' fieldName];
end
end

function value = get_optional_cost(parent, fieldName)
value = [];
if isstruct(parent) && isfield(parent, fieldName)
    value = parent.(fieldName);
end
end

function value = multiply_if_available(a, b)
if is_missing_number(a) || is_missing_number(b)
    value = [];
else
    value = a * b;
end
end

function value = subtract_if_available(a, b)
if is_missing_number(a) || is_missing_number(b)
    value = [];
else
    value = a - b;
end
end

function value = sum_required_costs(values)
value = 0;
for i = 1:numel(values)
    if is_missing_number(values{i})
        value = [];
        return;
    end
    value = value + values{i};
end
end

function value = zero_if_empty(value)
if isempty(value)
    value = 0;
end
end

function tf = is_missing_number(value)
tf = isempty(value) || ~(isnumeric(value) && isscalar(value)) || isnan(value);
end

function status = note_missing_if_empty(status, value, name, message)
if isempty(value)
    status.notCalculated{end+1, 1} = name;
    status.messages{end+1, 1} = message;
end
end

function [value, status] = divide_if_available(num, den, name, status)
if is_missing_number(num) || is_missing_number(den) || den == 0
    value = [];
    status.notCalculated{end+1, 1} = name;
    status.messages{end+1, 1} = [name ' not calculated because the denominator is missing, zero, or NaN.'];
else
    value = num / den;
end
end

function [baselineCases, interventionCases] = baseline_intervention_cases(results)
baselineCases = [];
interventionCases = [];
if ~isfield(results, 'raw') || ~istable(results.raw)
    return;
end
vars = results.raw.Properties.VariableNames;
if ~ismember('nActiveBy20y', vars) || ~ismember('nPreventedActiveTB', vars)
    return;
end
baselineVec = results.raw.nActiveBy20y;
interventionVec = results.raw.nActiveBy20y - results.raw.nPreventedActiveTB;
baselineCases = median_omit_nan(baselineVec);
interventionCases = median_omit_nan(interventionVec);
end

function value = median_omit_nan(x)
x = double(x(:));
x = x(~isnan(x));
if isempty(x)
    value = [];
else
    value = median(x);
end
end

function summaryTable = build_summary_table(econ)
Metric = [ ...
    "testingCost"; ...
    "treatmentCost"; ...
    "falsePositiveIncrementalCost"; ...
    "programSetupCost"; ...
    "programRunningCost"; ...
    "baselineTBDiseaseCost"; ...
    "interventionTBDiseaseCost"; ...
    "tbDiseaseCostsAverted"; ...
    "totalProgramCost"; ...
    "netCostVsBaseline"; ...
    "costPerInfectionCured"; ...
    "costPerTBCasePrevented" ...
    ];

Value = [ ...
    numeric_or_nan(econ.costs.testingCost); ...
    numeric_or_nan(econ.costs.treatmentCost); ...
    numeric_or_nan(econ.costs.falsePositiveIncrementalCost); ...
    numeric_or_nan(econ.costs.programSetupCost); ...
    numeric_or_nan(econ.costs.programRunningCost); ...
    numeric_or_nan(econ.costs.baselineTBDiseaseCost); ...
    numeric_or_nan(econ.costs.interventionTBDiseaseCost); ...
    numeric_or_nan(econ.costs.tbDiseaseCostsAverted); ...
    numeric_or_nan(econ.costs.totalProgramCost); ...
    numeric_or_nan(econ.costs.netCostVsBaseline); ...
    numeric_or_nan(econ.costEffectiveness.costPerInfectionCured); ...
    numeric_or_nan(econ.costEffectiveness.costPerTBCasePrevented) ...
    ];

Calculated = ~isnan(Value);
CurrencyCode = repmat(string(econ.metadata.currencyCode), numel(Metric), 1);
PriceYear = repmat(numeric_or_nan(econ.metadata.priceYear), numel(Metric), 1);
summaryTable = table(Metric, Value, Calculated, CurrencyCode, PriceYear);
end

function value = numeric_or_nan(value)
if is_missing_number(value)
    value = NaN;
end
end
