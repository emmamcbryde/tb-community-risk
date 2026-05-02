function saveInfo = save_scenario_v9(config, filename, economicsConfig)
%SAVE_SCENARIO_V9 Save APY v9 scenario config to JSON.

if nargin < 2 || isempty(filename)
    error('save_scenario_v9 requires a filename.');
end
if nargin < 3
    economicsConfig = [];
end

payload = struct();
payload.contractVersion = 'apy_v9_scenario_v1';
payload.modelVersion = 'v9';
payload.savedAt = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
payload.scenarioLabel = get_scenario_label(config);
payload.config = config_to_json_struct(config);
if ~economics_config_is_blank(economicsConfig)
    payload.economics = economics_config_to_json_struct(economicsConfig);
end

jsonText = jsonencode(payload, 'PrettyPrint', true);
fid = fopen(filename, 'w');
if fid < 0
    error('Could not open scenario file for writing: %s', filename);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonText);

saveInfo = struct();
saveInfo.filename = filename;
saveInfo.contractVersion = payload.contractVersion;
saveInfo.scenarioLabel = payload.scenarioLabel;
saveInfo.hasEconomics = isfield(payload, 'economics');
end

function label = get_scenario_label(config)
if isstruct(config) && isfield(config, 'scenarioLabel') && ~isempty(config.scenarioLabel)
    label = char(string(config.scenarioLabel));
else
    label = 'APY scenario';
end
end

function out = economics_config_to_json_struct(economicsConfig)
out = economicsConfig;
end

function tf = economics_config_is_blank(economicsConfig)
tf = true;
if isempty(economicsConfig) || ~isstruct(economicsConfig)
    return;
end
tf = struct_is_blank(economicsConfig);
end

function tf = struct_is_blank(s)
tf = true;
fields = fieldnames(s);
for i = 1:numel(fields)
    value = s.(fields{i});
    if isstruct(value)
        if ~struct_is_blank(value)
            tf = false;
            return;
        end
    elseif ~isempty(value)
        tf = false;
        return;
    end
end
end

function out = config_to_json_struct(config)
out = config;
if isfield(out, 'usesDefaults')
    if ~isfield(out, 'useDefaults') || isempty(out.useDefaults)
        out.useDefaults = out.usesDefaults;
    end
    out = rmfield(out, 'usesDefaults');
end
if isfield(out, 'ageDistributionTable') && istable(out.ageDistributionTable)
    out.ageDistributionTableRows = table2struct(out.ageDistributionTable, 'ToScalar', false);
    out = rmfield(out, 'ageDistributionTable');
else
    out.ageDistributionTableRows = struct.empty(0, 1);
    if isfield(out, 'ageDistributionTable')
        out = rmfield(out, 'ageDistributionTable');
    end
end
end
