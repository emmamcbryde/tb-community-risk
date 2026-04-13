function saveInfo = save_scenario_v9(config, filename)
%SAVE_SCENARIO_V9 Save APY v9 scenario config to JSON.

if nargin < 2 || isempty(filename)
    error('save_scenario_v9 requires a filename.');
end

payload = struct();
payload.contractVersion = 'apy_v9_scenario_v1';
payload.modelVersion = 'v9';
payload.savedAt = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
payload.scenarioLabel = get_scenario_label(config);
payload.config = config_to_json_struct(config);

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
end

function label = get_scenario_label(config)
if isstruct(config) && isfield(config, 'scenarioLabel') && ~isempty(config.scenarioLabel)
    label = char(string(config.scenarioLabel));
else
    label = 'APY scenario';
end
end

function out = config_to_json_struct(config)
out = config;
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
