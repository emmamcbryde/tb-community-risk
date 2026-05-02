function econConfig = build_economics_preset_kwab150_v9()
%BUILD_ECONOMICS_PRESET_KWAB150_V9 Return KWAB150 Australia 2019 AUD costs.
%
% Values are populated from the repository data/costs.csv cost parameter
% names used for the KWAB150 migrant LTBI cost-effectiveness model.

econConfig = build_default_economics_config_v9();
econConfig.metadata.currencyCode = 'AUD';
econConfig.metadata.priceYear = 2019;
econConfig.metadata.locationLabel = 'Australia';
econConfig.metadata.programCostBasis = 'total';
econConfig.metadata.sourceNotes = ['KWAB150 preset populated from local data/costs.csv: ' ...
    'cscreenqft, cscreentst, ctreat*, and ctb mid values. ' ...
    'Program setup/running and false-positive incremental costs are not specified in this preset.'];

econConfig.costs.test.IGRA = 113.48; % cscreenqft
econConfig.costs.test.TST = 116.07; % cscreentst

econConfig.costs.regimen.x3HP = 165.5072; % ctreat3HP
econConfig.costs.regimen.x4R = 123.3172; % ctreat4R
econConfig.costs.regimen.x3HR = 134.2272; % ctreat3HR
econConfig.costs.regimen.x6H = 187.7508; % ctreat6H
econConfig.costs.regimen.x9H = 254.8544; % ctreat9H

econConfig.costs.activeTBDiseasePerCase = 19079.6; % ctb
end
