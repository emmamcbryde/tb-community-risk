function econConfig = build_default_economics_config_v9()
%BUILD_DEFAULT_ECONOMICS_CONFIG_V9 Return a blank v9 economics config.

econConfig = struct();

econConfig.metadata = struct();
econConfig.metadata.currencyCode = '';
econConfig.metadata.priceYear = [];
econConfig.metadata.locationLabel = '';
econConfig.metadata.sourceNotes = '';
econConfig.metadata.programCostBasis = '';

econConfig.costs = struct();

econConfig.costs.test = struct();
econConfig.costs.test.IGRA = [];
econConfig.costs.test.TST = [];

econConfig.costs.regimen = struct();
econConfig.costs.regimen.x3HP = [];
econConfig.costs.regimen.x4R = [];
econConfig.costs.regimen.x3HR = [];
econConfig.costs.regimen.x6H = [];
econConfig.costs.regimen.x9H = [];

econConfig.costs.falsePositiveIncrementalPerPerson = [];
econConfig.costs.activeTBDiseasePerCase = [];
econConfig.costs.programSetupTotal = [];
econConfig.costs.programRunningTotal = [];
end
