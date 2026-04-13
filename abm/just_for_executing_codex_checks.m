clear functions
rehash
appState = app_startup_v9();

disp(isfield(appState, 'Schema'))
disp(isfield(appState, 'CurrentConfig'))
disp(isfield(appState, 'LastValidationReport'))
disp(appState.ResultsAreStale)
disp(appState.LastRunSucceeded)


%% test 2

appState = app_startup_v9();
appState = app_validate_v9(appState);

disp(~isempty(appState.LastValidationReport))
disp(appState.ResultsAreStale)

%% test 3

appState = app_startup_v9();
appState.CurrentUIState.N = 200;
appState.CurrentUIState.nReps = 10;

appState = app_run_v9(appState);

disp(appState.LastRunSucceeded)
disp(appState.ResultsAreStale)
disp(~isempty(appState.LastResults))
disp(~isempty(appState.LastBundle))

%% test 4

appState.CurrentUIState.screenCoverage = 0.40;
[appState, report] = app_input_changed_v9(appState, appState.CurrentUIState);

disp(appState.LastRunSucceeded)
disp(appState.ResultsAreStale)
disp(~isempty(appState.LastResults))
disp(~isempty(appState.LastBundle))
disp(isempty(appState.LastExports))
disp(isempty(appState.LastDoNothing))
disp(isempty(appState.LastAttributableRisk))

%% test 5
clear functions
rehash
appState = app_run_do_nothing_v9(appState);
disp(~isempty(appState.LastDoNothing))
disp(appState.LastBundle.doNothing.available)

%%
clear functions
rehash

schema = build_ui_schema_v9();
cfg = build_default_config_v9();

uiState = config_to_ui_state_v9(cfg, schema);
cfg2 = ui_state_to_config_v9(uiState, schema, cfg);

disp(cfg.targetAgeOR)
disp(cfg2.targetAgeOR)
disp(strcmp(cfg.testType, cfg2.testType))
disp(strcmp(cfg.regimen, cfg2.regimen))
disp(strcmp(cfg.screeningStrategy, cfg2.screeningStrategy))

disp(schema_field_key_v9('diseaseOR.alcohol'))
disp(schema_component_name_v9('screenCoverage', 'number'))
disp(schema_component_name_v9('testType', 'select'))
disp(schema_component_name_v9('ageDistributionTable', 'table'))
disp(schema_component_name_v9('useDefaults', 'boolean'))