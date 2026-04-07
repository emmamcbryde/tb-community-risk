clear functions
rehash

results = run_tb_screening_user_options_v9;

dn   = run_tb_screening_do_nothing_v9(results);
attr = run_tb_screening_reactivation_attributable_v9(results, 50000);
out  = run_tb_screening_natural_history_addons_v9(results, 50000);
%%
comparison = run_tb_screening_compare_strategies_v9;
%%
gradient   = run_tb_screening_targeted_gradient_v9;
profile    = run_tb_screening_targeting_profile_v9;
optima     = run_tb_screening_targeting_optima_v9;
charts     = run_tb_screening_igra_charts_v9;
addons     = run_tb_screening_igra_addons_v9;


results = run_tb_screening_user_options_v9;
dn = run_tb_screening_do_nothing_v9(results);
attr = run_tb_screening_reactivation_attributable_v9(results, 50000);

