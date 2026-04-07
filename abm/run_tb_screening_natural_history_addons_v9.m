function out = run_tb_screening_natural_history_addons_v9(resultsOrCsv, Nattr)
% Convenience wrapper for the natural-history / attributable-risk add-ons.
%
% If a v9 RESULTS struct is supplied, the do-nothing summary is derived from
% that run and the attributable-risk table is estimated from a larger
% no-intervention cohort using the same settings.

if nargin < 2 || isempty(Nattr)
    Nattr = 50000;
end

if nargin < 1 || isempty(resultsOrCsv)
    results = run_tb_screening_example_v9;
else
    results = resultsOrCsv;
    if ~isstruct(resultsOrCsv)
        results = tb_screening_mc_model_v9(resultsOrCsv, 'nReps', 2000);
    end
end

out = struct();
out.doNothing = run_tb_screening_do_nothing_v9(results);
out.attributable = run_tb_screening_reactivation_attributable_v9(results, Nattr);
end
