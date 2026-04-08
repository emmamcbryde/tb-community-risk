function out = run_tb_screening_do_nothing_v9(results)
% Summarise the natural-history ("do nothing") scenario from a v9 run.
%
% If RESULTS is omitted, the default v9 example is run.
%
% Outputs:
%   out.summary         table with medians and 95% intervals
%   out.derived         replicate-level derived metrics
%   out.resultsUsed     original results struct
%
% A CSV file is also written to the current folder:
%   tb_do_nothing_summary_v9.csv
outDir = get_output_dir_v9();
if nargin < 1 || isempty(results)
    results = run_tb_screening_example_v9;
end

validate_v9_results(results);
raw = results.raw;

postActive20 = raw.nActiveBy20y - raw.nPreventedActiveTB;
relReduction20 = safe_fraction_vec(raw.nPreventedActiveTB, raw.nActiveBy20y);
ltbiPrev = safe_fraction_vec(raw.nInfected, results.settings.N);
activePrev2y = safe_fraction_vec(raw.nActiveBy2y, results.settings.N);
activePrev20y = safe_fraction_vec(raw.nActiveBy20y, results.settings.N);
postActivePrev20y = safe_fraction_vec(postActive20, results.settings.N);

Metric = [ ...
    "LTBI prevalence at baseline (do nothing)"; ...
    "Active TB cases by 2 years (do nothing)"; ...
    "Active TB prevalence by 2 years (do nothing)"; ...
    "Active TB cases by 20 years (do nothing)"; ...
    "Active TB prevalence by 20 years (do nothing)"; ...
    "Active TB cases by 20 years after current strategy"; ...
    "Active TB prevalence by 20 years after current strategy"; ...
    "Active TB cases prevented by 20 years"; ...
    "Relative reduction in 20-year active TB burden"; ...
    "NNS to prevent one active TB case (current strategy)"; ...
    "NNS to cure one infection (current strategy)" ...
    ];

X = { ...
    ltbiPrev, ...
    raw.nActiveBy2y, ...
    activePrev2y, ...
    raw.nActiveBy20y, ...
    activePrev20y, ...
    postActive20, ...
    postActivePrev20y, ...
    raw.nPreventedActiveTB, ...
    relReduction20, ...
    raw.NNS_preventActiveTB, ...
    raw.NNS_cureInfection ...
    };

Median = zeros(numel(Metric),1);
Low95  = zeros(numel(Metric),1);
High95 = zeros(numel(Metric),1);
for i = 1:numel(X)
    [Median(i), Low95(i), High95(i)] = summarise_numeric(X{i});
end

summary = table(Metric, Median, Low95, High95);

% Add user-friendly percentage formatting for prevalence/reduction rows.
IsPercent = contains(Metric, "prevalence") | contains(Metric, "reduction");
summary.DisplayScale = repmat("count", height(summary), 1);
summary.DisplayScale(IsPercent) = "proportion";

out = struct();
out.summary = summary;
out.derived = table(raw.nInfected, raw.nActiveBy2y, raw.nActiveBy20y, postActive20, ...
    raw.nPreventedActiveTB, relReduction20, ltbiPrev, activePrev2y, activePrev20y, postActivePrev20y, ...
    'VariableNames', {'nInfected','nActiveBy2y_DoNothing','nActiveBy20y_DoNothing', ...
    'nActiveBy20y_AfterStrategy','nActiveBy20y_Prevented','relReduction20y', ...
    'ltbiPrev_DoNothing','activePrev2y_DoNothing','activePrev20y_DoNothing','activePrev20y_AfterStrategy'});
out.resultsUsed = results;

try
    writetable(summary, fullfile(outDir, 'tb_do_nothing_summary_v9.csv'));
catch
end
end

function validate_v9_results(results)
if ~isstruct(results) || ~isfield(results, 'raw') || ~istable(results.raw)
    error('Input must be a v9 results struct containing a raw results table.');
end
req = {'nInfected','nActiveBy2y','nActiveBy20y','nPreventedActiveTB','NNS_preventActiveTB','NNS_cureInfection'};
for i = 1:numel(req)
    if ~ismember(req{i}, results.raw.Properties.VariableNames)
        error('Results table is missing required column: %s', req{i});
    end
end
if ~isfield(results, 'settings') || ~isfield(results.settings, 'N')
    error('Results struct is missing settings.N');
end
end

function y = safe_fraction_vec(num, den)
if isscalar(den)
    den = repmat(den, size(num));
end
y = num ./ den;
y(den == 0) = NaN;
end

function [med, low, high] = summarise_numeric(x)
x = x(:);
x = x(~isnan(x));
if isempty(x)
    med = NaN; low = NaN; high = NaN; return;
end
x = sort(x);
med = median(x);
low = quantile_sorted(x, 0.025);
high = quantile_sorted(x, 0.975);
end

function q = quantile_sorted(x, p)
% x must be sorted column vector.
n = numel(x);
if n == 1
    q = x; return;
end
pos = 1 + (n - 1) * p;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    q = x(lo);
else
    q = x(lo) + (pos - lo) * (x(hi) - x(lo));
end
end
