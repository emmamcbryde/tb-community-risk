function out = run_tb_screening_igra_addons_v9(regimen, nReps, coverageGrid, Nprofile)
% Convenience driver for the v9 add-ons.
% Runs:
%   1) risk-factor composition profiling for cure/prevent targeting
%   2) IGRA example charts comparing targeted strategies with random

if nargin < 1 || isempty(regimen)
    regimen = '3HP';
end
if nargin < 2 || isempty(nReps)
    nReps = 2000;
end
if nargin < 3 || isempty(coverageGrid)
    coverageGrid = [0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.80 1.00];
end
if nargin < 4 || isempty(Nprofile)
    Nprofile = 50000;
end

profile = run_tb_screening_targeting_profile_v9(regimen, coverageGrid, Nprofile, 1);
charts = run_tb_screening_igra_charts_v9(regimen, nReps, coverageGrid);

out = struct();
out.profile = profile;
out.charts = charts;
end
