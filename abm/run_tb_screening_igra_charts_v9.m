function out = run_tb_screening_igra_charts_v9(regimen, nReps, coverageGrid)
% Create IGRA example charts comparing targeted strategies with random.
% Default regimen is now 3HP.
%
% Figures written:
%   1. proportion treated vs mean cured infection (random vs cure-targeted)
%   2. proportion treated vs mean prevented active TB (random vs prevent-targeted)
%   3. proportion treated vs pooled NNS to cure infection (random vs cure-targeted)
%   4. proportion treated vs pooled NNS to prevent active TB (random vs prevent-targeted)
%
% Outputs:
%   out.plotData          table used for the plots
%   out.dataCsv           CSV path for the plot data
%   out.figureFiles       struct of PNG paths

thisFile = mfilename('fullpath');
[thisDir, ~, ~] = fileparts(thisFile);
outDir = get_output_dir_v9();
csvFile = fullfile(thisDir, 'default_data.csv');

if nargin < 1 || isempty(regimen)
    regimen = '3HP';
end
if nargin < 2 || isempty(nReps)
    nReps = 2000;
end
if nargin < 3 || isempty(coverageGrid)
    coverageGrid = [0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.80 1.00];
end

N = 1500;
strategies = {'random','cure','prevent'};
rows = repmat(struct( ...
    'screeningStrategy', '', ...
    'screenCoverage', NaN, ...
    'meanScreened', NaN, ...
    'meanStartedTPT', NaN, ...
    'propTreated', NaN, ...
    'meanCuredInfection', NaN, ...
    'meanPreventedActiveTB', NaN, ...
    'pooledNNS_cure', NaN, ...
    'pooledNNS_prevent', NaN), numel(strategies)*numel(coverageGrid), 1);

row = 0;
for s = 1:numel(strategies)
    for c = 1:numel(coverageGrid)
        row = row + 1;
        res = tb_screening_mc_model_v9(csvFile, ...
            'N', N, ...
            'nReps', nReps, ...
            'screenWindow', 2, ...
            'followHorizon', 20, ...
            'testType', 'IGRA', ...
            'regimen', regimen, ...
            'screeningStrategy', strategies{s}, ...
            'screenCoverage', coverageGrid(c), ...
            'seed', 1);

        rows(row).screeningStrategy = strategies{s};
        rows(row).screenCoverage = coverageGrid(c);
        rows(row).meanScreened = mean(res.raw.nScreened);
        rows(row).meanStartedTPT = mean(res.raw.nStartTPT);
        rows(row).propTreated = mean(res.raw.nStartTPT) / N;
        rows(row).meanCuredInfection = mean(res.raw.nCuredInfection);
        rows(row).meanPreventedActiveTB = mean(res.raw.nPreventedActiveTB);
        rows(row).pooledNNS_cure = sum(res.raw.nScreened) / sum(res.raw.nCuredInfection);
        rows(row).pooledNNS_prevent = sum(res.raw.nScreened) / sum(res.raw.nPreventedActiveTB);
    end
end

plotData = struct2table(rows, 'AsArray', true);
plotData = sortrows(plotData, {'screeningStrategy','screenCoverage'});

dataCsv = fullfile(outDir, sprintf('igra_%s_strategy_plot_data_v9.csv', lower(regimen)));
writetable(plotData, dataCsv);

% Subset for the direct comparisons requested.
plotCure = plotData(ismember(plotData.screeningStrategy, {'random','cure'}), :);
plotPrevent = plotData(ismember(plotData.screeningStrategy, {'random','prevent'}), :);

fig1 = figure('Name', 'IGRA cure yield', 'Color', 'w');
hold on
make_strategy_plot(plotCure, 'meanCuredInfection');
xlabel('Proportion treated')
ylabel('Mean infections cured / protected per cohort')
title(sprintf('IGRA + %s: treated proportion vs infections cured', regimen))
grid on
legend('Location','best')
ylim([0 inf])
png1 = fullfile(outDir, sprintf('igra_%s_prop_treated_vs_cured_v9.png', lower(regimen)));
saveas(fig1, png1);

fig2 = figure('Name', 'IGRA prevention yield', 'Color', 'w');
hold on
make_strategy_plot(plotPrevent, 'meanPreventedActiveTB');
xlabel('Proportion treated')
ylabel('Mean active TB cases prevented per cohort')
title(sprintf('IGRA + %s: treated proportion vs active TB prevented', regimen))
grid on
legend('Location','best')
ylim([0 inf])
png2 = fullfile(outDir, sprintf('igra_%s_prop_treated_vs_prevented_v9.png', lower(regimen)));
saveas(fig2, png2);

fig3 = figure('Name', 'IGRA cure NNS', 'Color', 'w');
hold on
make_strategy_plot(plotCure, 'pooledNNS_cure');
xlabel('Proportion treated')
ylabel('Pooled NNS to cure / protect one infection')
title(sprintf('IGRA + %s: treated proportion vs NNS to cure infection', regimen))
grid on
legend('Location','best')
ylim([0 inf])
png3 = fullfile(outDir, sprintf('igra_%s_prop_treated_vs_nns_cure_v9.png', lower(regimen)));
saveas(fig3, png3);

fig4 = figure('Name', 'IGRA prevent NNS', 'Color', 'w');
hold on
make_strategy_plot(plotPrevent, 'pooledNNS_prevent');
xlabel('Proportion treated')
ylabel('Pooled NNS to prevent one active TB case')
title(sprintf('IGRA + %s: treated proportion vs NNS to prevent active TB', regimen))
grid on
legend('Location','best')
ylim([0 inf])
png4 = fullfile(outDir, sprintf('igra_%s_prop_treated_vs_nns_prevent_v9.png', lower(regimen)));
saveas(fig4, png4);

out = struct();
out.plotData = plotData;
out.dataCsv = dataCsv;
out.figureFiles = struct( ...
    'propTreatedVsCured', png1, ...
    'propTreatedVsPrevented', png2, ...
    'propTreatedVsNNSCure', png3, ...
    'propTreatedVsNNSPrevent', png4);

fprintf('Wrote IGRA plot data to: %s\n', dataCsv);
fprintf('Saved figure: %s\n', png1);
fprintf('Saved figure: %s\n', png2);
fprintf('Saved figure: %s\n', png3);
fprintf('Saved figure: %s\n', png4);
end

function make_strategy_plot(tbl, yField)
strategies = unique(tbl.screeningStrategy, 'stable');
for i = 1:numel(strategies)
    idx = strcmp(tbl.screeningStrategy, strategies{i});
    x = tbl.propTreated(idx);
    y = tbl.(yField)(idx);
    y(~isfinite(y)) = NaN;
    plot(x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', strategies{i});
end
end
