function exported = export_results_v9(results, baseName, outDir)
%EXPORT_RESULTS_V9 Export v9 results tables to the standard output folder.

if nargin < 2 || isempty(baseName)
    baseName = 'tb_screening_v9';
end
if nargin < 3 || isempty(outDir)
    outDir = get_output_dir_v9();
end

summary = summarise_results_v9(results);

summaryCsv = fullfile(outDir, sprintf('%s_summary.csv', baseName));
keyCsv = fullfile(outDir, sprintf('%s_key_metrics.csv', baseName));

writetable(summary.summaryTable, summaryCsv);
writetable(summary.keyMetrics, keyCsv);

exported = struct();
exported.summaryCsv = summaryCsv;
exported.keyMetricsCsv = keyCsv;
end
