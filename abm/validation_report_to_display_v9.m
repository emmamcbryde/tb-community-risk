function displayModel = validation_report_to_display_v9(report)
%VALIDATION_REPORT_TO_DISPLAY_V9 Convert validation report to UI-friendly display model.

displayModel = struct();
displayModel.isValid = false;
displayModel.hasWarnings = false;
displayModel.errorCount = 0;
displayModel.warningCount = 0;
displayModel.infoCount = 0;
displayModel.errorRows = struct.empty(0, 1);
displayModel.warningRows = struct.empty(0, 1);
displayModel.infoRows = struct.empty(0, 1);
displayModel.summaryText = '';

if isempty(report) || ~isstruct(report)
    displayModel.summaryText = 'No validation report available.';
    return;
end

displayModel.isValid = logical(report.isValid);
displayModel.hasWarnings = logical(report.hasWarnings);
displayModel.errorRows = issue_array_or_empty(report, 'errors');
displayModel.warningRows = issue_array_or_empty(report, 'warnings');
displayModel.infoRows = issue_array_or_empty(report, 'infos');
displayModel.errorCount = numel(displayModel.errorRows);
displayModel.warningCount = numel(displayModel.warningRows);
displayModel.infoCount = numel(displayModel.infoRows);
displayModel.summaryText = sprintf('Valid: %d | Errors: %d | Warnings: %d | Infos: %d', ...
    displayModel.isValid, displayModel.errorCount, displayModel.warningCount, displayModel.infoCount);
end

function rows = issue_array_or_empty(report, fieldName)
if isfield(report, fieldName) && isstruct(report.(fieldName))
    rows = report.(fieldName);
else
    rows = struct.empty(0, 1);
end
end
