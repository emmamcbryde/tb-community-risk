function displayModel = results_bundle_to_display_v9(bundle)
%RESULTS_BUNDLE_TO_DISPLAY_V9 Convert results bundle to UI-friendly display model.

displayModel = struct();
displayModel.metadata = get_or_empty(bundle, 'metadata');
displayModel.validation = get_or_empty(bundle, 'validation');
displayModel.headline = get_or_empty(bundle, 'headline');
displayModel.technical = get_or_empty(bundle, 'technical');
displayModel.targetedVsRandom = get_or_empty(bundle, 'targetedVsRandom');
displayModel.doNothing = get_or_empty(bundle, 'doNothing');
displayModel.attributableRisk = get_or_empty(bundle, 'attributableRisk');
displayModel.charts = get_or_empty(bundle, 'charts');
displayModel.economics = economics_display_or_empty(bundle);
displayModel.downloads = get_or_empty(bundle, 'downloads');

displayModel.visibleSections = {};
sections = fieldnames(displayModel);
for i = 1:numel(sections)
    name = sections{i};
    section = displayModel.(name);
    if isstruct(section) && isfield(section, 'available') && logical(section.available)
        displayModel.visibleSections{end+1} = name; %#ok<AGROW>
    end
end
end

function value = get_or_empty(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = struct();
end
end

function value = economics_display_or_empty(bundle)
value = get_or_empty(bundle, 'economics');
if ~isstruct(value) || ~isfield(value, 'available') || ~logical(value.available)
    return;
end
if isfield(value, 'summaryTable') && istable(value.summaryTable)
    value.summaryRows = table2struct(value.summaryTable, 'ToScalar', false);
end
end
