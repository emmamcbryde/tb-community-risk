function fieldName = regimen_cost_field_v9(regimen)
%REGIMEN_COST_FIELD_V9 Map regimen labels to MATLAB-safe cost field names.

regimenText = char(string(regimen));
switch upper(strtrim(regimenText))
    case '3HP'
        fieldName = 'x3HP';
    case '4R'
        fieldName = 'x4R';
    case '3HR'
        fieldName = 'x3HR';
    case '6H'
        fieldName = 'x6H';
    case '9H'
        fieldName = 'x9H';
    otherwise
        error('regimen_cost_field_v9:InvalidRegimen', ...
            'Unsupported regimen for economics costs: %s', regimenText);
end
end
