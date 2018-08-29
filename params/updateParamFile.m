%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function updateParamFile(P, filename)
    % Modify the parameter file using the variables in the P struct

    csLines = file2cellstr_(filename); %read to cell string
    csLines_var = first_string_(csLines);

    paramKeys = fieldnames(P);
    paramValues = cellfun(@(field) P.(field), paramKeys, 'UniformOutput', 0);
    for i = 1:numel(paramKeys)
        vcName = paramKeys{i}; % find field name with
        if isstruct(paramValues{i})
            continue;
        end % do not write struct

        vcValue = field2str_(paramValues{i});
        iLine = find(strcmpi(csLines_var, vcName));

        if numel(iLine)>1 % more than one variable found
            error(['updateParamFile: Multiple copies of variables found: ' vcName]);
        elseif isempty(iLine) %did not find, append
            csLines{end+1} = sprintf('%s = %s;', vcName, vcValue);
        else
            vcComment = getCommentExpr_(csLines{iLine});
            csLines{iLine} = sprintf('%s = %s;\t\t\t%s', vcName, vcValue, vcComment);
        end
    end
    cellstr2file_(filename, csLines);
end % function
