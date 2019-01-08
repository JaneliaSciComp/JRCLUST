%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function edit_prm_file_(P, vcFile_prm)
    % Modify the parameter file using the variables in the P struct

    csLines = file2cellstr_(vcFile_prm); %read to cell string
    csLines_var = first_string_(csLines);

    csName = fieldnames(P);
    csValue = cellfun(@(vcField)P.(vcField), csName, 'UniformOutput',0);
    for i=1:numel(csName)
        vcName = csName{i}; %find field name with
        if isstruct(csValue{i}), continue; end %do not write struct
        vcValue = jrclust.utils.field2str(csValue{i});
        iLine = find(strcmpi(csLines_var, vcName));
        if numel(iLine)>1 % more than one variable found
            error(['edit_prm_file_: Multiple copies of variables found: ' vcName]);
        elseif isempty(iLine) %did not find, append
            csLines{end+1} = sprintf('%s = %s;', vcName, vcValue);
        else
            vcComment = getCommentExpr_(csLines{iLine});
            csLines{iLine} = sprintf('%s = %s;\t\t\t%s', vcName, vcValue, vcComment);
        end
    end
    cellstr2file_(vcFile_prm, csLines);
end %func
