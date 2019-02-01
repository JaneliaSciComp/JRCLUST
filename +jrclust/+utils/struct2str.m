function strstruct = struct2str(S)
    %STRUCT2STR Convert a struct to a string
    strstruct = '';
    if ~isstruct(S)
        return;
    end

    fieldNames = fieldnames(S);
    for iField = 1:numel(fieldNames)
        strstruct = sprintf('%s%s = %s;', strstruct, fieldNames{iField}, jrclust.utils.field2str(S.(fieldNames{iField})));
        if iField < numel(fieldNames)
            strstruct = sprintf('%s\n', strstruct);
        end
    end %for
end
