function S1 = mergeStructs(S1, S2, fieldNames)
    %MERGESTRUCTS Merge a pair of structs
    if isempty(S1)
        S1 = S2;
        return;
    elseif isempty(S2)
        return;
    end

    if nargin < 3
        fieldNames = fieldnames(S2);
    else % fieldNames can be a cell or a char here
        fieldNames = intersect(fieldNames, fieldnames(S2));
    end

    for iField = 1:numel(fieldNames)
        S1.(fieldNames{iField}) = S2.(fieldNames{iField});
    end
end
