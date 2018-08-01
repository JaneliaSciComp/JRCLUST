%--------------------------------------------------------------------------
% 12/28/17 JJJ: originalStruct can be empty (v3.2.1)
% 8/4/17 JJJ: selective struct merge
% 7/31/17 JJJ: documentation and testing
function finalStruct = mergeStructs(originalStruct, structToMerge, fieldsToMerge)
    % Merge second struct to first one
    % originalStruct = mergeStructs(originalStruct, P_append)
    % originalStruct = mergeStructs(originalStruct, P_append, var_list) : only update list of variable names
    if isempty(originalStruct)
        finalStruct = structToMerge;
    elseif isempty(structToMerge)
        finalStruct = originalStruct;
    else
        finalStruct = originalStruct;
        if nargin < 3 % no fields specified, merge all fields
            fieldsToMerge = fieldnames(structToMerge);
        end

        if ischar(fieldsToMerge) % single string
            fieldsToMerge = { fieldsToMerge };
        end

        for iField = 1:numel(fieldsToMerge)
            fieldName = fieldsToMerge{iField};
            if isfield(structToMerge, fieldName)
                finalStruct.(fieldName) = structToMerge.(fieldName);
            end
        end
    end
end %func
