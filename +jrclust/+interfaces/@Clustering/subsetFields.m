function subsetFields(obj, keepMe)
    %SUBSETFIELDS Subset all data fields, taking only those indices we want
    % to keep, prior to rearranging or deleting
    vectorFields = obj.unitFields.vectorFields;
    hFun = @(vals, indices) vals(indices);
    for i = 1:numel(vectorFields)
        fn = vectorFields{i};
        if ~isempty(obj.(fn))
            obj.(fn) = hFun(obj.(fn), keepMe);
        end
    end

    % see JRCLUST/json/Clustering.json for subsetting functions for
    % non-vector fields
    otherFields = obj.unitFields.otherFields;
    otherFieldNames = fieldnames(otherFields);
    for i = 1:numel(otherFieldNames)
        fn = otherFieldNames{i};
        if ~isempty(obj.(fn))
            hFun = eval(otherFields.(fn).subset);
            obj.(fn) = hFun(obj.(fn), keepMe);
        end
    end
end
