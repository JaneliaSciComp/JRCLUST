function inconsistencies = inconsistentFields(obj)
    %INCONSISTENTFIELDS Check all fields have the correct sizes
    inconsistencies = {};

    vectorFields = obj.unitFields.vectorFields;
    hFun = @(vals) numel(vals) == obj.nClusters;
    for i = 1:numel(vectorFields)
        fn = vectorFields{i};
        if ~isempty(obj.(fn)) && ~hFun(obj.(fn))
            inconsistencies{end+1} = vectorFields{i}; %#ok<AGROW>
        end
    end

    % see JRCLUST/json/Clustering.json for consistency functions for
    % non-vector fields
    otherFields = obj.unitFields.otherFields;
    otherFieldNames = fieldnames(otherFields);
    for i = 1:numel(otherFieldNames)
        fn = otherFieldNames{i};
        hFun = eval(otherFields.(fn).consistent);
        if ~isempty(obj.(fn)) && ~hFun(obj.(fn), obj.nClusters)
            inconsistencies{end+1} = vectorFields{i}; %#ok<AGROW>
        end
    end
end