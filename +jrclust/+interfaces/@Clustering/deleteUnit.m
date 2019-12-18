function res = deleteUnit(obj, spikeClusters, unitID, metadata)
    %DELETEUNIT Speculatively delete a unit, returning a snapshot of the
    %spike table and metadata fields.
    if nargin < 4
        metadata = struct();
    end
    res = struct('spikeClusters', [], ...
                 'metadata', []);

    indices = find(spikeClusters == unitID);
    goodUnits = unique(spikeClusters(spikeClusters > 0));
    nClusters = numel(goodUnits);
    
    if isempty(indices)
        return;
    end

    % subset fields
    keepSubset = setdiff(goodUnits, unitID); % units to keep, i.e., not this one
    isConsistent = 1;

    vectorFields = obj.unitFields.vectorFields;
    hFunSubset = @(vals, indices) vals(indices);
    hFunConsistent = @(vals) numel(vals) == nClusters - 1; % ensure we have the expected number of clusters after a deletion

    for i = 1:numel(vectorFields)
        fn = vectorFields{i};
        if ~isfield(metadata, fn)
            metadata.(fn) = obj.(fn);
        end

        if ~isempty(metadata.(fn))
            metadata.(fn) = hFunSubset(metadata.(fn), keepSubset);

            if ~hFunConsistent(metadata.(fn))
                isConsistent = 0;
                break;
            end
        end
    end

    if isConsistent
        otherFields = obj.unitFields.otherFields;
        otherFieldNames = fieldnames(otherFields);

        for i = 1:numel(otherFieldNames)
            fn = otherFieldNames{i};
            if ~isfield(metadata, fn)
                metadata.(fn) = obj.(fn);
            end

            if ~isempty(metadata.(fn))
                hFunSubset = eval(otherFields.(fn).subset);
                hFunConsistent = eval(otherFields.(fn).consistent); % see /json/Clustering.json for subsetting functions for non-vector fields
                metadata.(fn) = hFunSubset(metadata.(fn), keepSubset);

                if ~hFunConsistent(metadata.(fn), nClusters - 1)
                    isConsistent = 0;
                    break;
                end
            end
        end
    end

    if isConsistent
        spikeClusters(indices) = -1;
        
        % side effect: shift all larger units down by unity
        mask = (spikeClusters >= unitID);
        spikeClusters(mask) = spikeClusters(mask) - 1;

        res.spikeClusters = spikeClusters;
        res.metadata = metadata;
    end
end
