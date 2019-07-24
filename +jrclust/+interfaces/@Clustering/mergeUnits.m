function res = mergeUnits(obj, spikeClusters, mergeTargets, mergingUnits, metadata)
    %MERGEUNITS Speculatively merge a pair of units, returning a snapshot of
    %the spike table and metadata fields.
    if nargin < 5
        metadata = struct();
    end
    res = struct('spikeClusters', [], ...
                 'metadata', []);

    % don't try to merge into ourselves
    intersection = ismember(mergingUnits, mergeTargets);
    if any(intersection)
        mergingUnits(intersection) = [];
        mergeTargets(intersection) = [];
        if isempty(mergingUnits)
            warning('merge targets found among merging units; none left to merge!')
            return;
        else
            warning('merge targets found among merging units; removing duplicates');
        end
    end

    nMerging = numel(mergingUnits);
    indices = find(ismember(spikeClusters, mergingUnits));
    goodUnits = unique(spikeClusters(spikeClusters > 0));
    nClusters = numel(goodUnits);
    
    if isempty(indices)
        return;
    end

    % flag these units to update later
    metadata.unitCount = obj.unitCount;
    metadata.unitCount(unique(mergeTargets)) = nan;

    % subset fields
    keepSubset = setdiff(goodUnits, mergingUnits); % units to keep, i.e., not the ones that will be merged
    isConsistent = 1;

    vectorFields = obj.unitFields.vectorFields;
    hFunSubset = @(vals, indices) vals(indices);
    hFunConsistent = @(vals) numel(vals) == nClusters - nMerging; % ensure we have the expected number of clusters after a deletion

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

                if ~hFunConsistent(metadata.(fn), nClusters - nMerging)
                    isConsistent = 0;
                    break;
                end
            end
        end
    end

    if isConsistent
        for i = 1:numel(mergingUnits)
            spikeClusters(spikeClusters == mergingUnits(i)) = mergeTargets(i);
        end

        % side effect: shift all larger units down by unity
        mergingUnits = sort(mergingUnits, 'descend');
        for i = mergingUnits
            mask = (spikeClusters > i);
            spikeClusters(mask) = spikeClusters(mask) - 1;
        end

        res.spikeClusters = spikeClusters;
        res.metadata = metadata;
    end
end
