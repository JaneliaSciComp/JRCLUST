function res = mergeUnits(obj, spikeClusters, targetUnit, mergingUnits, metadata)
    %MERGEUNITS Speculatively merge a group of units, returning a snapshot of the
    %spike table and metadata fields, along with a diff.
    if nargin < 5
        metadata = struct();
    end
    res = struct('spikeClusters', [], ...
                 'diffTable', [], ...
                 'metadata', []);

    % don't try to merge into ourselves
    mergingUnits = setdiff(mergingUnits, targetUnit);

    nMerging = numel(mergingUnits);
    indices = find(ismember(spikeClusters, mergingUnits));
    goodUnits = unique(spikeClusters(spikeClusters > 0));
    nClusters = numel(goodUnits);
    
    if isempty(indices)
        return;
    end

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
        % first row: indices of spikes to merge
        % second row: old unit IDs of spikes to merge
        % third row: new unit IDs of spikes to merge
        diffTable = [indices(:)'; spikeClusters(indices)'; zeros(1, numel(indices))+targetUnit];
        spikeClusters(indices) = targetUnit;

        % side effect: shift all larger units down by unity
        mergingUnits = sort(mergingUnits, 'descend');
        for i = mergingUnits
            mask = (spikeClusters > i);
            spikeClusters(mask) = spikeClusters(mask) - 1;
        end

        % update spikesByCluster since we already have this information
        if isfield(metadata, 'spikesByCluster')
            metadata.spikesByCluster{targetUnit} = sort([metadata.spikesByCluster{targetUnit}; diffTable(1, :)']);
        end

        res.spikeClusters = spikeClusters;
        res.diffTable = diffTable;
        res.metadata = metadata;
    end
end
