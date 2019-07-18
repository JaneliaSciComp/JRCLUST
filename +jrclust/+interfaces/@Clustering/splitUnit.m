function res = splitUnit(obj, spikeClusters, unitID, unitPart, metadata)
    %SPLITUNIT Speculatively split a unit, returning a snapshot of the
    %spike table and metadata fields, along with a diff.
    if nargin < 5
        metadata = struct();
    end
    res = struct('spikeClusters', [], ...
                 'diffTable', [], ...
                 'metadata', []);

    nSplits = numel(unitPart); % number of *new* units created after a split
    unitIndices = find(spikeClusters == unitID);
    goodUnits = unique(spikeClusters(spikeClusters > 0));
    nClusters = numel(goodUnits);
    
    if isempty(unitIndices)
        return;
    end

    unitPart = cellfun(@(x) x(:)', unitPart, 'UniformOutput', 0);

    % check for duplicates
    for i = 1:numel(unitPart)
        iPart = unitPart{i};
        if numel(unique(iPart)) < numel(iPart)
            error('duplicate indices in partition %d', i);
        end
        for j = i+1:numel(unitPart)
            jPart = unitPart{j};
            if ~isempty(intersect(iPart, jPart))
                error('duplicate indices between partition %d and partition %d', i, j);
            end
        end
    end

    % replace local (in-unit) indices with global (spike-table) indices
    splitOff = [unitPart{:}]';
    if max(splitOff) > numel(unitIndices) || min(splitOff) < 1
        error('split indices exceed bounds');
    end
    unitPart = cellfun(@(i) unitIndices(i)', unitPart, 'UniformOutput', 0);

    % augment fields
    isConsistent = 1;

    vectorFields = obj.unitFields.vectorFields;
    hFunConsistent = @(vals) numel(vals) == nClusters + nSplits; % ensure we have the expected number of clusters after a deletion

    for i = 1:numel(vectorFields)
        fn = vectorFields{i};
        if ~isfield(metadata, fn)
            metadata.(fn) = obj.(fn);
        end

        if ~isempty(metadata.(fn))
            if isrow(metadata.(fn))
                hFunAug = @(vals, augmentAfter) [vals(1:augmentAfter) empty(vals, [1 nSplits]) vals(augmentAfter+1:end)];
            else
                hFunAug = @(vals, augmentAfter) [vals(1:augmentAfter); empty(vals, [nSplits 1]); vals(augmentAfter+1:end)];
            end
            metadata.(fn) = hFunAug(metadata.(fn), unitID);

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
                switch fn
                    case {'templateSim', 'waveformSim'}
                        augShape = [nClusters + nSplits, nSplits];
                        
                    case 'clusterCentroids'
                        augShape = [nSplits, 2];

                    otherwise
                        clusterDims = otherFields.(fn).cluster_dims;
                        % get remaining (non-cluster-indexed) dimensions
                        shape = size(metadata.(fn));
                        augShape = [shape(1:clusterDims-1) nSplits shape(clusterDims+1:end)];
                end
                filler = empty(metadata.(fn), augShape);

                hFunAug = eval(otherFields.(fn).augment);
                hFunConsistent = eval(otherFields.(fn).consistent); % see /json/Clustering.json for subsetting functions for non-vector fields
                metadata.(fn) = hFunAug(metadata.(fn), unitID, filler);

                if ~hFunConsistent(metadata.(fn), nClusters + nSplits)
                    isConsistent = 0;
                    break;
                end
            end
        end
    end

    if isConsistent
        % first row: indices of spikes to split off
        % second row: old unit ID of spikes to split off
        % third row: new unit IDs of spikes to split off
        indices = [unitPart{:}];
        newUnits = arrayfun(@(i) zeros(1, numel(unitPart{i}))+unitID+i, 1:nSplits, ...
                            'UniformOutput', 0);
        newUnits = [newUnits{:}];
        diffTable = [indices; zeros(1, numel(indices))+unitID; newUnits];
        
        % side effect: shift all larger units up by unity
        mask = (spikeClusters > unitID);
        spikeClusters(mask) = spikeClusters(mask) + nSplits;

        % AFTER making room for new units, assign split off spikes to their
        % new units
        spikeClusters(diffTable(1, :)) = diffTable(3, :);

        % update spikesByCluster since we already have this information
        if isfield(metadata, 'spikesByCluster')
            metadata.spikesByCluster{unitID} = setdiff(metadata.spikesByCluster{unitID}, diffTable(1, :));
            for unit = unique(diffTable(3, :))
                metadata.spikesByCluster{unit} = diffTable(1, diffTable(3, :) == unit)';
            end
        end

        res.spikeClusters = spikeClusters;
        res.diffTable = diffTable;
        res.metadata = metadata;
    end
end

function filler = empty(vals, shape)
    cls = class(vals);
    switch cls
        case 'cell' % expecting a homogeneous cell array, e.g., of char
            filler = cell(shape);
            if ~isempty(vals)
                filler = cellfun(@(x) cast(x, 'like', vals{1}), filler, 'UniformOutput', 0);
            end

        case 'char'
            filler = char(shape);

        otherwise
            filler = zeros(shape, cls);
    end
end
