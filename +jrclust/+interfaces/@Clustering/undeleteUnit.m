function res = undeleteUnit(obj, spikeClusters, unitID, indices, metadata)
    %UNDELETEUNITS Speculatively undelete a unit, returning a snapshot of
    %the spike table and metadata fields, along with a diff.
    if nargin < 5
        metadata = struct();
    end
    res = struct('spikeClusters', [], ...
                 'diffTable', [], ...
                 'metadata', []);
    
    if isempty(indices)
        return;
    end

    goodUnits = unique(spikeClusters(spikeClusters > 0));
    nClusters = numel(goodUnits);

    % ensure all indices are in bounds
    if any((indices < 1) | (indices > numel(spikeClusters)))
        error('indices exceed bounds');
    end

    % check for spikes which haven't actually been deleted
    if any(spikeClusters(indices) > 0)
        error('non-deleted spikes in deleted unit');
    end

    % augment fields
    isConsistent = 1;

    vectorFields = obj.unitFields.vectorFields;
    hFunConsistent = @(vals) numel(vals) == nClusters + 1; % ensure we have the expected number of clusters after an augment

    for i = 1:numel(vectorFields)
        fn = vectorFields{i};
        if ~isfield(metadata, fn)
            metadata.(fn) = obj.(fn);
        end

        if ~isempty(metadata.(fn))
            if isrow(metadata.(fn))
                hFunAug = @(vals, augmentAfter) [vals(1:augmentAfter) empty(vals, 1) vals(augmentAfter+1:end)];
            else
                hFunAug = @(vals, augmentAfter) [vals(1:augmentAfter); empty(vals, 1); vals(augmentAfter+1:end)];
            end
            metadata.(fn) = hFunAug(metadata.(fn), unitID-1);

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
                        augShape = [nClusters + 1, 1];
                        
                    case 'clusterCentroids'
                        augShape = [1, 2];

                    otherwise
                        clusterDims = otherFields.(fn).cluster_dims;
                        % get remaining (non-cluster-indexed) dimensions
                        shape = size(metadata.(fn));
                        augShape = [shape(1:clusterDims-1) 1 shape(clusterDims+1:end)];
                end
                filler = empty(metadata.(fn), augShape);

                hFunAug = eval(otherFields.(fn).augment);
                hFunConsistent = eval(otherFields.(fn).consistent); % see /json/Clustering.json for subsetting functions for non-vector fields
                metadata.(fn) = hFunAug(metadata.(fn), unitID-1, filler);

                if ~hFunConsistent(metadata.(fn), nClusters + 1)
                    isConsistent = 0;
                    break;
                end
            end
        end
    end

    if isConsistent
        % first row: indices of spikes to undelete
        % second row: old unit ID of spikes to undelete
        % third row: new unit IDs of spikes to undelete
        spikeClusters = spikeClusters(:);
        diffTable = [indices(:)'; spikeClusters(indices)'; zeros(1, numel(indices))+unitID];

        % side effect: shift all larger units up by unity
        mask = (spikeClusters >= unitID);
        spikeClusters(mask) = spikeClusters(mask) + 1;

        % AFTER making room for new units, assign split off spikes to their
        % new units
        spikeClusters(diffTable(1, :)) = unitID;

        % update spikesByCluster since we already have this information
        if isfield(metadata, 'spikesByCluster')
            metadata.spikesByCluster{unitID} = double(diffTable(1, :)');
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
