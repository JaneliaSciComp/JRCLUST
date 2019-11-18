function success = revert(obj, revertTo)
    %REVERT Delete history
    success = 0;
    if revertTo < 0 || revertTo > obj.nEdits
        return;
    end

    % load history from file
    if exist(obj.hCfg.histFile, 'file') ~= 2
        error('history file %s not found!', obj.hCfg.histFile);
    end

    fidHist = fopen(obj.hCfg.histFile, 'r');
    checkInt = fread(fidHist, 1, 'int32');
    fRes = -isempty(checkInt); % -1 if no checkInt read, 0 if okay
    while fRes > -1 && ~isempty(checkInt) && checkInt ~= revertTo
        fRes = fseek(fidHist, 4*obj.nSpikes, 'cof');
        checkInt = fread(fidHist, 1, 'int32');
    end

    if isempty(checkInt) || checkInt ~= revertTo
        fclose(fidHist);
        warning('failed to revert: entry %d not found in history file', revertTo);
    end

    spikeClusters = fread(fidHist, obj.nSpikes, 'int32');
    fclose(fidHist);

    spikesByCluster = arrayfun(@(iC) find(spikeClusters == iC), 1:max(spikeClusters), 'UniformOutput', 0);
    flagged = 1:numel(spikesByCluster);
    metadata = struct();

    isConsistent = 1;

    vectorFields = obj.unitFields.vectorFields;
    hFunConsistent = @(vals) numel(vals) == numel(spikesByCluster); % ensure we have the expected number of clusters after a deletion

    for i = 1:numel(vectorFields)
        fn = vectorFields{i};
        metadata.(fn) = obj.(fn);

        if ~isempty(metadata.(fn))
            if numel(spikesByCluster) < obj.nClusters % losing clusters
                hFun = @(vals, indices) vals(indices);
                metadata.(fn) = hFun(metadata.(fn), flagged);
            elseif numel(spikesByCluster) > obj.nClusters && isrow(metadata.(fn)) % gaining clusters
                hFun = @(vals, nAugs) [vals empty(vals, [1 nAugs])];
                metadata.(fn) = hFun(metadata.(fn), numel(spikesByCluster) - obj.nClusters);
            elseif numel(spikesByCluster) > obj.nClusters
                hFun = @(vals, nAugs) [vals; empty(vals, [nAugs 1])];
                metadata.(fn) = hFun(metadata.(fn), numel(spikesByCluster) - obj.nClusters);
            end

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
            metadata.(fn) = obj.(fn);

            if ~isempty(metadata.(fn))
                if numel(spikesByCluster) < obj.nClusters % losing clusters
                    hFunSubset = eval(otherFields.(fn).subset);
                    metadata.(fn) = hFunSubset(metadata.(fn), flagged);
                elseif numel(spikesByCluster) > obj.nClusters
                    nExtra = numel(spikesByCluster) - obj.nClusters;
                    switch fn
                        case {'templateSim', 'waveformSim'}
                            augShape = [numel(spikesByCluster), nExtra];

                        case 'clusterCentroids'
                            augShape = [nExtra, 2];

                        otherwise
                            clusterDims = otherFields.(fn).cluster_dims;
                            % get remaining (non-cluster-indexed) dimensions
                            shape = size(metadata.(fn));
                            augShape = [shape(1:clusterDims-1) nExtra shape(clusterDims+1:end)];
                    end
                    filler = empty(metadata.(fn), augShape);

                    hFunAug = eval(otherFields.(fn).augment);
                    metadata.(fn) = hFunAug(metadata.(fn), 0, filler);
                end
                hFunConsistent = eval(otherFields.(fn).consistent); % see /json/Clustering.json for subsetting functions for non-vector fields

                if ~hFunConsistent(metadata.(fn), numel(spikesByCluster))
                    isConsistent = 0;
                    break;
                end
            end
        end
    end

    if ~isConsistent
        warning('failed to revert: inconsistent fields');
        return;
    end

    % create a backup of changed fields
    backup = struct('spikeClusters', obj.spikeClusters);
    fieldNames = fieldnames(metadata);
    for i = 1:numel(fieldNames)
        fn = fieldNames{i};
        if isprop(obj, fn)
            backup.(fn) = obj.(fn);
        end
    end

    try
        obj.spikeClusters = spikeClusters;

        for i = 1:numel(fieldNames)
            fn = fieldNames{i};
            if isprop(obj, fn)
                obj.(fn) = metadata.(fn);
            end
        end

        histkeys = keys(obj.history);
        remove(obj.history, num2cell(setdiff([histkeys{:}], 1:revertTo))); %% FIXME

        obj.syncHistFile();
    catch ME
        % restore spike table and metadata entries from backup
        fieldNames = fieldnames(backup);
        for i = 1:numel(fieldNames)
            fn = fieldNames{i};
            obj.(fn) = backup.(fn);
        end

        warning('failed to revert: %s', ME.message);
        return;
    end

    obj.editPos = revertTo;

    if ~isempty(flagged)
        obj.spikesByCluster(flagged) = arrayfun(@(iC) find(obj.spikeClusters == iC), flagged, 'UniformOutput', 0);
        % update count of spikes per unit
        obj.unitCount(flagged) = cellfun(@numel, obj.spikesByCluster(flagged));

        % update cluster sites
%         if ~isempty(obj.spikeSites)
            obj.clusterSites(flagged) = double(arrayfun(@(iC) mode(obj.spikeSites(obj.spikesByCluster{iC})), flagged));
%         end

        % don't recompute mean waveforms (and their consequences) if we didn't have them already
%         if ~isempty(obj.meanWfGlobal)
            % update mean waveforms for flagged units
            obj.updateWaveforms(flagged);

            % update unit positions
            obj.computeCentroids(flagged);

            % compute quality scores for altered units
            obj.computeQualityScores(flagged);
%         end
    end

    success = 1;
end

%% LOCAL FUNCTIONS
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
