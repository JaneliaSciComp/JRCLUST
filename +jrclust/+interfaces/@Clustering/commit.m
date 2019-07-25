function success = commit(obj, spikeClusters, metadata, msg)
    %COMMIT Commit a modification of clustering to history log
    success = 0;
    if obj.nEdits > obj.editPos % not at tip of edit history, back out
        warning('cannot branch from history; use revert() first');
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

    % flag units which need to have fields recomputed
    if ~isfield(metadata, 'unitCount')
        obj.spikesByCluster = cell(obj.nClusters, 1);
        obj.clusterNotes = cell(obj.nClusters, 1);
        obj.unitCount = zeros(obj.nClusters, 1);
        flagged = 1:obj.nClusters;
    else
        flagged = find(isnan(metadata.unitCount));
    end
    
    % apply new spike table and metadata entries, update history
    try
        obj.spikeClusters = spikeClusters;

        for i = 1:numel(fieldNames)
            fn = fieldNames{i};
            if isprop(obj, fn)
                obj.(fn) = metadata.(fn);
            end
        end

        nEntries = int32(size(obj.history, 1) + 1);
        obj.history(nEntries) = msg;

        % open history file and write new spike table to it
        if exist(obj.hCfg.histFile, 'file') ~= 2
            error('history file %s does not exist', obj.hCfg.histFile);
        end

        fidHist = fopen(obj.hCfg.histFile, 'a');
        fwrite(fidHist, int32(nEntries), 'int32');
        fwrite(fidHist, int32(spikeClusters), 'int32');
        fclose(fidHist);
    catch ME
        % restore spike table and metadata entries from backup
        fieldNames = fieldnames(backup);
        for i = 1:numel(fieldNames)
            fn = fieldNames{i};
            obj.(fn) = backup.(fn);
        end

        warning('failed to commit: %s', ME.message);
        return;
    end

    % update edit position
    if isempty(obj.editPos)
        obj.editPos = 1;
    else
        obj.editPos = obj.editPos + 1;
    end

    if ~isempty(flagged)
        obj.spikesByCluster(flagged) = arrayfun(@(iC) find(obj.spikeClusters == iC), flagged, 'UniformOutput', 0);
        % update count of spikes per unit
        obj.unitCount(flagged) = cellfun(@numel, obj.spikesByCluster(flagged));

        % update cluster sites
        obj.clusterSites(flagged) = double(arrayfun(@(iC) mode(obj.spikeSites(obj.spikesByCluster{iC})), flagged));

        % update mean waveforms for flagged units
        obj.updateWaveforms(flagged);

        % update unit positions
        obj.computeCentroids(flagged);

        % compute quality scores for altered units
        obj.computeQualityScores(flagged);
    end

    success = 1;
end