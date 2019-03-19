function success = partitionCluster(obj, iCluster, assignPart)
    %PARTITIONCLUSTER Break up a cluster
    success = 0;
    if iCluster < 1 || iCluster > obj.nClusters
        return;
    end

    initialSpikes = find(obj.spikeClusters == iCluster);
    clustersBak = obj.spikeClusters;

    for iSplit = 1:numel(assignPart)-1
        if iSplit > 1
            iSpikes = find(obj.spikeClusters == iCluster);
        else
            iSpikes = initialSpikes;
        end

        splitOff = initialSpikes(assignPart{iSplit});
        retain = iSpikes(~ismember(iSpikes, splitOff));

        iSite = mode(obj.spikeSites(retain));
        jSite = mode(obj.spikeSites(splitOff));

        % swap retain and split-off spikes; always make splitOff the next
        % cluster
        if iSite > jSite
            splitOff = retain;
        end

        %%% augment fields by one
        fieldNames = fieldnames(obj);
        % augment vector fields
        vecFields = fieldNames(cellfun(@(fn) isvector(obj.(fn)) && numel(obj.(fn)) == obj.nClusters, fieldNames));
        for iField = 1:numel(vecFields)
            fn = vecFields{iField};
            fd = obj.(fn);

            if iscell(fd)
                fd{end+1} = []; %#ok<*AGROW>
            else
                fd(end+1) = 0;
            end
            obj.(fn) = fd;
        end

        % augment matrix fields
        if ~isempty(obj.clusterCentroids) % nClusters x 2
            obj.clusterCentroids(end+1, :) = zeros(1, 2, 'like', obj.clusterCentroids);
        end
        if ~isempty(obj.waveformSim) % nClusters x nClusters
            obj.waveformSim = [obj.waveformSim zeros(obj.nClusters, 1, 'like', obj.waveformSim); ...
                            zeros(1, obj.nClusters + 1, 'like', obj.waveformSim)];
        end

        % augment tensor fields
        tenFields = fieldNames(cellfun(@(fn) ndims(obj.(fn)) == 3 && size(obj.(fn), 3) == obj.nClusters, fieldNames));
        for iField = 1:numel(tenFields)
            fn = tenFields{iField};
            fd = obj.(fn);
            fd(:, :, end+1) = zeros(size(fd, 1), size(fd, 2), 'like', fd);
            obj.(fn) = fd;
        end

        % make room for new cluster
        mask = (obj.spikeClusters > iCluster + iSplit - 1);
        obj.spikeClusters(mask) = obj.spikeClusters(mask) + 1;
        obj.spikeClusters(splitOff) = iCluster + iSplit;

        % swap fields
        obj.subsetFields([1:(iCluster+iSplit-1) obj.nClusters (iCluster+iSplit):obj.nClusters-1]);
    end

    if isempty(obj.inconsistentFields())
        success = 1;
        obj.postOp(iCluster:(iCluster+numel(assignPart)-1));
        %obj.orderClusters('clusterSites'); % can be done manually if the user really wants
    else
        warning('Cluster data is inconsistent after splitting %d', iCluster);
        obj.spikeClusters = clustersBak;
        success = 0;
        obj.subsetFields([1:iCluster (iCluster+numel(assignPart)):obj.nClusters]);
        obj.refresh(1, []);
    end
end