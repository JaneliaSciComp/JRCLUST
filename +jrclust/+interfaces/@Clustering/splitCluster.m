function [success, retained] = splitCluster(obj, iCluster, retained)
    %SPLITCLUSTER Split a cluster
    success = 0;
    if iCluster < 1 || iCluster > obj.nClusters
        return;
    end

    iSpikes = find(obj.spikeClusters == iCluster);
    if ~all(ismember(retained, iSpikes))
        return;
    end

    clustersBak = obj.spikeClusters;
    splitOff = iSpikes(~ismember(iSpikes, retained));

    % swap retained and splitOff if iSite > jSite
    iSite = mode(obj.spikeSites(retained));
    jSite = mode(obj.spikeSites(splitOff));
    if iSite > jSite
        [retained, splitOff] = deal(splitOff, retained);
    end

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
    mask = (obj.spikeClusters > iCluster);
    obj.spikeClusters(mask) = obj.spikeClusters(mask) + 1;
    obj.spikeClusters(splitOff) = iCluster + 1;

    % swap fields
    obj.subsetFields([1:iCluster obj.nClusters iCluster+1:obj.nClusters-1]);

    if isempty(obj.inconsistentFields())
        success = 1;
        obj.postOp([iCluster, iCluster + 1]);
        %obj.orderClusters('clusterSites'); % can be done manually if the user really wants
    else
        warning('Cluster data is inconsistent after splitting %d', iCluster);
        obj.spikeClusters = clustersBak;
        success = 0;
        obj.subsetFields([1:iCluster iCluster+1:obj.nClusters+1]);
        obj.refresh(1, []);
    end
end