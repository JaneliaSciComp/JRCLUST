function success = deleteClusters(obj, deleteMe)
    %DELETECLUSTERS Delete clusters
    success = 1;

    clustersBak = obj.spikeClusters; % in case we need to restore

    nClustersOld = numel(obj.spikesByCluster);
    keepMe = setdiff(1:nClustersOld, deleteMe);
    obj.subsetFields(keepMe);

    garbageCluster = min(obj.spikeClusters) - 1;
    if garbageCluster == 0 % reserved for noise
        garbageCluster = -1;
    end

    deleteSpikes = ismember(obj.spikeClusters, deleteMe);
    obj.spikeClusters(deleteSpikes) = garbageCluster;
    nClusters_ = numel(keepMe);

    good = (obj.spikeClusters > 0);
    mapFrom = zeros(1, nClustersOld);
    mapFrom(keepMe) = 1:nClusters_;

    obj.spikeClusters(good) = mapFrom(obj.spikeClusters(good));

    if ~isempty(obj.inconsistentFields())
        warning('Cluster data is inconsistent after deleting %d', deleteMe);
        obj.spikeClusters = clustersBak;
        success = 0;
        obj.refresh(1, []);
    end
end