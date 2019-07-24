function reassign(obj, recompute)
    %REASSIGN Reassign clusters, e.g., after a change of parameters
    if obj.nEdits > obj.editPos % not at tip of edit history, back out
        warning('cannot branch from history; use revert() first');
        return;
    end

    obj.sRes = jrclust.sort.assignClusters(obj.dRes, obj.sRes, obj.hCfg);

    % these fields are mutable so we need to store copies in obj
    obj.spikeClusters = obj.initialClustering;
    if isfield(obj.sRes, 'clusterCenters')
        obj.clusterCenters = obj.sRes.clusterCenters;
    else
        obj.clusterCenters = [];
    end
    obj.clusterCentroids = [];

    if recompute
        obj.updateWaveforms();
        obj.computeCentroids();
        obj.computeQualityScores([]);
    end

    obj.orderClusters('clusterSites');
end