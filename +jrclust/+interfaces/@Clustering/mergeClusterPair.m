function success = mergeClusterPair(obj, iCluster, jCluster)
    %MERGECLUSTERPAIR Merge a pair of clusters
    success = 0;
    if iCluster == jCluster || numel(intersect([iCluster, jCluster], obj.spikeClusters)) < 2
        return;
    end

    % keep the smaller of the two and shift others left by one
    [iCluster, jCluster] = deal(min(iCluster, jCluster), max(iCluster, jCluster));

    [clustersOld, clustersNew] = deal(obj.spikeClusters);
    
    clustersNew(obj.spikeClusters == jCluster) = iCluster;
    % shift all clusters larger than j down by 1
    clustersNew(clustersNew > jCluster) = clustersNew(clustersNew > jCluster) - 1;
    obj.spikeClusters = clustersNew;
    obj.subsetFields([1:(jCluster-1) (jCluster+1):(obj.nClusters+1)]);

    if ~isempty(obj.inconsistentFields())
        warning('Cluster data is inconsistent after merging %d and %d', iCluster, jCluster);
        obj.spikeClusters = clustersOld;
        success = 0;
        obj.refresh(1, []);
    else
        obj.postOp(iCluster);
        obj.orderClusters('clusterSites');
        success = 1;
    end
end