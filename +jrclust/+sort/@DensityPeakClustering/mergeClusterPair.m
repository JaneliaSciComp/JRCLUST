function success = mergeClusterPair(obj, iCluster, jCluster)
%MERGECLUSTERPAIR Merge a pair of clusters
    % take the denser of the two centers, using delta as tiebreaker
    iCenter = obj.clusterCenters(iCluster);
    jCenter = obj.clusterCenters(jCluster);
    if obj.spikeRho(iCenter) < obj.spikeRho(jCenter)
        obj.clusterCenters(iCluster) = jCenter;
    elseif obj.spikeRho(iCenter) == obj.spikeRho(jCenter) && obj.spikeDelta(iCenter) < obj.spikeDelta(jCenter)
        obj.clusterCenters(iCluster) = jCenter;
    end % otherwise, keep iCenter

    success = mergeClusterPair@jrclust.interfaces.Clustering(obj, iCluster, jCluster);
end

