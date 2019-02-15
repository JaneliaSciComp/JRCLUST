function [centeredSpikes, whichCentered] = getCenteredSpikes(obj, iCluster)
    %GETCENTEREDSPIKES Return subset of spikes which occur on center site of cluster
    iClusterSite = obj.clusterSites(iCluster);
    centeredSpikes = obj.spikesByCluster{iCluster};
    clusterSites_ = obj.spikeSites(centeredSpikes);
    whichCentered = clusterSites_ == iClusterSite;
    centeredSpikes = centeredSpikes(whichCentered);
end