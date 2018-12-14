function [clusterSpikes, neighbors] = subsampleCenteredSpikes(hClust, iCluster)
    %SUBSAMPLECENTEREDSPIKES Subsample spikes from the requested cluster centered at the center site and mid-time range (drift)
    nSamplesMax = 1000;

    % subselect based on the center site
    clusterSpikes = hClust.spikesByCluster{iCluster};
    if isempty(clusterSpikes)
        neighbors = [];
        return;
    end

    iClusterSite = hClust.clusterSites(iCluster);
    centeredSpikes = (hClust.spikeSites(clusterSpikes) == iClusterSite);
    neighbors = hClust.hCfg.siteNeighbors(:, iClusterSite);
    clusterSpikes = clusterSpikes(centeredSpikes);
    if isempty(clusterSpikes)
        return;
    end

    clusterSpikes = jrclust.utils.subsample(clusterSpikes, nSamplesMax);
end