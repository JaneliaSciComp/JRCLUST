%function [S_clu, vlKeep_clu] = S_clu_refresh_(S_clu, fRemoveEmpty)
function [clusterData, vlKeep_clu] = S_clu_refresh_(clusterData, spikeSites, fRemoveEmpty)
    %
    if nargin < 3
        fRemoveEmpty = true;
    end

    nClusters = double(max(clusterData.spikeClusters));
    clusterData.nClusters = nClusters;

    %cviSpk_clu
    clusterData.spikesByCluster = arrayfun(@(iC) find(clusterData.spikeClusters == iC), 1:nClusters, 'UniformOutput', 0);
    %vnSpk_clu
    clusterData.clusterCounts = cellfun(@numel, clusterData.spikesByCluster);
    %viSite_clu
    clusterData.clusterSites = double(arrayfun(@(iC) mode(spikeSites(clusterData.spikesByCluster{iC})), 1:nClusters));

    if fRemoveEmpty
        [clusterData, vlKeep_clu] = S_clu_remove_empty_(clusterData);
    end
end
