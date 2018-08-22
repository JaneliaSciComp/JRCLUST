%--------------------------------------------------------------------------
function [S_clu, vlKeep_clu] = S_clu_refresh_(S_clu, fRemoveEmpty)

    if nargin<2, fRemoveEmpty=1; end
    nClu = double(max(S_clu.spikeClusters));
    S_clu.nClusters = nClu;
    spikeSites = get0_('spikeSites');

    % agglomerate spikes by cluster
    S_clu.spikesByCluster = arrayfun(@(iClu)find(S_clu.spikeClusters==iClu), 1:nClu, 'UniformOutput', 0);
    % count them
    S_clu.nSpikesPerCluster = cellfun(@numel, S_clu.spikesByCluster);
    % assign most-frequently (per cluster) occurring site to cluster
    S_clu.clusterSites = double(arrayfun(@(iClu)mode(spikeSites(S_clu.spikesByCluster{iClu})), 1:nClu));
    if fRemoveEmpty, [S_clu, vlKeep_clu] = S_clu_remove_empty_(S_clu, 1); end
end %func
