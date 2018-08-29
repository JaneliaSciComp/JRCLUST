%--------------------------------------------------------------------------
function [clusterSpikes, iCenteredSpikes] = getClusterCenteredSpikes(S_clu, clusterID, spikeSites)
    % get a subset of cluster that is centered
    % return only centered spikes
    % if nargin<2, S0 = get(0, 'UserData'); end
    % S_clu = S0.S_clu;
    if nargin < 3
        spikeSites = get0_('spikeSites');
    end

    centerSite = S_clu.clusterSites(clusterID);
    clusterSpikes = S_clu.spikesByCluster{clusterID};
    clusterSites = spikeSites(clusterSpikes);
    iCenteredSpikes = find(clusterSites == centerSite);
    clusterSpikes = clusterSpikes(iCenteredSpikes);
end % function
