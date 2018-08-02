%--------------------------------------------------------------------------
function [viSpk_clu1, viiSpk_clu1] = S_clu_viSpk_(S_clu, iClu1, spikeSites)
    % get a subset of cluster that is centered
    % return only centered spikes
    % if nargin<2, S0 = get(0, 'UserData'); end
    % S_clu = S0.S_clu;
    if nargin<3, spikeSites = get0_('spikeSites'); end
    iSite_clu1 = S_clu.clusterSites(iClu1);
    viSpk_clu1 = S_clu.cviSpk_clu{iClu1};
    clusterSites1 = spikeSites(viSpk_clu1);
    viiSpk_clu1 = find(clusterSites1 == iSite_clu1);
    viSpk_clu1 = viSpk_clu1(viiSpk_clu1);
end %func
