%--------------------------------------------------------------------------
function [viSpk_clu1, viiSpk_clu1] = S_clu_viSpk_(S_clu, iClu1, viSite_spk)
    % get a subset of cluster that is centered
    % return only centered spikes
    % if nargin<2, S0 = get(0, 'UserData'); end
    % S_clu = S0.S_clu;
    if nargin<3, viSite_spk = get0_('viSite_spk'); end
    iSite_clu1 = S_clu.viSite_clu(iClu1);
    viSpk_clu1 = S_clu.cviSpk_clu{iClu1};
    viSite_clu1 = viSite_spk(viSpk_clu1);
    viiSpk_clu1 = find(viSite_clu1 == iSite_clu1);
    viSpk_clu1 = viSpk_clu1(viiSpk_clu1);
end %func
