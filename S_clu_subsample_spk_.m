%--------------------------------------------------------------------------
function [viSpk_clu1, viSite_clu1, vlSpk_clu1] = S_clu_subsample_spk_(S_clu, iClu, S0)
    % subsample spikes from the requested cluster centered at the center site and mid-time range (drift)

    fSelect_mid = 0; % TW
    nSamples_max = 1000;
    if nargin<3, S0 = get(0, 'UserData'); end

    % [P, viSite_spk] = get0_('P', 'viSite_spk'); end
    [viSpk_clu1, viSite_clu1, vlSpk_clu1] = deal([]);
    % Subselect based on the center site
    viSpk_clu1 = S_clu.cviSpk_clu{iClu}; %
    if isempty(viSpk_clu1), return; end
    iSite_clu1 = S_clu.viSite_clu(iClu);
    vlSpk_clu1 = iSite_clu1 == S0.viSite_spk(viSpk_clu1);
    viSite_clu1 = S0.P.miSites(:,iSite_clu1);
    viSpk_clu1 = viSpk_clu1(vlSpk_clu1);
    if isempty(viSpk_clu1), return; end

    if fSelect_mid
        viSpk_clu1 = spk_select_mid_(viSpk_clu1, S0.viTime_spk, S0.P);
    end
    viSpk_clu1 = subsample_vr_(viSpk_clu1, nSamples_max);
end %func
