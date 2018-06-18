%--------------------------------------------------------------------------
function S_clu = post_recenter_(S_clu, P) % recenters the features and
    % trFet_spk, viSite_spk, viSite2_spk, cviSpk_site, cviSpk2_site
    global trFet_spk %nFet/site x nFet x nSpk
    [viSite_spk, viSite2_spk] = get0_('viSite_spk', 'viSite2_spk');
    nSites = numel(P.viSite2Chan);
    % vpCentered_clu = arrayfun(@(iClu)mean(viSite_spk(S_clu.cviSpk_clu{iClu}) == S_clu.viSite_clu(iClu) | viSite2_spk(S_clu.cviSpk_clu{iClu}) == S_clu.viSite_clu(iClu)), 1:S_clu.nClu);
    % find spikes not centered
    for iClu = 1:S_clu.nClu
        viSpk_clu1 = S_clu.cviSpk_clu{iClu};
        vlSpk_swap_clu1 = viSite2_spk(viSpk_clu1) == S_clu.viSite_clu(iClu);
        vi_swap_ = viSpk_clu1(vlSpk_swap_clu1);
        % swap indices
        %     [viSite1_spk, viSite2_spk] = swap_vr_(viSite_spk, viSite2_spk, viSpk_swap_clu1);
        [viSite_spk(vi_swap_), viSite2_spk(vi_swap_)] = deal(viSite2_spk(vi_swap_), viSite_spk(vi_swap_));
        trFet_spk_clu1 = trFet_spk(:,:,vi_swap_);
        trFet_spk(:,1,vi_swap_) = trFet_spk_clu1(:,2,:);
        trFet_spk(:,2,vi_swap_) = trFet_spk_clu1(:,1,:);
    end %for
    cviSpk_site = arrayfun(@(iSite)find(viSite_spk == iSite), 1:nSites, 'UniformOutput', 0);
    cviSpk2_site = arrayfun(@(iSite)find(viSite2_spk == iSite), 1:nSites, 'UniformOutput', 0);
    S0 = set0_(viSite_spk, viSite2_spk, cviSpk_site, cviSpk2_site);
end %func
