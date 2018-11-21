%--------------------------------------------------------------------------
function viSpk_clu2 = spk_select_mid_(viSpk_clu1, viTime_spk, P)
    % viTime_spk = get0_('viTime_spk');
    iSpk_mid = round(numel(viTime_spk)/2);
    viSpk_clu1_ord = jrclust.utils.rankorder(abs(viSpk_clu1 - iSpk_mid), 'ascend');
    nSpk1_max = round(numel(viSpk_clu1) / P.nTime_clu);
    viSpk_clu2 = viSpk_clu1(viSpk_clu1_ord <= nSpk1_max);
end %func
