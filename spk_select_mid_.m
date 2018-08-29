%--------------------------------------------------------------------------
function viSpk_clu2 = spk_select_mid_(viSpk_clu1, spikeTimes, P)
    % spikeTimes = get0_('spikeTimes');
    iSpk_mid = round(numel(spikeTimes)/2);
    viSpk_clu1_ord = rankorder_(abs(viSpk_clu1 - iSpk_mid), 'ascend');
    nSpk1_max = round(numel(viSpk_clu1) / P.nTime_clu);
    viSpk_clu2 = viSpk_clu1(viSpk_clu1_ord <= nSpk1_max);
end % function
