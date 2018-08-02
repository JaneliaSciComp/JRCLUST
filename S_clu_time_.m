%--------------------------------------------------------------------------
function [viTime1, viSpk1, viSpk2] = S_clu_time_(S_clu, iClu)
    % return time of cluster time in adc sample index unit.
    % viSpk2: spike indices of centered spikes in cluster iClu

    S0 = get(0, 'UserData');
    if isfield(S_clu, 'cviSpk_clu')
        viSpk1 = S_clu.cviSpk_clu{iClu};
        viTime1 = S0.spikeTimes(viSpk1);
    else
        viSpk1 = find(S_clu.spikeClusters == iClu);
        viTime1 = S0.spikeTimes(viSpk1);
    end
    if nargout>=3
        iSite1 = S_clu.clusterSites(iClu);
        viSpk2 = viSpk1(S0.viSite_spk(viSpk1) == iSite1);
    end
end %func
