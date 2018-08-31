%--------------------------------------------------------------------------
function [viTime_clu1, viSpk_clu1] = clu_time_(iClu1)
    % returns time in sec
    [S_clu, spikeTimes] = get0_('S_clu', 'spikeTimes');
    viSpk_clu1 = S_clu.spikesByCluster{iClu1};
    viTime_clu1 = spikeTimes(S_clu.spikesByCluster{iClu1});
end % function
