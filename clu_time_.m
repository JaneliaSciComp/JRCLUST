%--------------------------------------------------------------------------
function [viTime_clu1, viSpk_clu1] = clu_time_(iClu1)
    % returns time in sec
    [S_clu, viTime_spk] = get0_('S_clu', 'viTime_spk');
    viSpk_clu1 = S_clu.cviSpk_clu{iClu1};
    viTime_clu1 = viTime_spk(S_clu.cviSpk_clu{iClu1});
end %func
