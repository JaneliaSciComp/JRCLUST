%--------------------------------------------------------------------------
function export_tmrWav_clu_(hObject, event)
    % exports spike waveforms by clusters
    [S_clu, P] = get0_('S_clu', 'P');
    tmrWav_clu = S_clu.tmrWav_clu;
    assignWorkspace_(tmrWav_clu);
end % function
