%--------------------------------------------------------------------------
function raw_waveform_(hMenu)
    figure_wait_(1); drawnow;
    [P, S_clu] = get0_('P', 'S_clu');
    if get_(P, 'fWav_raw_show')
        P.fWav_raw_show = 0;
    else
        P.fWav_raw_show = 1;
    end

    if isempty(get_(S_clu, 'tmrWav_raw_clu'))
        S_clu = S_clu_wav_(S_clu);
        S0 = set0_(P, S_clu);
    else
        S0 = set0_(P);
    end
    set(hMenu, 'Checked', ifeq_(P.fWav_raw_show, 'on', 'off'));
    % redraw windows
    plot_FigWav_(S0);
    button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0);
    figure_wait_(0);
end %func
