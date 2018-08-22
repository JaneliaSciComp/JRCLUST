%--------------------------------------------------------------------------
function rescale_FigWav_(event, S0, P)
    % figure_wait_(1);
    set(0, 'UserData', S0);

    [S_fig, maxAmp_prev, hFigWav] = set_fig_maxAmp_('FigWav', event);
    set_fig_(hFigWav, plot_tnWav_clu_(S_fig, P));
    multiplot(S0.hCopy, S_fig.maxAmp);
    if ~isempty(S0.secondarySelectedCluster)
        multiplot(S0.hPaste, S_fig.maxAmp);
    end
    rescale_spikes_(S_fig.hSpkAll, maxAmp_prev, P);
    title_(S_fig.hAx, sprintf(S_fig.vcTitle, S_fig.maxAmp)); %update scale
end
