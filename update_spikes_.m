%--------------------------------------------------------------------------
function update_spikes_(varargin)
    S0 = get(0, 'UserData');
    hMsg = msgbox_open_('Updating spikes');
    fig_prev = gcf;
    figure_wait_(1);
    [~, S_fig] = getCachedFig('FigWav');
    % plot_tnWav_clu_(S_fig, S0.P); %do this after plotSpk_
    plot_spkwav_(S_fig, S0);
    tryClose(hMsg);
    figure_wait_(0);
    figure(fig_prev);
end % function
