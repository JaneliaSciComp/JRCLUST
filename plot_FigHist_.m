%--------------------------------------------------------------------------
function plot_FigHist_(S0)

    if nargin<1, S0 = get(0, 'UserData'); end
    S_clu = S0.S_clu; P = S0.P;
    [hFig, S_fig] = get_fig_cache_('FigHist');

    nBins_hist = 50; % @TODO: put this in param file

    vrX = logspace(0, 4, nBins_hist);
    vrY1 = isi_hist_(S0.iCluCopy, vrX);
    vcTitle = sprintf('Cluster %d', S0.iCluCopy);

    % draw
    if isempty(S_fig) %first time the iCluPaste is always empty
        S_fig.hAx = axes_new_(hFig);
        S_fig.hPlot1 = stairs(S_fig.hAx, nan, nan, 'k');
        S_fig.hPlot2 = stairs(S_fig.hAx, nan, nan, 'r');
        xlim_(S_fig.hAx, [1 10000]); %in msec
        grid(S_fig.hAx, 'on');
        xlabel(S_fig.hAx, 'ISI (ms)');
        ylabel(S_fig.hAx, 'Prob. Density');
        set(S_fig.hAx, 'XScale', 'log');
    end
    update_plot_(S_fig.hPlot1, vrX, vrY1);
    if ~isempty(S0.iCluPaste)
        vrY2 = isi_hist_(S0.iCluPaste, vrX);
        vcTitle = sprintf('Cluster %d (black) vs %d (red)', S0.iCluCopy, S0.iCluPaste);
        update_plot_(S_fig.hPlot2, vrX, vrY2);
    else
        update_plot_(S_fig.hPlot2, nan, nan);
    end
    title_(S_fig.hAx, vcTitle);

    set(hFig, 'UserData', S_fig);
end %func
