%--------------------------------------------------------------------------
function plot_FigHist_(S0)

    if nargin<1, S0 = get(0, 'UserData'); end
    S_clu = S0.S_clu; P = S0.P;
    [hFig, S_fig] = getCachedFig('FigHist');

    nBins_hist = 50; % @TODO: put this in param file

    vrX = logspace(0, 4, nBins_hist);
    vrY1 = isi_hist_(S0.primarySelectedCluster, vrX);
    vcTitle = sprintf('Cluster %d', S0.primarySelectedCluster);

    % draw
    if isempty(S_fig) %first time the secondarySelectedCluster is always empty
        S_fig.hAx = newAxes(hFig);
        S_fig.hPlotFG = stairs(S_fig.hAx, nan, nan, 'k');
        S_fig.hPlotFG2 = stairs(S_fig.hAx, nan, nan, 'r');
        xlim_(S_fig.hAx, [1 10000]); %in msec
        grid(S_fig.hAx, 'on');
        xlabel(S_fig.hAx, 'ISI (ms)');
        ylabel(S_fig.hAx, 'Prob. Density');
        set(S_fig.hAx, 'XScale', 'log');
    end
    updatePlot(S_fig.hPlotFG, vrX, vrY1);
    if ~isempty(S0.secondarySelectedCluster)
        vrY2 = isi_hist_(S0.secondarySelectedCluster, vrX);
        vcTitle = sprintf('Cluster %d (black) vs %d (red)', S0.primarySelectedCluster, S0.secondarySelectedCluster);
        updatePlot(S_fig.hPlotFG2, vrX, vrY2);
    else
        updatePlot(S_fig.hPlotFG2, nan, nan);
    end
    title_(S_fig.hAx, vcTitle);

    set(hFig, 'UserData', S_fig);
end % function
