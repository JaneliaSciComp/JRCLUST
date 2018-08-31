%--------------------------------------------------------------------------
function updateFigTime(hFig, figData)

    S0 = get(0, 'UserData');
    if ~isVisible_(figData.hAx)
        return;
    end

    [vrFet0, vrTime0, vcYlabel] = getFigTimeFeatures(figData.primarySite, [], S0);
    if ~isfield(figData, 'fPlot0')
        figData.fPlot0 = 1;
    end
    toggleVisible_(figData.hPlotBG, figData.fPlot0);
    updatePlot(figData.hPlotBG, vrTime0, vrFet0);

    set(figData.hPlotFG, 'YData', getFigTimeFeatures(figData.primarySite, S0.primarySelectedCluster, S0));

    if ~isempty(S0.secondarySelectedCluster)
        set(figData.hPlotFG2, 'YData', getFigTimeFeatures(figData.primarySite, S0.secondarySelectedCluster, S0));
    else
        clearPlots(figData.hPlotFG2);
    end

    autoscale_pct = getOr(S0.P, 'autoscale_pct', 99.5);
    limits = [0 1.5] * quantile([figData.hPlotFG.YData figData.hPlotFG2.YData], autoscale_pct/100);
    figData.maxAmp = limits(2);
    set(hFig, 'UserData', figData);
    ylim_(figData.hAx, limits);
    imrect_set_(figData.hRect, [], limits);

    grid(figData.hAx, 'on');
    ylabel(figData.hAx, vcYlabel);
end % function
