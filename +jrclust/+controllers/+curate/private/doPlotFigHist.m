function hFigHist = doPlotFigHist(hFigHist, hClust, hCfg, selected)
    %DOPLOTFIGHIST Plot ISI histogram
    if numel(selected) == 1
        iCluster = selected(1);
        jCluster = iCluster;
    else
        iCluster = selected(1);
        jCluster = selected(2);
    end

    nBinsHist = 50; %TODO: put this in param file (maybe)

    XData = logspace(0, 4, nBinsHist);
    YData1 = getISIHistogram(iCluster, XData, hClust, hCfg);

    % draw the plot
    if isempty(hFigHist.figData)
        hFigHist.axes();
        hFigHist.addPlot('hPlot1', @stairs, nan, nan, 'Color', hCfg.mrColor_proj(2, :));
        hFigHist.addPlot('hPlot2', @stairs, nan, nan, 'Color', hCfg.mrColor_proj(3, :));
        hFigHist.axApply(@set, 'XLim', [1 10000], 'XScale', 'log'); % ms
        hFigHist.axApply(@grid, 'on');
        hFigHist.axApply(@xlabel, 'ISI (ms)');
        hFigHist.axApply(@ylabel, 'Prob. Density');
    end

    hFigHist.updatePlot('hPlot1', XData, YData1);

    if iCluster ~= jCluster
        YData2 = getISIHistogram(jCluster, XData, hClust, hCfg);
        hFigHist.axApply(@title, sprintf('Cluster %d (black) vs. %d (red)', iCluster, jCluster), 'Interpreter', 'none', 'FontWeight', 'normal');
        hFigHist.updatePlot('hPlot2', XData, YData2);
    else
        hFigHist.axApply(@title, sprintf('Cluster %d', iCluster), 'Interpreter', 'none', 'FontWeight', 'normal');
        hFigHist.hidePlot('hPlot2');
    end
end

%% LOCAL FUNCTIONS
function YData = getISIHistogram(iCluster, XData, hClust, hCfg)
    %GETISIHISTOGRAM Get a histogram of ISI values
    %   TODO: replace hist call with histogram
    clusterTimes = double(hClust.spikeTimes(hClust.spikesByCluster{iCluster}))/hCfg.sampleRate;
    YData = hist(diff(clusterTimes)*1000, XData);
    YData(end) = 0;
    YData = YData./sum(YData); % normalize
end
