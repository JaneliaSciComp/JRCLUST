function hFigPreview = doPlotFigPreview(hFigPreview, figData, fKeepView, hCfg)
    %DOPLOTFIGPREVIEW
    hWait = jrclust.utils.qMsgBox('Plotting...', 0, 1);

    nSites = size(figData.tracesFilt, 2);

    XDataSamp = figData.windowBounds(1):figData.windowBounds(2);
    XData = XDataSamp / hCfg.sampleRate;
    tlimSecs = (figData.windowBounds + [-1 1]) / hCfg.sampleRate;

    %% Mean plot
    if ~hFigPreview.hasPlot('hPlotMean')
        hFigPreview.addPlot('hPlotMean', @plot, hFigPreview.hAxes('hAxMean'), nan, nan, 'k');
        hFigPreview.addPlot('hPlotMeanThresh', @plot, hFigPreview.hAxes('hAxMean'), nan, nan, 'r');
    end

    if strcmp(figData.refView, 'original')
        YData = figData.tracesCAR(XDataSamp);
    else % binned
        YData = figData.channelMeansMAD(XDataSamp);
    end

    hFigPreview.plotApply('hPlotMean', @set, 'XData', XData, 'YData', YData);
    hFigPreview.axApply('hAxMean', @set, 'YLim', [0, 100]);
    hFigPreview.axApply('hAxMean', @xlabel, 'Time (sec)');
    hFigPreview.axApply('hAxMean', @ylabel, sprintf('Common ref. (MAD, %s)', figData.refView));

    fThreshRef = strcmpi(figData.refView, 'binned') && ~isempty(figData.blankThresh) && figData.blankThresh ~= 0;
    if fThreshRef
        hFigPreview.plotApply('hPlotMeanThresh', @set, 'XData', XData([1, end]), 'YData', repmat(figData.blankThresh, [1, 2]));
    else
        hFigPreview.hidePlot('hPlotMeanThresh');
    end
    hFigPreview.axApply('hAxMean', @grid, jrclust.utils.ifEq(figData.fGrid, 'on', 'off'));

    %% Traces plot
    if ~hFigPreview.hasPlot('hPlotTraces')
        hFigPreview.addPlot('hPlotTraces', @plot, hFigPreview.hAxes('hAxTraces'), nan, nan, 'Color', [1 1 1]*.5);
        hFigPreview.addPlot('hPlot_traces_spk', @plot, hFigPreview.hAxes('hAxTraces'), nan, nan, 'm.-', 'LineWidth', 1.5);
        hFigPreview.addPlot('hPlot_traces_spk1', @plot, hFigPreview.hAxes('hAxTraces'), nan, nan, 'ro');
        hFigPreview.addPlot('hPlotTracesThresh', @plot, hFigPreview.hAxes('hAxTraces'), nan, nan, 'm:');
        hFigPreview.addPlot('hPlotTracesBad', @plot, hFigPreview.hAxes('hAxTraces'), nan, nan, 'r');
    end

    if figData.fFilter
        hCfg.setTemporaryParams('filterType', figData.filterType);
        traces = jrclust.utils.bit2uV(figData.tracesFilt(XDataSamp, :), hCfg);
        hCfg.resetTemporaryParams();
    else
        traces = jrclust.utils.meanSubtract(single(figData.tracesClean(XDataSamp, :))*hCfg.bitScaling);
    end

    hFigPreview.multiplot('hPlotTraces', figData.maxAmp, XData, traces, 1:nSites);

    % plot bad sites in red
    if ~isempty(figData.ignoreSites)
        hFigPreview.multiplot('hPlotTracesBad', figData.maxAmp, XData, traces(:, figData.ignoreSites), figData.ignoreSites);
    else
        hFigPreview.hidePlot('hPlotTracesBad');
    end

    if figData.fShowSpikes
        inBounds = figData.spikeTimes >= figData.windowBounds(1) & figData.spikeTimes <= figData.windowBounds(end);
        spikeTimes = single(figData.spikeTimes(inBounds) - figData.windowBounds(1)+1);
        spikeTimesSec = single(figData.spikeTimes(inBounds)) / hCfg.sampleRate;
        spikeSites = single(figData.spikeSites(inBounds));
    else
        spikeTimes = [];
    end

    if isempty(spikeTimes)
        hFigPreview.hidePlot('hPlot_traces_spk1');
    else
        hFigPreview.multiplot('hPlot_traces_spk1', figData.maxAmp, spikeTimesSec, jrclust.utils.rowColSelect(traces, spikeTimes, spikeSites), spikeSites, 1);
    end

    hCfg.setTemporaryParams('filterType', figData.filterType);
    siteThreshuV = jrclust.utils.bit2uV(-figData.siteThresh(:), hCfg);
    hCfg.resetTemporaryParams();
    siteThreshuV(figData.ignoreSites) = nan;

    if figData.fShowThresh && figData.fFilter
        hFigPreview.multiplot('hPlotTracesThresh', figData.maxAmp, XData([1,end,end])', repmat(siteThreshuV, [1, 3])');
        hFigPreview.multiplot('hPlot_traces_spk', figData.maxAmp, XData, ...
            mrSet(traces, ~figData.isThreshCrossing(XDataSamp, :), nan)); % show spikes
    else
        hFigPreview.hidePlot('hPlotTracesThresh');
        hFigPreview.hidePlot('hPlot_traces_spk');
    end

    hFigPreview.axApply('hAxTraces', @ylabel, 'Site #');
    filterLabel = jrclust.utils.ifEq(figData.fFilter, sprintf('Filter=%s', figData.filterType), 'Filter off');
    hFigPreview.figApply(@set, 'Name', sprintf('%s; %s; CommonRef=%s', hCfg.configFile, filterLabel, figData.CARMode));

    hFigPreview.axApply('hAxTraces', @title, sprintf('Scale: %0.1f uV', figData.maxAmp));

    if ~fKeepView
        hFigPreview.axApply('hAxTraces', @set, 'YTick', 1:nSites, 'YLim', figData.siteLim + [-1, 1], 'XLim', tlimSecs);
    end
    hFigPreview.axApply('hAxTraces', @grid, jrclust.utils.ifEq(figData.fGrid, 'on', 'off'));

    %% Site plot
    if ~hFigPreview.hasPlot('hPlotSite')
        hFigPreview.addPlot('hPlotSite', @barh, hFigPreview.hAxes('hAxSites'), nan, nan, 1);
        hFigPreview.addPlot('hPlotSiteBad', @barh, hFigPreview.hAxes('hAxSites'), nan, nan, 1, 'r');
        hFigPreview.addPlot('hPlotSiteThresh', @plot, hFigPreview.hAxes('hAxSites'), nan, nan, 'r');
    end

    switch figData.siteView
        case 'Site correlation'
            YData = figData.maxCorrSite;

        case 'Spike threshold'
            YData = single(figData.siteThresh);

        case 'Event rate (Hz)'
            YData = figData.siteEventRate;

        case 'Event SNR (median)'
            YData = figData.siteEventSNR;
    end % switch

    hFigPreview.plotApply('hPlotSite', @set, 'XData', 1:nSites, 'YData', YData);
    hFigPreview.axApply('hAxSites', @set, 'YLim', figData.siteLim + [-1,1]);
    hFigPreview.axApply('hAxPSD', @grid, jrclust.utils.ifEq(figData.fGrid, 'on', 'off'));

    if isempty(figData.siteCorrThresh) || ~strcmpi(figData.siteView, 'Site correlation')
        hFigPreview.hidePlot('hPlotSiteThresh');
    else
        hFigPreview.plotApply('hPlotSiteThresh', @set, 'XData', figData.siteCorrThresh *[1,1], 'YData', [0, nSites+1]);
    end
    if ~isempty(figData.ignoreSites)
        YData(~figData.ignoreMe) = 0;
        hFigPreview.plotApply('hPlotSiteBad', @set, 'XData', 1:nSites, 'YData', YData); %switch statement
    else
        hFigPreview.hidePlot('hPlotSiteBad');
    end

    hFigPreview.axApply('hAxSites', @xlabel, figData.siteView);
    hFigPreview.axApply('hAxSites', @ylabel, 'Site #');
    hFigPreview.axApply('hAxSites', @title, sprintf('siteCorrThresh=%0.4f', figData.siteCorrThresh));

    %% PSD plot
    if ~hFigPreview.hasPlot('hPlotPSD')
        hFigPreview.addPlot('hPlotPSD', @plot, hFigPreview.hAxes('hAxPSD'), nan, nan, 'k');
        hFigPreview.addPlot('hPlotCleanPSD', @plot, hFigPreview.hAxes('hAxPSD'), nan, nan, 'g');
        hFigPreview.addPlot('hPlotPSDThresh', @plot, hFigPreview.hAxes('hAxPSD'), nan, nan, 'r');
    end

    hFigPreview.plotApply('hPlotPSD', @set, 'XData', figData.psdFreq, 'YData', figData.psdPower);
    hFigPreview.plotApply('hPlotCleanPSD', @set, 'XData', figData.psdFreq, 'YData', figData.psdPowerClean);

    hFigPreview.axApply('hAxPSD', @set, 'XLim', [0, hCfg.sampleRate/2]);
    hFigPreview.axApply('hAxPSD', @grid, jrclust.utils.ifEq(figData.fGrid, 'on', 'off'));

    hFigPreview.axApply('hAxPSD', @xlabel, 'Frequency (Hz)');
    hFigPreview.axApply('hAxPSD', @ylabel, 'Power [dB]');
    hFigPreview.axApply('hAxPSD', @title, sprintf('fftThresh=%s', num2str(figData.fftThresh)));


    %% finish up
    hFigPreview.wait(0);
    jrclust.utils.tryClose(hWait);
end

%% LOCAL FUNCTIONS
function mr = mrSet(mr, ml, val)
    mr(ml) = val;
end
