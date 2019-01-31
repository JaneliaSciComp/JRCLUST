function hFigPos = doPlotFigPos(hFigPos, hClust, hCfg, selected, maxAmp)
    %DOPLOTFIGPOS Plot position of cluster on probe
    c1Data = hClust.exportUnitInfo(selected(1));
    if numel(selected) > 1
        c2Data = hClust.exportUnitInfo(selected(2));
    end

    if ~hFigPos.hasAxes('default')
        hFigPos.addAxes('default');
    else
        hFigPos.axApply('default', @cla);
    end

    plotPosUnit(c1Data, hFigPos, hCfg, 0, maxAmp);

    clusterPos = hClust.clusterCentroids(c1Data.cluster, :)/hCfg.umPerPix;
    nSpikes = hClust.unitCount(c1Data.cluster);

    if numel(selected) == 1
        figTitle = sprintf('Unit %d: %d spikes; (X=%0.1f, Y=%0.1f) [um]', c1Data.cluster, nSpikes, clusterPos);
        try
            figTitle = sprintf('%s\n%0.1fuVmin, %0.1fuVpp, SNR:%0.1f ISI%%:%2.3f IsoDist:%0.1f L-rat:%0.1f', ...
                figTitle, c1Data.peaksRaw, c1Data.vpp, c1Data.SNR, c1Data.ISIRatio*100, c1Data.IsoDist,  c1Data.LRatio);
        catch
        end
    else
        nSpikes2 = hClust.unitCount(c2Data.cluster);
        clusterPos2 = hClust.clusterCentroids(c2Data.cluster, :)/hCfg.umPerPix;
        plotPosUnit(c2Data, hFigPos, hCfg, 1, maxAmp);

        figTitle = sprintf('Unit %d(black)/%d(red); (%d/%d) spikes\n(X=%0.1f/%0.1f, Y=%0.1f/%0.1f) [um]', ...
            c1Data.cluster, c2Data.cluster, nSpikes, nSpikes2, ...
            [clusterPos(1), clusterPos2(1), clusterPos(2), clusterPos2(2)]);
    end
    hFigPos.axApply('default', @title, figTitle);
end

%% LOCAL FUNCTIONS
function plotPosUnit(cData, hFigPos, hCfg, fSecondary, maxAmp)
    if isempty(cData)
        return;
    end
    if fSecondary
        cmapMean = hCfg.colorMap(3, :); % red
    else
        cmapMean = hCfg.colorMap(2, :); % black
    end

    % plot individual unit
    nSamples = size(cData.meanWf, 1);
    XBase = (1:nSamples)'/nSamples;
    XBase([1, end]) = nan; % line break

    if fSecondary
        sampleWf = zeros(1, 1, 0);
    else
        sampleWf = cData.sampleWf;
    end

    siteXData = hCfg.siteLoc(cData.neighbors, 1) / hCfg.umPerPix;
    siteYData = hCfg.siteLoc(cData.neighbors, 2) / hCfg.umPerPix;

    % show example traces
    for iWav = size(sampleWf, 3):-1:0
        if iWav == 0 % plot the cluster mean waveform
            YData = cData.meanWf/maxAmp;
            lineWidth = 1.5;
            cmap = cmapMean;
        else
            YData = sampleWf(:,:,iWav) / maxAmp;
            lineWidth = 0.5;
            cmap = 0.5*[1, 1, 1];
        end

        YData = bsxfun(@plus, YData, siteYData');
        XData = bsxfun(@plus, repmat(XBase, [1, size(YData, 2)]), siteXData');
        hFigPos.addPlot(sprintf('neighbor%d%d', double(fSecondary) + 1, iWav), @line, XData(:), YData(:), ...
                        'Color', cmap, 'LineWidth', lineWidth);
    end

    hFigPos.axApply('default', @xlabel, 'X pos [pix]');
    hFigPos.axApply('default', @ylabel, 'Z pos [pix]');
    hFigPos.axApply('default', @grid, 'on');
    hFigPos.axApply('default', @set, 'XLim', [min(XBase(:)), max(XBase(:))] + median(siteXData));
    hFigPos.axApply('default', @set, 'YLim', [floor(min(YData(:))-1), ceil(max(YData(:))+1)]);
end
