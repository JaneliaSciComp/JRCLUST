function hFigPos = doPlotFigPos(hFigPos, hClust, hCfg, selected, maxAmp)
    %DOPLOTFIGPOS Plot position of cluster on probe
    c1Data = hClust.exportUnitInfo(selected(1));
    if numel(selected) > 1
        c2Data = hClust.exportUnitInfo(selected(2));
    end

    if isempty(hFigPos.figData)
        hFigPos.axes();
    else
        hFigPos.cla();
        hFigPos.hold('on');
    end

    plotPosUnit(c1Data, hFigPos, hCfg, false, maxAmp);

    clusterPos = hClust.clusterCentroids(c1Data.cluster, :)/hCfg.um_per_pix;
    nSpikes = hClust.clusterCounts(c1Data.cluster);

    if numel(selected) == 1
        figTitle = sprintf('Unit %d: %d spikes; (X=%0.1f, Y=%0.1f) [um]', c1Data.cluster, nSpikes, clusterPos);
        try
            figTitle = sprintf('%s\n%0.1fuVmin, %0.1fuVpp, SNR:%0.1f ISI%%:%2.3f IsoDist:%0.1f L-rat:%0.1f', ...
                figTitle, c1Data.peaksRaw, c1Data.vpp, c1Data.SNR, c1Data.ISIRatio*100, c1Data.IsoDist,  c1Data.LRatio);
        catch
        end
    else
        nSpikes2 = hClust.clusterCounts(c2Data.cluster);
        clusterPos2 = hClust.clusterCentroids(c2Data.cluster, :)/hCfg.um_per_pix;
        plotPosUnit(c2Data, hFigPos, hCfg, true, maxAmp);

        figTitle = sprintf('Unit %d(black)/%d(red); (%d/%d) spikes\n(X=%0.1f/%0.1f, Y=%0.1f/%0.1f) [um]', ...
            c1Data.cluster, c2Data.cluster, nSpikes, nSpikes2, ...
            [clusterPos(1), clusterPos2(1), clusterPos(2), clusterPos2(2)]);
    end
    hFigPos.title(figTitle);
end

%% LOCAL FUNCTIONS
function plotPosUnit(cData, hFigPos, hCfg, fSecondary, maxAmp)
    if isempty(cData)
        return;
    end
    if fSecondary
        cmapMean = hCfg.mrColor_proj(3, :); % red
    else
        cmapMean = hCfg.mrColor_proj(2, :); % black
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

        siteXData = hCfg.siteLoc(cData.neighbors, 1) / hCfg.um_per_pix;
        siteYData = hCfg.siteLoc(cData.neighbors, 2) / hCfg.um_per_pix;
        YData = bsxfun(@plus, YData, siteYData');
        XData = bsxfun(@plus, repmat(XBase, [1, size(YData, 2)]), siteXData');
        hFigPos.addLine(sprintf('neighbor%d', iWav), XData(:), YData(:), ...
                        'Color', cmap, 'LineWidth', lineWidth);
    end

    hFigPos.xlabel('X pos [pix]');
    hFigPos.ylabel('Z pos [pix]');
    hFigPos.grid('on');
    hFigPos.axSet('XLim', [min(XBase(:)), max(XBase(:))]);
    hFigPos.axSet('YLim', [floor(min(YData(:))-1), ceil(max(YData(:))+1)]);
end
