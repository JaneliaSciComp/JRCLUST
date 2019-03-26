function hFigWav = plotMeanWaveforms(hFigWav, hClust, hCfg, maxAmp)
    %PLOTMEANWAVEFORMS Plot mean cluster waveforms in FigWav
    if hCfg.showRaw
        waveforms = hClust.meanWfGlobalRaw;
    else
        waveforms = hClust.meanWfGlobal;
    end

    [nSamples, nSites, nClusters] = size(waveforms);
    nSitesShow = size(hCfg.siteNeighbors, 1);

    % determine x
    if hCfg.showRaw
        xOffset = hCfg.evtWindowRawSamp(2)/(diff(hCfg.evtWindowRawSamp) + 1);
    else
        xOffset = hCfg.evtWindowSamp(2)/(diff(hCfg.evtWindowSamp) + 1);
    end
    XData = (1:nSamples*nClusters)/nSamples + xOffset;

    % breaks between clusters
    XData(1:nSamples:end) = nan;
    XData(nSamples:nSamples:end) = nan;
    waveforms = waveforms/maxAmp;

    XData = repmat(XData(:), [1, nSitesShow]);
    XData = reshape(XData, [nSamples, nClusters, nSitesShow]);
    XData = reshape(permute(XData, [1 3 2]), [nSamples*nSitesShow, nClusters]);

    YData = zeros(nSamples * nSitesShow, nClusters, 'single');
    for iCluster = 1:nClusters
        iSites = hCfg.siteNeighbors(:, hClust.clusterSites(iCluster));
        iWaveforms = waveforms(:, iSites, iCluster);
        iWaveforms = bsxfun(@plus, iWaveforms, single(iSites'));
        YData(:, iCluster) = iWaveforms(:);
    end

    if ~hFigWav.hasPlot('hGroup1')
        plotGroup(hFigWav, XData, YData, 'LineWidth', hCfg.getOr('LineWidth', 1));
    else
        %updateGroup(hFigWav, XData, YData); % this is broken
        plotKeys = keys(hFigWav.hPlots);
        hGroup = plotKeys(startsWith(plotKeys, 'hGroup'));
        if ~isempty(hGroup)
            cellfun(@(plotKey) hFigWav.rmPlot(plotKey), hGroup);
        end
        plotGroup(hFigWav, XData, YData, 'LineWidth', hCfg.getOr('LineWidth', 1));
    end

    hFigWav.axApply('default', @set, 'YTick', 1:nSites, 'XTick', 1:nClusters);
end

%% LOCAL FUNCTIONS
function updateGroup(hFig, XData, YData)
    %UPDATE Update group-plotted data
    nGroups = sum(cellfun(@(c) ~isempty(c), regexp(keys(hFig.hPlots), '^hGroup\d')));

    for iGroup = 1:numel(nGroups)
        iXData = XData(:, iGroup:nGroups:end);
        iYData = YData(:, iGroup:nGroups:end);
        hFig.updatePlot(sprintf('hGroup%d', iGroup), iXData(:), iYData(:));
    end
end

function plotGroup(hFig, XData, YData, varargin)
    %PLOTGROUP Plot XData and YData colored by groups
    colorMap = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]';
    nGroups = min(size(colorMap, 2), size(XData, 2));
    colorMap = colorMap(:, 1:nGroups);

    hFig.axApply('default', @hold, 'on');
    for iGroup = 1:nGroups
        iXData = XData(:, iGroup:nGroups:end);
        iYData = YData(:, iGroup:nGroups:end);
        hFig.addPlot(sprintf('hGroup%d', iGroup), iXData(:), iYData(:), varargin{:}, 'Color', colorMap(:, iGroup)');
    end
end
