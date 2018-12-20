function hFigWav = plotMeanWaveforms(hFigWav, hClust, hCfg, maxAmp)
    %PLOTMEANWAVEFORMS Plot mean cluster waveforms in the main view
    if hCfg.showRaw
        waveforms = hClust.meanWfGlobalRaw;
    else
        waveforms = hClust.meanWfGlobal;
    end
    
    [nSamples, nSites, nClusters] = size(waveforms);
    nSitesShow = size(hCfg.siteNeighbors, 1);

    % determine x
    xOffset = hCfg.evtWindowSamp(2)/(diff(hCfg.evtWindowSamp) + 1); % same for raw and filt
    xData = (1:nSamples*nClusters)/nSamples + xOffset;
    xData(1:nSamples:end) = nan;
    xData(nSamples:nSamples:end) = nan;
    waveforms = waveforms/maxAmp;

    xData = repmat(xData(:), [1, nSitesShow]);
    xData = reshape(xData, [nSamples, nClusters, nSitesShow]);
    xData = reshape(permute(xData, [1 3 2]), [nSamples*nSitesShow, nClusters]);

    yData = zeros(nSamples * nSitesShow, nClusters, 'single');
    for iCluster = 1:nClusters
        iSites = hCfg.siteNeighbors(:, hClust.clusterSites(iCluster))   ;
        iWaveforms = waveforms(:, iSites, iCluster);
        iWaveforms = bsxfun(@plus, iWaveforms, single(iSites'));
        yData(:, iCluster) = iWaveforms(:);
    end

    if ~hFigWav.hasPlot('Group1')
        plotGroup(hFigWav, xData, yData, 'LineWidth', hCfg.getOr('LineWidth', 1));
    else
        updateGroup(hFigWav, xData, yData);
    end

    hFigWav.axSet('YTick', 1:nSites, 'XTick', 1:nClusters);
end

%% LOCAL FUNCTIONS
function updateGroup(hFig, xData, yData)
    %UPDATE Update group-plotted data
    iGroup = 1;

    while hFig.hasPlot(sprintf('hGroup%d', iGroup))
        iXData = xData(:, iPlot:nPlots:end);
        iYData = yData(:, iPlot:nPlots:end);
        hFig.updatePlot(sprintf('hGroup%d', iGroup), iXData, iYData);

        iGroup = iGroup + 1;
    end
end

function plotGroup(hFig, xData, yData, varargin)
    %PLOTGROUP Plot xData and yData colored by groups
    colorMap = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]';
    nGroups = min(size(colorMap, 2), size(xData, 2));
    colorMap = colorMap(:, 1:nGroups);

    hFig.hold('on');
    for iGroup = 1:nGroups
        iXData = xData(:, iGroup:nGroups:end);
        iYData = yData(:, iGroup:nGroups:end);
        hFig.addPlot(sprintf('hGroup%d', iGroup), iXData(:), iYData(:), varargin{:}, 'Color', colorMap(:, iGroup)');
    end
end
