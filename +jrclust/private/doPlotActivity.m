function doPlotActivity(spikeData, hCfg)
    %DOPLOTACTIVITY Plot activity as a function of depth and time
    tBin = 10; % activity every 10 sec
    % plot activity as a function of time
    recDur = double(max(spikeData.spikeTimes)) / hCfg.sampleRate; % in sec
    nBins = ceil(recDur / tBin);

    % 90th percentile of amplitudes per time bin per site
    amp90 = zeros(nBins, hCfg.nSites);

    binLims = [1, tBin * hCfg.sampleRate];
    for iSite = 1:hCfg.nSites
        iSpikes = find(spikeData.spikeSites == iSite);
        iAmps = spikeData.spikeAmps(iSpikes);
        if isempty(iAmps)
            continue;
        end

        iTimes = spikeData.spikeTimes(iSpikes);
        for iBin = 1:nBins
            bounds = binLims + (iBin-1) * binLims(2);
            inBounds = iTimes >= bounds(1) & iTimes <= bounds(2);
            if ~any(inBounds)
                continue;
            end

            binAmps = iAmps(inBounds);
            amp90(iBin, iSite) = quantile(abs(binAmps), .9);
        end
    end % for
    amp90 = amp90';

    leftEdge = hCfg.siteLoc(:, 1) == 0; % the leftmost column on the probe
    YLocs = hCfg.siteLoc(:, 2);

    hFigActivity = jrclust.views.Figure('FigActivity', [0 0 .5 1], hCfg.sessionName, 1, 1);
    hFigActivity.addSubplot('hActivity', 1, 2);

    hFigActivity.subplotApply('hActivity', 1, @imagesc, amp90(leftEdge, :), 'XData', (1:nBins) * tBin, 'YData', YLocs(leftEdge));
    hFigActivity.subplotApply('hActivity', 1, @axis, 'xy');
    hFigActivity.subplotApply('hActivity', 1, @xlabel, 'Time');
    hFigActivity.subplotApply('hActivity', 1, @ylabel, 'Sites');
    hFigActivity.subplotApply('hActivity', 1, @title, 'Left edge sites');
    hFigActivity.subplotApply('hActivity', 2, @imagesc, amp90(~leftEdge, :), 'XData', (1:nBins) * tBin, 'YData', YLocs(leftEdge));
    hFigActivity.subplotApply('hActivity', 2, @axis, 'xy');
    hFigActivity.subplotApply('hActivity', 2, @xlabel, 'Time');
    hFigActivity.subplotApply('hActivity', 2, @ylabel, 'Sites');
    hFigActivity.subplotApply('hActivity', 2, @title, 'Right edge sites');

    [~, meanSite] = max(mean(amp90, 2));
    meanSiteNeighbors = intersect(1:max(meanSite)+2, meanSite + (-2:2));
    amp90Center = amp90(meanSiteNeighbors, :);
    centroid = bsxfun(@rdivide, sum(bsxfun(@times, amp90Center.^2, YLocs(meanSiteNeighbors))), sum(amp90Center.^2));
    if numel(centroid) == 1
        LineStyle = 'r*';
    else
        LineStyle = 'r-';
    end
    hFigActivity.subplotApply('hActivity', 1, @hold, 'on');
    hFigActivity.subplotApply('hActivity', 1, @plot, (1:nBins) * tBin, centroid, LineStyle);
    hFigActivity.subplotApply('hActivity', 2, @hold, 'on');
    hFigActivity.subplotApply('hActivity', 2, @plot, (1:nBins) * tBin, centroid, LineStyle);
end
