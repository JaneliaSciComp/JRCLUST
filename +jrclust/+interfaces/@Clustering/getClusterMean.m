function [clusterMean, siteNeighbors, clusterMeanLow, clusterMeanHigh] = getClusterMean(obj, spikeWindows, iCluster)
    %GETCLUSTERMEAN Get mean cluster waveform, optionally low and high
    [clusterMean, clusterMeanLow, clusterMeanHigh] = deal([]);
    iSite = obj.clusterSites(iCluster);
    siteNeighbors = obj.hCfg.siteNeighbors(:, iSite);

    clusterSpikes = obj.spikesByCluster{iCluster};
    clusterSites_ = obj.spikeSites(clusterSpikes);

    if isempty(clusterSpikes)
        return;
    end

    if ~obj.hCfg.driftMerge || isempty(obj.spikePositions)
        middlemost = getMidmostSpikes(clusterSpikes, obj.spikeTimes, obj.hCfg.nClusterIntervals);
        clusterMean = mean(single(spikeWindows(:, :, middlemost)), 3);
        clusterMean = jrclust.utils.meanSubtract(clusterMean);
        return;
    end

    yPos = obj.spikePositions(clusterSpikes, 2); % position based quantile
    yLims = quantile(yPos, [0, 1, 2, 3]/3);

    [selectedSpikes, selectedSites] = getSpikesInBounds(clusterSpikes, clusterSites_, yPos, yLims(2:3));
    clusterMean = nanmeanInt16(spikeWindows(:, :, selectedSpikes), iSite, selectedSites, obj.hCfg); % * hCfg.uV_per_bit;

    % get mean waveforms at high and low positions on the probe (raw traces
    % only)
    if nargout > 2
        [selectedSpikes, selectedSites] = getSpikesInBounds(clusterSpikes, clusterSites_, yPos, yLims(1:2));
        clusterMeanLow = nanmeanInt16(spikeWindows(:, :, selectedSpikes), iSite, selectedSites, obj.hCfg);

        [selectedSpikes, selectedSites] = getSpikesInBounds(clusterSpikes, clusterSites_, yPos, yLims(3:4));
        clusterMeanHigh = nanmeanInt16(spikeWindows(:, :, selectedSpikes), iSite, selectedSites, obj.hCfg);
    end
end
