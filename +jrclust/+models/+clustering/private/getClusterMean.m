function [clusterMean, siteNeighbors, clusterMeanLow, clusterMeanHigh] = getClusterMean(hClust, spikeWindows, iCluster)
    %GETCLUSTERMEAN Get mean cluster waveform, optionally low and high
    [clusterMean, clusterMeanLow, clusterMeanHigh] = deal([]);
    iSite = hClust.clusterSites(iCluster);
    siteNeighbors = hClust.hCfg.siteNeighbors(:, iSite);

    clusterSpikes = hClust.spikesByCluster{iCluster};
    clusterSites_ = hClust.spikeSites(clusterSpikes);

    if isempty(clusterSpikes)
        return;
    end

    if ~hClust.hCfg.driftMerge || isempty(hClust.spikePositions)
        middlemost = spk_select_mid_(clusterSpikes, hClust.spikeTimes, hClust.hCfg.nClusterIntervals);
        clusterMean = mean(single(spikeWindows(:, :, middlemost)), 3);
        clusterMean = jrclust.utils.meanSubtract(clusterMean);
        return;
    end

    yPos = hClust.spikePositions(clusterSpikes, 2); % position based quantile
    yLims = quantile(yPos, [0, 1, 2, 3]/3);

    [selectedSpikes, selectedSites] = selectSpikesInBounds(clusterSpikes, clusterSites_, yPos, yLims(2:3));
    clusterMean = nanmeanInt16(spikeWindows(:, :, selectedSpikes), iSite, selectedSites, hClust.hCfg); % * hCfg.uV_per_bit;

    if nargout > 2
        [selectedSpikes, selectedSites] = selectSpikesInBounds(clusterSpikes, clusterSites_, yPos, yLims(1:2));
        clusterMeanLow = nanmeanInt16(spikeWindows(:, :, selectedSpikes), iSite, selectedSites, hClust.hCfg);

        [selectedSpikes, selectedSites] = selectSpikesInBounds(clusterSpikes, clusterSites_, yPos, yLims(3:4));
        clusterMeanHigh = nanmeanInt16(spikeWindows(:, :, selectedSpikes), iSite, selectedSites, hClust.hCfg);
    end
end

%% LOCAL FUNCTIONS
function [spikesOut, sitesOut] = selectSpikesInBounds(spikesIn, sitesIn, yPos, yLims)
    %SELECTSPIKESINBOUNDS Select spikes by whether their y-positions fall within given limits
    nSamplesMax = 1000;
    inBounds = yPos >= yLims(1) & yPos < yLims(2);

    if ~any(inBounds)
        spikesOut = spikesIn;
        sitesOut = sitesIn;
        return;
    end

    spikesOut = jrclust.utils.subsample(spikesIn(inBounds), nSamplesMax);
    sitesOut = jrclust.utils.subsample(sitesIn(inBounds), nSamplesMax);
end

function mrWav_clu1 = nanmeanInt16(spikeWindows, iSite, sites, hCfg)
    iSiteNeighbors = hCfg.siteNeighbors(:, iSite);
    trWav = nan([size(spikeWindows, 1), numel(iSiteNeighbors), numel(sites)], 'single');
    uniqueSites = unique(sites);
    nUniqueSites = numel(uniqueSites);
    uniqueNeighbors = hCfg.siteNeighbors(:, uniqueSites);

    for jSite = 1:nUniqueSites
        iSiteUnique = uniqueSites(jSite);
        viSpk_ = find(sites == iSiteUnique);

        [~, viSite1a_, viSite1b_] = intersect(iSiteNeighbors, uniqueNeighbors(:, jSite));
        if isempty(viSite1a_)
            continue;
        end

        trWav(:, viSite1a_, viSpk_) = spikeWindows(:, viSite1b_, viSpk_);
    end

    mrWav_clu1 = nanmean(trWav, 3);
    mrWav_clu1 = jrclust.utils.meanSubtract(mrWav_clu1); %122717 JJJ
end

function subSpikes = spk_select_mid_(spikes, spikeTimes, nClusterIntervals)
    % viTime_spk = get0_('viTime_spk');
    iSpikeMid = round(numel(spikeTimes)/2); % index of the middlest spike
    nearestToCenter = jrclust.utils.rankorder(abs(spikes - iSpikeMid), 'ascend');
    nSpikesInterval = round(numel(spikes) / nClusterIntervals);
    subSpikes = spikes(nearestToCenter <= nSpikesInterval);
end
