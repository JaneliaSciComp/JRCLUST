function computeMeanWaveforms(obj, updateMe, useRaw)
    %COMPUTEMEANWAVEFORMS Compute the mean waveform for each cluster
    if nargin < 2
        updateMe = [];
    end
    if nargin < 3
        useRaw = 1;
    end

    if isempty(obj.spikesFilt)
        return;
    end

    useRaw = useRaw && ~isempty(obj.spikesRaw);
    obj.hCfg.updateLog('meanWf', 'Computing cluster mean waveforms', 1, 0);

    [nSamples, nSitesEvt, ~] = size(obj.spikesFilt);

    meanWfLocal_ = zeros(nSamples, nSitesEvt, obj.nClusters, 'single');
    meanWfGlobal_ = zeros(nSamples, obj.hCfg.nSites, obj.nClusters, 'single');

    if useRaw
        nSamplesRaw = size(obj.spikesRaw, 1);
        meanWfLocalRaw_ = zeros(nSamplesRaw, nSitesEvt, obj.nClusters, 'single');
        meanWfGlobalRaw_ = zeros(nSamplesRaw, obj.hCfg.nSites, obj.nClusters, 'single');
        meanWfRawLow_ = zeros(nSamplesRaw, obj.hCfg.nSites, obj.nClusters, 'single');
        meanWfRawHigh_ = zeros(nSamplesRaw, obj.hCfg.nSites, obj.nClusters, 'single');
    else
        [meanWfLocalRaw_, meanWfGlobalRaw_, meanWfRawLow_, meanWfRawHigh_] = deal([]);
    end

    % we have specific clusters to update
    if ~isempty(updateMe) && ~isempty(obj.meanWfLocal)
        % visit all clusters explicitly requested or not previously seen
        visitMe = false(obj.nClusters, 1);
        visitMe(updateMe) = 1;
        visitMe((1:obj.nClusters) > size(obj.meanWfLocal, 3)) = 1;

        meanWfLocal_ = obj.meanWfLocal;
        meanWfGlobal_ = obj.meanWfGlobal;

        meanWfLocalRaw_ = obj.meanWfLocalRaw;
        meanWfGlobalRaw_ = obj.meanWfGlobalRaw;

        meanWfRawLow_ = obj.meanWfRawLow;
        meanWfRawHigh_ = obj.meanWfRawHigh;
    else % nothing specific requested, update all
        visitMe = true(obj.nClusters, 1);
    end

    for iCluster = 1:obj.nClusters
        if visitMe(iCluster)
            [clusterMean, siteNeighbors] = obj.getClusterMean(obj.spikesFilt, iCluster);

            if isempty(clusterMean)
                continue;
            end

            clusterMean = jrclust.utils.bit2uV(clusterMean, obj.hCfg);
            meanWfLocal_(:, :, iCluster) = clusterMean;
            meanWfGlobal_(:, siteNeighbors, iCluster) = clusterMean;
        end
    end

    if useRaw
        for iCluster = 1:obj.nClusters
            if visitMe(iCluster)
                [clusterMean, siteNeighbors, clusterMeanLow, clusterMeanHigh] = obj.getClusterMean(obj.spikesRaw, iCluster);
                if isempty(clusterMean)
                    continue;
                end

                clusterMean = jrclust.utils.meanSubtract(clusterMean)*obj.hCfg.bitScaling;
                meanWfGlobalRaw_(:, siteNeighbors, iCluster) = clusterMean;
                meanWfLocalRaw_(:, :, iCluster) = clusterMean;

                if isempty(clusterMeanLow) || isempty(clusterMeanHigh)
                    meanWfRawLow_(:, siteNeighbors, iCluster) = zeros(nSamplesRaw, numel(siteNeighbors));
                    meanWfRawHigh_(:, siteNeighbors, iCluster) = zeros(nSamplesRaw, numel(siteNeighbors));
                else
                    meanWfRawLow_(:,siteNeighbors,iCluster) = jrclust.utils.meanSubtract(clusterMeanLow)*obj.hCfg.bitScaling;
                    meanWfRawHigh_(:,siteNeighbors,iCluster) = jrclust.utils.meanSubtract(clusterMeanHigh)*obj.hCfg.bitScaling;
                end
            end
        end
    end

    % measure waveforms
    [unitPeaks_, unitPeakSites_] = min(permute(min(meanWfLocal_), [2, 3, 1]), [], 1);
    obj.unitPeaks = abs(unitPeaks_(:));
    obj.unitPeakSites = unitPeakSites_(:);

    % collect computed values
    obj.meanWfLocal = meanWfLocal_;
    obj.meanWfGlobal = meanWfGlobal_;
    obj.meanWfLocalRaw = meanWfLocalRaw_;
    obj.meanWfGlobalRaw = meanWfGlobalRaw_;
    obj.meanWfRawLow = meanWfRawLow_;
    obj.meanWfRawHigh = meanWfRawHigh_;

    obj.hCfg.updateLog('meanWf', 'Finished computing cluster mean waveforms', 0, 1);
end