function results = doComputeMeanWaveforms(hClust, updateMe, useRaw)
    %DOCOMPUTEMEANWAVEFORMS Compute the mean waveform for each cluster
    results = struct();
    if isempty(hClust.spikesFilt)
        return;
    end

    useRaw = useRaw && ~isempty(hClust.spikesRaw);
    verbose = hClust.hCfg.verbose && isempty(updateMe);
    if verbose
        fprintf('Calculating cluster mean waveform.\n\t');
        t = tic;
    end

    nSites = numel(hClust.hCfg.siteMap);
    [nSamples, nSitesEvt, ~] = size(hClust.spikesFilt);

    meanWfLocal_ = zeros(nSamples, nSitesEvt, hClust.nClusters, 'single');
    meanWfGlobal_ = zeros(nSamples, nSites, hClust.nClusters, 'single');

    if useRaw
        nSamplesRaw = size(hClust.spikesRaw, 1);
        meanWfLocalRaw_ = zeros(nSamplesRaw, nSitesEvt, hClust.nClusters, 'single');
        meanWfGlobalRaw_ = zeros(nSamplesRaw, nSites, hClust.nClusters, 'single');
        meanWfRawLow_ = zeros(nSamplesRaw, nSites, hClust.nClusters, 'single');
        meanWfRawHigh_ = zeros(nSamplesRaw, nSites, hClust.nClusters, 'single');
    else
        [meanWfLocalRaw_, meanWfGlobalRaw_, meanWfRawLow_, meanWfRawHigh_] = deal([]);
    end

    % we have specific clusters to update
    if ~isempty(updateMe) && ~isempty(hClust.meanWfLocal)
        % visit all clusters explicitly requested or not previously seen
        visitMe = false(hClust.nClusters, 1);
        visitMe(updateMe) = 1;
        visitMe((1:hClust.nClusters) > size(hClust.meanWfLocal, 3)) = 1;

        meanWfLocal_ = hClust.meanWfLocal;
        meanWfGlobal_ = hClust.meanWfGlobal;

        meanWfLocalRaw_ = hClust.meanWfLocalRaw;
        meanWfGlobalRaw_ = hClust.meanWfGlobalRaw;

        meanWfRawLow_ = hClust.tmrWav_raw_lo_clu;
        meanWfRawHigh_ = hClust.tmrWav_raw_hi_clu;
    else % nothing specific requested, update all
        visitMe = true(hClust.nClusters, 1);
    end

    for iCluster = 1:hClust.nClusters
        if visitMe(iCluster)
            [clusterMean, siteNeighbors] = getClusterMean(hClust, hClust.spikesFilt, iCluster);

            if isempty(clusterMean)
                continue;
            end

            clusterMean = jrclust.utils.bit2uV(clusterMean, hClust.hCfg);
            meanWfLocal_(:, :, iCluster) = clusterMean;
            meanWfGlobal_(:, siteNeighbors, iCluster) = clusterMean;
        end

        if verbose
            fprintf('.');
        end
    end

    if useRaw
        for iCluster = 1:hClust.nClusters
            if visitMe(iCluster)
                [clusterMean, siteNeighbors, clusterMeanLow, clusterMeanHigh] = getClusterMean(hClust, hClust.spikesRaw, iCluster);
                if isempty(clusterMean)
                    continue;
                end

                clusterMean = jrclust.utils.meanSubtract(clusterMean)*hClust.hCfg.bitScaling;
                meanWfGlobalRaw_(:, siteNeighbors, iCluster) = clusterMean;
                meanWfLocalRaw_(:, :, iCluster) = clusterMean;

                if isempty(clusterMeanLow) || isempty(clusterMeanHigh)
                    meanWfRawLow_(:, siteNeighbors, iCluster) = zeros(nSamplesRaw, numel(siteNeighbors));
                    meanWfRawHigh_(:, siteNeighbors, iCluster) = zeros(nSamplesRaw, numel(siteNeighbors));
                else
                    meanWfRawLow_(:,siteNeighbors,iCluster) = jrclust.utils.meanSubtract(clusterMeanLow)*hClust.hCfg.bitScaling;
                    meanWfRawHigh_(:,siteNeighbors,iCluster) = jrclust.utils.meanSubtract(clusterMeanHigh)*hClust.hCfg.bitScaling;
                end
            end

            if verbose
                fprintf('.');
            end
        end
    end

    % measure waveforms
    [unitPeaks_, unitPeakSites_] = min(permute(min(meanWfLocal_), [2, 3, 1]), [], 1);
    results.unitPeaks = abs(unitPeaks_(:));
    results.unitPeakSites = unitPeakSites_(:);

    % collect computed values
    results.meanWfLocal = meanWfLocal_;
    results.meanWfGlobal = meanWfGlobal_;
    results.meanWfLocalRaw = meanWfLocalRaw_;
    results.meanWfGlobalRaw = meanWfGlobalRaw_;
    results.meanWfRawLow = meanWfRawLow_;
    results.meanWfRawHigh = meanWfRawHigh_;

    if verbose
        fprintf('\n\ttook %0.1fs\n', toc(t));
    end
end