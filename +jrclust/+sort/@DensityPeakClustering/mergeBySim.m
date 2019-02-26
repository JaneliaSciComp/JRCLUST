function nMerged = mergeBySim(obj)
    %MERGEBYSIM Automatically merge clusters by similarity score
    obj.hCfg.updateLog('mergeBySim', 'Merging clusters by similarity', 1, 0);
    waveformSim_ = obj.waveformSim;
    nClustersPrev = size(waveformSim_, 2);

    % identify clusters to remove, update and same (no change), disjoint sets
    waveformSim_(tril(true(nClustersPrev))) = 0; % ignore bottom half
    [scoresMax, mapTo] = max(waveformSim_);

    % keep clusters whose maximum similarity to some other cluster
    % is LESS than our threshold
    keepMe_ = scoresMax < obj.hCfg.maxUnitSim;

    if all(keepMe_)
        nMerged = 0;
        obj.hCfg.updateLog('mergeBySim', 'No clusters to merge', 0, 1);
        return;
    end

    minScore = min(scoresMax(~keepMe_));
    keepMe = find(keepMe_);
    mapTo(keepMe_) = keepMe; % map units to keep to themselves

    keepMe = setdiff(keepMe, mapTo(~keepMe_));
    removeMe = setdiff(1:nClustersPrev, mapTo);
    updateMe = setdiff(setdiff(1:nClustersPrev, keepMe), removeMe);

    spikeClusters_ = obj.spikeClusters;
    good = spikeClusters_ > 0;
    spikeClusters_(good) = int32(mapTo(spikeClusters_(good))); % translate cluster number

    obj.subsetFields(union(keepMe, updateMe));
    [~, ~, obj.spikeClusters(good)] = unique(spikeClusters_(good)); % remap to 1:nNewClusters

    updateMe = arrayfun(@(i) find(union(keepMe, updateMe) == i), updateMe);
    obj.refresh(0, updateMe); % recount

    arrayfun(@obj.rmRefracSpikes, updateMe); % remove refrac spikes
    obj.removeEmptyClusters();

    % update cluster waveforms and distance
    obj.computeMeanWaveforms(updateMe);
    obj.computeWaveformSim(updateMe);

    nMerged = nClustersPrev - obj.nClusters;
    obj.hCfg.updateLog('mergeBySim', sprintf('Finished merging clusters (was %d, now %d: %d merged; minimum score: %0.3f)', ...
                                             nClustersPrev, obj.nClusters, nMerged, minScore), 0, 1);
end