function nMerged = mergeBySim(obj)
    %MERGEBYSIM Automatically merge clusters by similarity score
    obj.hCfg.updateLog('mergeBySim', 'Merging clusters by similarity', 1, 0);
    waveformSim = obj.waveformSim;
    nClustersPrev = size(waveformSim, 2);

    % ensure targets of merges are smaller than units merging into them
    waveformSim(tril(true(nClustersPrev))) = 0;
    [maxScores, mergeTargets] = max(waveformSim);

    % keep clusters whose maximum similarity to some other cluster
    % is LESS than our threshold
    underThresh = maxScores < obj.hCfg.maxUnitSim;

    if all(underThresh)
        nMerged = 0;
        obj.hCfg.updateLog('mergeBySim', 'No clusters to merge', 0, 1);
        return;
    end

    mergingUnits = find(~underThresh);

    mergeTargets = mergeTargets(mergingUnits);
    % intersection represents merge targets which are also set to merge
    % into some other unit, and thus shouldn't be merged down (but can
    % still be merged into)
    intersection = ismember(mergingUnits, mergeTargets);
    mergingUnits(intersection) = [];
    mergeTargets(intersection) = [];

    nMerged = numel(mergingUnits);

    res = obj.mergeUnits(obj.spikeClusters, mergeTargets, mergingUnits);
    msg = sprintf('merge %d units by similarity', nMerged);
    obj.commit(res.spikeClusters, res.metadata, msg);
    
    obj.hCfg.updateLog('mergeBySim', sprintf('Finished merging clusters (was %d, now %d: %d merged; minimum score: %0.3f)', ...
                                             nClustersPrev, obj.nClusters, nMerged, min(maxScores(~underThresh))), 0, 1);
end