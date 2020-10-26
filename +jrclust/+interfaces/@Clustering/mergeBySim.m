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

[mergingUnits, argsort] = sort(mergingUnits, 'descend');
mergeTargets = mergeTargets(argsort);
for i=1:numel(mergingUnits)
    unit = max(mergingUnits(i), mergeTargets(i));
    target = min(mergingUnits(i), mergeTargets(i));

    obj.mergeMultiple([target, unit]);

    mask = (1:numel(mergingUnits) > i);

    % adjust larger units to account for shift down
    mergingUnitsMask = mask & mergingUnits > unit;
    mergingUnits(mergingUnitsMask) = mergingUnits(mergingUnitsMask) - 1;

    mergeTargetsMask = mask & mergeTargets > unit;
    mergeTargets(mergeTargetsMask) = mergeTargets(mergeTargetsMask) - 1;
end

nMerged = numel(mergingUnits);

obj.hCfg.updateLog('mergeBySim', sprintf('Finished merging clusters (was %d, now %d: %d merged; minimum score: %0.3f)', ...
                                         nClustersPrev, obj.nClusters, nMerged, min(maxScores(~underThresh))), 0, 1);
end % fun