function success = commit(obj, spikeClusters, metadata, msg)
%COMMIT Commit a modification of clustering to history log
success = 0;

if ~isempty(obj.recompute)
    obj.spikesByCluster(flagged) = arrayfun(@(iC) find(obj.spikeClusters == iC), flagged, 'UniformOutput', 0);
    % update count of spikes per unit
    obj.unitCount(flagged) = cellfun(@numel, obj.spikesByCluster(flagged));

    % update cluster sites
    obj.clusterSites(flagged) = double(arrayfun(@(iC) mode(obj.spikeSites(obj.spikesByCluster{iC})), flagged));

    % update mean waveforms for flagged units
    obj.updateWaveforms(flagged);

    % update unit positions
    obj.computeCentroids(flagged);

    % compute quality scores for altered units
    obj.computeQualityScores(flagged);
end

success = 1;
end