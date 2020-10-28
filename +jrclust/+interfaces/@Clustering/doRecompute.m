function success = doRecompute(obj)
%DORECOMPUTE Recompute derivative properties of the spike table after an
%edit.
%   Mean waveforms, centroids, quality scores, and aggregated indices are
%   all computed here.
success = 1;

recompute = obj.recompute(:);
if numel(obj.recompute) == 0
    recompute = (1:obj.nClusters)';

    obj.clusterNotes = cell(obj.nClusters, 1);
end

% take a census of each cluster
try
    spikesByCluster = arrayfun(@(iC) find(obj.spikeClusters == iC), recompute, 'UniformOutput', 0); 
    clusterSites = double(cellfun(@(sbc) mode(obj.spikeSites(sbc)), spikesByCluster));
catch ME
    warning(ME.message);
    success = 0;
    return;
end

obj.spikesByCluster(recompute) = spikesByCluster;
obj.clusterSites(recompute) = clusterSites;

% compute quality scores
try
    obj.computeQualityScores(recompute);
catch ME
    warning(ME.message);
    success = 0;
end

% reset notes only from affected clusters
obj.clusterNotes(recompute) = arrayfun(@(i) '', recompute, 'UniformOutput', 0);
end % fun

