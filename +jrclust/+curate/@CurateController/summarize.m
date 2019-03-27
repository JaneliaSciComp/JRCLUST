function summarize(obj)
    %SUMMARIZE Summarize a JRCLUST session's results
    summaryText = cell(5, 1);

    nSpikes = numel(obj.hClust.spikeTimes);
    countsPerSite = cellfun(@numel, obj.hClust.spikesBySite);
    [minCount, argMin] = min(countsPerSite);
    [maxCount, argMax] = max(countsPerSite);
    medCount = median(countsPerSite);

    summaryText{1} = sprintf('Spike count: %d', nSpikes);
    summaryText{2} = sprintf('Spike counts per site: min %d (site %d), max %d (site %d), median %i', ...
                             minCount, argMin, maxCount, argMax, medCount);

    nClustersInitial = max(obj.hClust.initialClustering);
    nClusters = numel(obj.hClust.unitCount);

    if jrclust.utils.isEqual(obj.hClust.initialClustering, obj.hClust.spikeClusters)
        summaryText{3} = sprintf('Clusters: %d (no merges)', nClusters);
    else
        summaryText{3} = sprintf('Clusters: %d (initial), %d (current)', nClustersInitial, nClusters);
    end

    % print spike counts
    [minCount, argMin] = min(obj.hClust.unitCount);
    [maxCount, argMax] = max(obj.hClust.unitCount);
    medCount = median(obj.hClust.unitCount);
    summaryText{4} = sprintf('Spike count per cluster: min %d (cluster %d), max %d (cluster %d), median %i', ...
                             minCount, argMin, maxCount, argMax, medCount);

    % print site counts
    nSitesByCluster = cellfun(@(x) numel(unique(obj.hClust.spikeSites(x))), obj.hClust.spikesByCluster);
    [minCount, argMin] = min(nSitesByCluster);
    [maxCount, argMax] = max(nSitesByCluster);
    medCount = median(nSitesByCluster);
    summaryText{5} = sprintf('Site count per cluster: min %d (cluster %d), max %d (cluster %d), median %i', ...
                             minCount, argMin, maxCount, argMax, medCount);

    jrclust.utils.qMsgBox(strjoin(summaryText, '\n'));
end

