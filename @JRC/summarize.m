function summarize(obj)
    %SUMMARIZE Summarize a JRCLUST session's results
    if ~isfield(obj.res, 'spikeTimes')
        return;
    end

    if obj.isDetect
        nSpikes = numel(obj.res.spikeTimes);
        countsPerSite = cellfun(@numel, obj.res.spikesBySite);
        [minCount, argMin] = min(countsPerSite);
        [maxCount, argMax] = max(countsPerSite);
        medCount = median(countsPerSite);

        fprintf('\n====DETECTION SUMMARY====\n');
        fprintf('Detection completed in %0.2f s\n', obj.res.detectTime);
        fprintf('Spike count: %d\n', nSpikes);
        fprintf('Spike counts per site: min %d (site %d), max %d (site %d), median %i\n', ...
                minCount, argMin, maxCount, argMax, medCount);
    end

    if obj.isSort && isfield(obj.res, 'hClust')
        nClustersInitial = max(obj.res.hClust.initialClustering);
        nClusters = numel(obj.res.hClust.unitCount);

        fprintf('\n====SORTING SUMMARY====\n');
        fprintf('Sorting completed in %0.2f s\n', obj.res.sortTime);
        if nClustersInitial == nClusters
            fprintf('Clusters: %d (no merges)\n', nClusters);
        else
            fprintf('Clusters: %d (initial) merged to %d (final)\n', nClustersInitial, nClusters);
        end

        % print spike counts
        [minCount, argMin] = min(obj.res.hClust.unitCount);
        [maxCount, argMax] = max(obj.res.hClust.unitCount);
        medCount = median(obj.res.hClust.unitCount);
        fprintf('Spike count per cluster: min %d (cluster %d), max %d (cluster %d), median %i\n', ...
                minCount, argMin, maxCount, argMax, medCount);

        % print site counts
        nSitesByCluster = cellfun(@(x) numel(unique(obj.res.spikeSites(x))), obj.res.hClust.spikesByCluster);
        [minCount, argMin] = min(nSitesByCluster);
        [maxCount, argMax] = max(nSitesByCluster);
        medCount = median(nSitesByCluster);
        fprintf('Site count per cluster: min %d (cluster %d), max %d (cluster %d), median %i\n', ...
                minCount, argMin, maxCount, argMax, medCount);
    end
end

