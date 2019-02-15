function refresh(obj, doRemoveEmpty, updateMe)
    %REFRESH Recount and store spikes by cluster, optionally removing empty clusters
    if nargin < 2
        doRemoveEmpty = 0;
    end
    if nargin < 3
        updateMe = [];
    end

    if isempty(updateMe)
        obj.spikesByCluster = arrayfun(@(iC) find(obj.spikeClusters == iC), 1:obj.nClusters, 'UniformOutput', 0);
        obj.unitCount = cellfun(@numel, obj.spikesByCluster);
        if ~isempty(obj.spikeSites)
            obj.clusterSites = double(arrayfun(@(iC) mode(obj.spikeSites(obj.spikesByCluster{iC})), 1:obj.nClusters));
        else
            obj.clusterSites = [];
        end
    else
        obj.spikesByCluster(updateMe) = arrayfun(@(iC) find(obj.spikeClusters == iC), updateMe, 'UniformOutput', 0);
        obj.unitCount(updateMe) = cellfun(@numel, obj.spikesByCluster(updateMe));
        if ~isempty(obj.spikeSites)
            obj.clusterSites(updateMe) = double(arrayfun(@(iC) mode(obj.spikeSites(obj.spikesByCluster{iC})), updateMe));
        else
            obj.clusterSites = [];
        end
    end

    if doRemoveEmpty
        obj.removeEmptyClusters();
    end
end