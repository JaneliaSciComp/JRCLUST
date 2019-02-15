function removeEmptyClusters(obj)
    %REMOVEEMPTYCLUSTERS Find and remove empty clusters
    keepClusters = obj.unitCount > 0;
    if all(keepClusters)
        return;
    end

    % subset all fields indexed by cluster
    obj.subsetFields(keepClusters);

    if min(obj.spikeClusters) < 1 % noise or garbage clusters
        obj.spikeClusters(obj.spikeClusters < 1) = 0;

        % renumber clusters 1:numel(unique(spikeClusters))
        [~, ~, obj.spikeClusters] = unique(obj.spikeClusters + 1);
        obj.spikeClusters = obj.spikeClusters - 1;
    else
        % renumber clusters 1:numel(unique(spikeClusters))
        [~, ~, obj.spikeClusters] = unique(obj.spikeClusters);
    end

    obj.spikeClusters = int32(obj.spikeClusters);
end