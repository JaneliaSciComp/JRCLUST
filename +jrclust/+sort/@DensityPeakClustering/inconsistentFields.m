function inconsistencies = inconsistentFields(obj)
    %INCONSISTENTFIELDS Check all fields have the correct sizes
    inconsistencies = inconsistentFields@jrclust.interfaces.Clustering(obj);

    % also check all clusterCenters are unique
    if ~ismember('clusterCenters', inconsistencies)
        % find all putative centers which don't actually belong to their
        % assigned cluster and reassign them
        dupes = find(obj.spikeClusters(obj.clusterCenters) ~= (1:obj.nClusters)');
        newCenters = zeros(numel(dupes), 1);
        for i = 1:numel(dupes)
            iDupe = dupes(i);
            iSpikes = find(obj.spikeClusters == iDupe);
            [~, iBest] = max(obj.spikeRho(iSpikes));
            newCenters(i) = iSpikes(iBest);
        end

        obj.clusterCenters(dupes) = newCenters;
    end
end
