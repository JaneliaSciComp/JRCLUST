function argsort = orderClusters(obj, by)
    %ORDERCLUSTERS Arrange cluster ID numbers by some criterion
    if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
        warning('cannot branch from history; use revert() first');
        return;
    end

    if nargin < 2 || isempty(by) || (~strcmp(by, 'Y + X') && ~isprop(obj, by))
        by = 'clusterSites';
    end

    if strcmpi(by, 'Y + X') && ~isempty(obj.clusterCentroids)
        [~, argsort] = sort(sum(obj.clusterCentroids, 2), 'ascend');
    elseif isprop(obj, by)
        [~, argsort] = sort(obj.(by), 'ascend');
    end

    if issorted(argsort)
        return;
    end

    %obj.spikeClusters = mapIndex_(obj.spikeClusters, argsort);
    map(argsort) = 1:numel(argsort);
    mask = (obj.spikeClusters > 0);
    obj.spikeClusters(mask) = map(obj.spikeClusters(mask)); % do not map zeros

    % reorder data fields
    obj.subsetFields(argsort);
%     obj.spikesByCluster = obj.spikesByCluster(argsort);
%     obj.clusterSites = obj.clusterSites(argsort);
%     obj.unitCount = obj.unitCount(argsort);
%     obj.clusterNotes = obj.clusterNotes(argsort);
%     if ~isempty(obj.clusterCentroids)
%         obj.clusterCentroids = obj.clusterCentroids(argsort, :);
%     end
end