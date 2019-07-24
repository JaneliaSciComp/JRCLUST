function argsort = orderClusters(obj, by)
    %ORDERCLUSTERS Arrange cluster ID numbers by some criterion
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

    metadata = struct();

    map(argsort) = 1:numel(argsort);
    spikeClusters = obj.spikeClusters;
    mask = spikeClusters > 0; % don't remap noise or deleted units
    spikeClusters(mask) = map(obj.spikeClusters(mask));

    % subset vector fields
    vectorFields = obj.unitFields.vectorFields;
    hFunSubset = @(vals, indices) vals(indices);

    for i = 1:numel(vectorFields)
        fn = vectorFields{i};
        if ~isfield(metadata, fn)
            metadata.(fn) = obj.(fn);
        end

        if ~isempty(metadata.(fn))
            metadata.(fn) = hFunSubset(metadata.(fn), argsort);
        end
    end

    % subset other (n > 1) fields
    otherFields = obj.unitFields.otherFields;
    otherFieldNames = fieldnames(otherFields);

    for i = 1:numel(otherFieldNames)
        fn = otherFieldNames{i};
        if ~isfield(metadata, fn)
            metadata.(fn) = obj.(fn);
        end

        if ~isempty(metadata.(fn))
            hFunSubset = eval(otherFields.(fn).subset);
            metadata.(fn) = hFunSubset(metadata.(fn), argsort);
        end
    end

    obj.commit(spikeClusters, metadata, sprintf('reorder clusters by %s', by));
end