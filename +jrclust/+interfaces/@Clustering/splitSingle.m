function success = splitSingle(obj, unitIds, partitioning)
%SPLITSINGLE Split a unit according to a given partitioning.
%   The partitioning is expected in terms of *ABSOLUTE INDICES*.
success = 0; %#ok<NASGU>

% unitIds needs to be unique
if numel(unique(unitIds)) ~= numel(unitIds)
    error('call to splitSingle contains duplicate units');
end

% unitIds does not contain at least 2 values, error
if numel(unitIds) < 2
    error('call to splitSingle with %d unique units (requires 2 or more)', numel(unitIds));
end

% noise or negative units, error
if any(unitIds < 1)
    error('call to splitSingle includes noise or deleted units');
end

% each id in unitIds needs to correspond to one partition
if numel(unitIds) ~= numel(partitioning)
    error('call to splitSingle with an incorrect label-partition correspondence');
end

% ensure unitIds are sorted
[unitIds, argsort] = sort(unitIds);
partitioning = partitioning(argsort);

% splitting unit is not found in the spike table, error
splitMask = ismember(obj.spikeClusters, unitIds(1));
if ~any(splitMask)
    error('splitting unit not found');
end

% unpack all spike indices in the partitioning and check we have the full
% set, error if not
departitioning = [];
for i = 1:numel(unitIds)
    part = partitioning{i};
    departitioning = [departitioning; part(:)];
end
departitioning = sort(departitioning(:));
if ~jrclust.utils.isEqual(find(splitMask), departitioning) %#ok<FNDSB>
    error('partitioning does not account for all spikes in unit to be split');
end

res = struct(); % speculative split, ensure consistency before committing
res.spikeClusters = obj.spikeClusters;
nUnitsAfter = obj.nClusters + numel(unitIds) - 1;

%% update spike table

% shift up larger units
for i = 2:numel(unitIds)
    shiftMask = res.spikeClusters >= unitIds(i);
    res.spikeClusters(shiftMask) = res.spikeClusters(shiftMask) + 1;
    res.spikeClusters(partitioning{i}) = unitIds(i);
end

isConsistent = 1;

%% augment metadata

% 1-D fields
fieldnames_ = obj.unitFields.vectorFields; % names of 1D metadata fields
for i = 1:numel(fieldnames_)
    fn = fieldnames_{i};
    val = obj.(fn);

    % no values to augment, nothing to do
    if isempty(val)
        continue;
    end

    % shift up entries for larger units
    filler = makeFiller(val, 1);

    for i = 2:numel(unitIds)
        uid = unitIds(i);
        if isrow(val)
            val = [val(1:uid-1) filler val(uid:end)];
        else
            val = [val(1:uid-1); filler; val(uid:end)];
        end
    end

    % brief sanity check
    if numel(val) ~= nUnitsAfter
        isConsistent = 0;
        break;
    end

    res.(fn) = val;
end

% consistency check failed on some vector field
if ~isConsistent
    success = 0;
    return;
end
% 2- and 3-D fields with special augmenting rules
% (see /json/Clustering.json for augmenting functions for non-vector
% fields)
otherFields = obj.unitFields.otherFields;
fieldnames_ = fieldnames(otherFields);
for i = 1:numel(fieldnames_)
    fn = fieldnames_{i};
    val = obj.(fn);

    % no values to augment, nothing to do
    if isempty(val)
        continue;
    end

    switch fn
        case {'templateSim', 'waveformSim'}
            augShape = [obj.nClusters + 1, 1];

        case 'clusterCentroids'
            augShape = [1, 2];

        otherwise
            clusterDims = otherFields.(fn).cluster_dims;
            % get remaining (non-cluster-indexed) dimensions
            shape = size(val);
            augShape = [shape(1:clusterDims-1) 1 shape(clusterDims+1:end)];
    end

    % evaluate augment and consistency check functions
    hFunAug = eval(otherFields.(fn).augment);
    hFunConsistent = eval(otherFields.(fn).consistent);

    % remove entries corresponding to deleted units
    filler = makeFiller(val, augShape);
    for i = 2:numel(unitIds)
        val = hFunAug(val, unitIds(i)-1, filler);
    end

    % brief sanity check
    if ~hFunConsistent(val, nUnitsAfter)
        isConsistent = 0;
        break;
    end

    res.(fn) = val;
end


% commit spike table and fields
backup = struct();

if isConsistent
    fieldnames_ = fieldnames(res);
    for i = 1:numel(fieldnames_)
        fn = fieldnames_{i};
        backup.(fn) = obj.(fn);

        try
            obj.(fn) = res.(fn);
        catch ME
            warning('Failed to update field: %s', ME.message);
            isConsistent = 0;
            break;
        end
    end

    % success! commit to history log
    try
        obj.history.optype{end+1} = 'split';
        obj.history.message{end+1} = sprintf('splitted %d -> %s', ...
            unitIds(1), jrclust.utils.field2str(unitIds));
        obj.history.indices{end+1} = {unitIds; partitioning}; % before, after

        % update units that need to be recomputed

        % shift up larger units
        for i = 2:numel(unitIds)
            shiftMask = obj.recompute >= unitIds(i);
            obj.recompute(shiftMask) = obj.recompute(shiftMask) + 1;
        end

        % we will need to recompute all splitted units
        obj.recompute = union(obj.recompute, unitIds);
    catch ME
        warning('Failed to commit: %s', ME.message);
        isConsistent = 0;
    end
end

% restore backed up fields
if ~isConsistent
    fieldnames_ = fieldnames(backup);
    for i = 1:numel(fieldnames_)
        fn = fieldnames_{i};
        obj.(fn) = backup.(fn);
    end
end

success = isConsistent;
end
