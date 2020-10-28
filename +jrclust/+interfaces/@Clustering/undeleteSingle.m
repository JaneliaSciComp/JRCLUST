function success = undeleteSingle(obj, deletedId, newId)
%UNDELETESINGLE Undelete a single deleted unit, given by `deletedId`, and
%insert it into the spike table at `newId`.
%   Shift values in the spike table up, restructure metadata, add an entry
%   to the history log.
success = 0; %#ok<NASGU>

% deletedId or newId is not a single value, error
if numel(deletedId) ~= 1 || numel(newId) ~= 1
    error('call to undeleteSingle with %d deleted units, %d new units', ...
          numel(deletedId), numel(newId));
end

% deletedId is nonnegative or not found in the spike table, error
deletedMask = ismember(obj.spikeClusters, deletedId);
if deletedId > -1 || ~any(deletedMask)
    error('call to undeleteSingle with an illegal deletedId %d', deletedId);
end

% newId is nonpositive, error
if newId < 1
    error('call to undeleteSingle with a nonpositive newId %d', newId);
end
% inserting newId into the spike table would create a gap, error
if newId - max(obj.spikeClusters) > 1
    error('call to undeleteSingle with a newId of %d (would create a gap of %d)', ...
          newId, newId - max(obj.spikeClusters));
end

res = struct(); % speculative undelete, ensure consistency before committing
res.spikeClusters = obj.spikeClusters;
shiftMask = res.spikeClusters >= newId; % all units to shift up
nUnitsAfter = obj.nClusters + 1;

%% update spike table

% shift larger units up by 1
res.spikeClusters(shiftMask) = res.spikeClusters(shiftMask) + 1;

% set indices to newId
res.spikeClusters(res.spikeClusters == deletedId) = newId;

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

    % shift all entries up by 1
    if isrow(val)
        val = [val(1:newId-1) makeFiller(val, 1) val(newId:end)];
    else
    	val = [val(1:newId-1); makeFiller(val, 1); val(newId:end)];
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
    val = hFunAug(val, newId-1, filler);

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
        obj.history.optype{end+1} = 'undelete';
        obj.history.message{end+1} = sprintf('undeleted %d', newId);
        obj.history.indices{end+1} = [deletedId, newId]; % [before, after]

        % update units that need to be recomputed

        % shift every unit above this up
        shiftMask = obj.recompute >= newId;
        obj.recompute(shiftMask) = obj.recompute(shiftMask) + 1;

        % we will need to recompute this unit as well
        if ~ismember(newId, obj.recompute)
            obj.recompute(end+1) = newId;
        end
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
end % func
