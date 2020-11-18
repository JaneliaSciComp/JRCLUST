function success = reorderMultiple(obj, beforeIds, afterIds)
%REORDERMULTIPLE Reorder units.
%   Rearrange values in the spike table and reorder metadata fields.
success = 0; %#ok<NASGU>

% before and after need to match up
if numel(beforeIds) ~= numel(afterIds)
    error('call to reorderMultiple with mismatched unit counts');
end

% each unit must be in both beforeIds and afterIds
if ~isempty(setdiff(beforeIds, afterIds)) || ~isempty(setdiff(afterIds, beforeIds))
    error('some units not accounted for in both sets')
end

% some ids are nonpositive or not found in the spike table
if any(beforeIds <= 0) || any(~ismember(beforeIds, obj.spikeClusters))
    error('call to reorderMultiple with noise or deleted units, or non-units');
end

res = struct(); % speculative reorder, ensure consistency before committing
res.spikeClusters = obj.spikeClusters;

subset = 1:obj.nClusters;
subset(beforeIds) = afterIds;
nUnitsAfter = numel(subset);

%% update spike table
for i = 1:numel(beforeIds)
    res.spikeClusters(obj.spikeClusters == beforeIds(i)) = afterIds(i);
end

isConsistent = 1;

%% subset metadata

% 1-D fields
fieldnames_ = obj.unitFields.vectorFields; % names of 1D metadata fields
for i = 1:numel(fieldnames_)
    fn = fieldnames_{i};
    val = obj.(fn);

    % no values to subset, nothing to do
    if isempty(val)
        continue;
    end

    % reorder entries
    val = val(subset);

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

% 2- and 3-D fields with special subsetting rules
% (see /json/Clustering.json for subsetting functions for non-vector
% fields)
otherFields = obj.unitFields.otherFields;
fieldnames_ = fieldnames(otherFields);
for i = 1:numel(fieldnames_)
    fn = fieldnames_{i};
    val = obj.(fn);

    % no values to subset, nothing to do
    if isempty(val)
        continue;
    end

    % evaluate subset and consistency check functions
    hFunSubs = eval(otherFields.(fn).subset);
    hFunConsistent = eval(otherFields.(fn).consistent);

    % reorder entries
    val = hFunSubs(val, subset);

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
        obj.history.optype{end+1} = 'reorder';
        obj.history.message{end+1} = 'reorder units';
        obj.history.indices{end+1} = [beforeIds(:) afterIds(:)]; % [before, after]

        % update units that need to be recomputed
        for i = 1:numel(beforeIds)
            obj.recompute(obj.recompute == beforeIds(i)) = afterIds(i);
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
end % fun
