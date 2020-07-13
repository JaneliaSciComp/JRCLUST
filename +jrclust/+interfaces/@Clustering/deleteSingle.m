function success = deleteSingle(obj, unitId)
%DELETESINGLE Delete a single unit.
%
% Shift values in the spike table down, restructure metadata.
success = 0; %#ok<NASGU>

% unitId not found in the spike table, nothing to do
unitMask = ismember(obj.spikeClusters, unitId);
if unitId < 1 || ~any(unitMask)
    success = 1; % vacuously true
    return;
end

res = struct(); % speculative delete, ensure consistency before committing

% do not delete in place!
res.spikeClusters = obj.spikeClusters;
shiftMask = res.spikeClusters > unitId; % all units to shift down

% indices of subset for metadata
subset = unique(res.spikeClusters((res.spikeClusters > 0) & (~unitMask)));
nUnitsAfter = numel(subset);

%% update spike table

% set indices to a negative value
deletedId = min(min(res.spikeClusters), 0) - 1; % save this for commit
res.spikeClusters(unitMask) = deletedId;

% shift larger units down by 1
res.spikeClusters(shiftMask) = res.spikeClusters(shiftMask) - 1;

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

    % remove entries corresponding to deleted units
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

    % remove entries corresponding to deleted units
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
        obj.history.message{end+1} = sprintf('deleted %d', unitId);
        obj.history.indices{end+1} = deletedId;

        % update units that need to be recomputed
        if ismember(unitId, obj.recompute)
            % don't need to recompute this unit as it's been deleted
            obj.recompute(obj.recompute == unitId) = [];
        end

        % shift every unit above this down
        shiftMask = obj.recompute > unitId;
        obj.recompute(shiftMask) = obj.recompute(shiftMask) - 1;
    catch ME
        warning('Failed to commit: %s', ME.message);
        isConsistent = 0;
    end
end

% restore backed up fields
if ~isConsistent
    fieldnames_ = fieldnames(backup);
    for i = 1:numel(fieldnames_)
        fn = fieldnames_(i);
        obj.(fn) = backup.(fn);
    end
end

success = isConsistent;
end % func

