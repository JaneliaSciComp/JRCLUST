function success = mergeMultiple(obj, unitIds)
%MERGEMULTIPLE Merge two or more units.
%   Shift values in the spike table down, restructure metadata, add an entry
%   to the history log.
success = 0; %#ok<NASGU>

% take just the unique values and ensure they're sorted
unitIds = sort(unique(unitIds));

% unitIds does not contain at least 2 values, error
if numel(unitIds) < 2
    error('call to mergeMultiple with %d unique units (requires 2 or more)', numel(unitIds));
end

% noise or negative units, error
if any(unitIds < 1)
    error('call to mergeMultiple includes noise or deleted units');
end

% some of the given units not found in the spike table, error
mergedMask = ismember(obj.spikeClusters, unitIds);
if ~all(ismember(unitIds, obj.spikeClusters))
    error('call to mergeMultiple includes nonexistent units');
end

res = struct(); % speculative merge, ensure consistency before committing
res.spikeClusters = obj.spikeClusters;

% indices of subset for metadata
subset = setdiff(1:obj.nClusters, unitIds(2:end));
nUnitsAfter = numel(subset);

%% update spike table

% keep a record of the indices used for merging so we can split later if
% need be
partitioning = arrayfun(@(iC) find(res.spikeClusters == iC), unitIds, 'UniformOutput', 0);
res.spikeClusters(mergedMask) = unitIds(1);

% units will have to shift down by this amount
shiftBy = arrayfun(@(i) sum(unitIds(2:end) <= i), 1:obj.nClusters);

for i = 1:obj.nClusters
    if shiftBy(i) == 0
        continue;
    end
    mask = (res.spikeClusters == i);
    res.spikeClusters(mask) = res.spikeClusters(mask) - shiftBy(i);
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
        obj.history.optype{end+1} = 'merge';
        obj.history.message{end+1} = sprintf('merged %s -> %d', ...
            jrclust.utils.field2str(unitIds), unitIds(1));
        obj.history.indices{end+1} = {unitIds; partitioning}; % before, after

        % update units that need to be recomputed

        % shift every unit above any merging units down by the correct
        % amount
        shiftBy = arrayfun(@(i) sum(unitIds(2:end) <= i), obj.recompute);
        obj.recompute = unique(obj.recompute - shiftBy);

        % we'll need to recompute the merged unit
        if ~ismember(unitIds(1), obj.recompute)
            obj.recompute(end+1) = unitIds(1);
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
end

