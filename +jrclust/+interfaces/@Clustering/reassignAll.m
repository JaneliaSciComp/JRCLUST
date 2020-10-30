function success = reassignAll(obj, beforeTable, afterTable)
%REASSIGNALL Reassign all spikes.
%   Force recompute all metadata fields.
success = 0; %#ok<NASGU>

% before and after need to match up
if numel(beforeTable) ~= numel(afterTable)
    error('call to reassignAll with mismatched unit counts');
end

res = struct(); % speculative reorder, ensure consistency before committing

%% update spike table
res.spikeClusters = afterTable;

isConsistent = 1;

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
        obj.history.optype{end+1} = 'reassign';
        obj.history.message{end+1} = 'reassign all spikes';
        obj.history.indices{end+1} = [beforeTable(:), afterTable(:)]; % [before, after]

        % update units that need to be recomputed
        recompute = unique(sort(afterTable));
        obj.recompute = recompute(recompute > 0);
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

end %fun