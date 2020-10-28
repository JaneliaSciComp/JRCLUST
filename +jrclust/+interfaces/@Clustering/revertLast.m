function success = revertLast(obj, n)
%REVERTLAST Revert the last `n` curation operations. 
%   Undoes curation operations.

% revert last 1 by default
if nargin < 1
    n = 1;
end

nops = obj.nEdits; % this changes size with each iteration, so fix it here

if n < 1 || n > nops
    error('call to revertLast with illegal argument %d', n);
end

%% perform the revert

for i = nops:-1:nops-n+1
    optype = obj.history.optype{i};
    indices = obj.history.indices{i};

    switch optype
        case 'delete'
            deletedId = indices(2); newId = indices(1);
            success = obj.undeleteSingle(deletedId, newId);
        case 'undelete'
            deletedId = indices(1); unitId = indices(2);
            success = obj.deleteSingle(unitId, deletedId);
        case 'merge'
            unitIds = indices{1}; partitioning = indices{2};
            success = obj.splitSingle(unitIds, partitioning);
        case 'split'
            unitIds = indices{1};
            success = obj.mergeMultiple(unitIds);
        case 'reorder'
            beforeIds = indices(:, 1);
            afterIds = indices(:, 2);
            success = obj.reorderMultiple(afterIds, beforeIds);
        case 'reassign'
            beforeTable = indices(:, 1);
            afterTable = indices(:, 2);
            success = obj.reassignAll(afterTable, beforeTable);
        otherwise
            success = 0;
    end
end

%% recompute all fields and truncate history entries
obj.doRecompute();

obj.history.optype(end - 2*n + 1:end) = [];
obj.history.message(end - 2*n + 1:end) = [];
obj.history.indices(end - 2*n + 1:end) = [];

end %fun

