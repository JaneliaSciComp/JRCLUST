function success = deleteMany(obj, unitIds)
%DELETEMANY Delete multiple units in one go.
%   Shift values in the spike table down, restructure metadata, add entries
%   to the history log.
success = 1;

% run through each unit in descending order to avoid reindexing
unitIds = sort(unitIds, 'descend');

for i = 1:numel(unitIds)
    success = obj.deleteSingle(unitIds(i));
    if ~success
        break;
    end
end
end
