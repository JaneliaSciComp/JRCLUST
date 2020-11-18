function convertHistory(obj)
%CONVERTHISTORY Convert history from histFile construct to struct.
if ~isfield(obj.res, 'hClust')
    return;
end

histFile = obj.hCfg.histFile;

if exist(histFile, 'file') ~= 2
    return;
end

% count the number of edits we've made
nSpikes = numel(obj.hClust.spikeTimes);

d = dir(histFile);
nEdits = d.bytes / (nSpikes + 1) / 4;

if nEdits == 0 % file is empty
    return;
end

%% record edits in history
history = struct('optype', cell(1), 'message', cell(1), 'indices', cell(1));

fid = fopen(histFile, 'r');
fread(fid, 1, 'int32'); % discard checkInt

beforeTable = fread(fid, nSpikes, 'int32');
for r = 2:nEdits
    fread(fid, 1, 'int32'); % discard checkInt
    afterTable = fread(fid, nSpikes, 'int32');

    if checkIsDelete(beforeTable, afterTable) % delete?
        [beforeIds, afterIds] = getDeletionMap(beforeTable, afterTable);
        for j = 1:numel(beforeIds)
            bid = beforeIds(j);
            aid = afterIds(j);
            history.optype{end+1} = 'delete';
            history.indices{end+1} = [bid, aid];
            history.message{end+1} = sprintf('deleted %d', bid);
        end
    elseif checkIsUndelete(beforeTable, afterTable) % undelete?
        [beforeIds, afterIds] = getUndeletionMap(beforeTable, afterTable);
        for j = 1:numel(beforeIds)
            bid = beforeIds(j);
            aid = afterIds(j);
            history.optype{end+1} = 'undelete';
            history.indices{end+1} = [bid, aid];
            history.message{end+1} = sprintf('undeleted %d', bid);
        end
    elseif checkIsMerge(beforeTable, afterTable) % merge?
        [unitIds, partitioning] = getMergeMap(beforeTable, afterTable);
        history.optype{end+1} = 'merge';
        history.indices{end+1} = {unitIds; partitioning};
        history.message{end+1} = sprintf('merged %s -> %d', ...
            jrclust.utils.field2str(unitIds), unitIds(1));
    elseif checkIsSplit(beforeTable, afterTable) % split?
        [unitIds, partitioning] = getSplitMap(beforeTable, afterTable);
        history.optype{end+1} = 'split';
        history.indices{end+1} = {unitIds; partitioning};
        history.message{end+1} = sprintf('splitted %d -> %s', unitIds(1), ...
            jrclust.utils.field2str(unitIds));
    else % something else
        maxBefore = max(beforeTable);
        beforeAfter = arrayfun(@(iC) unique(afterTable(beforeTable == iC)), ...
            1:maxBefore, 'UniformOutput', 0);

        if all(cellfun(@numel, beforeAfter) == 1) % a permutation
            history.optype{end+1} = 'reorder';
            history.indices{end+1} = [(1:maxBefore)' [beforeAfter{:}]'];

            % try to use the old commit message if we can
            if isa(obj.hClust.history, 'containers.Map') && obj.hClust.history.isKey(r)
                history.message{end+1} = obj.hClust.history(r);
            else
                history.message{end+1} = 'reorder units';
            end
        else
            history.optype{end+1} = 'reassign';
            history.indices{end+1} = [beforeTable(:) afterTable(:)];

            % try to use the old commit message if we can
            if isa(obj.hClust.history, 'containers.Map') && obj.hClust.history.isKey(r)
                history.message{end+1} = obj.hClust.history(r);
            else
                history.message{end+1} = 'reassign spikes';
            end
        end
    end

    beforeTable = afterTable;
end

fclose(fid);

% save res
obj.hClust.history = history;
try
    obj.saveRes(1);
catch ME
    warning('Unable to save history file: ''%s''. Bailing out.', ME.message);
    return;
end

% clean up
if obj.hCfg.getOr('testRun', 0)
    deleteFile(obj.hCfg.histFile); % don't prompt
else
    deleteFile(obj.hCfg.histFile, 'History has been converted. Delete the old history file?');
end

end %fun

%% LOCAL FUNCTIONS
function isDelete = checkIsDelete(beforeTable, afterTable)
%CHECKISDELETE Given two spike tables, a delete of one or more units
%implies the following:
%1. the maximum value of the after table is strictly less than that of the
%   before table, as all units must be shifted down.
%2. the minimum value of the after table is strictly less than that of the
%   before table, as all newly-deleted units must be negative.
%3. each new negative unit in the after table maps to exactly one positive
%   unit in the before table.
isDelete = 0;

if max(afterTable) >= max(beforeTable)
    return;
end

minAfter = min(afterTable);
minBefore = min(beforeTable);
if minAfter >= minBefore || minAfter >= 0
    return;
end

afterIds = minBefore-1:-1:minAfter;
beforeIds = zeros('like', afterIds);

% check that each deleted unit maps to exactly one prior, positive, unit
for i = numel(afterIds)
    uid = afterIds(i);
    mask = (afterTable == uid);

    % beforeEntry should be unique and positive
    beforeEntry = unique(beforeTable(mask));
    if numel(beforeEntry) > 1 || beforeEntry <= 0 || ismember(beforeEntry, beforeIds)
        return;
    end

    beforeIds(i) = beforeEntry;
end

isDelete = 1;
end %fun

function [beforeIds, afterIds] = getDeletionMap(beforeTable, afterTable)
%GETDELETIONMAP Assuming the mapping between beforeTable and afterTable is
%a delete, return the unique ids being deleted along with their (negative)
%ids as deleted units.
minAfter = min(afterTable);
minBefore = min(-1, min(beforeTable) - 1);

afterIds = (minBefore:-1:minAfter)';
beforeIds = zeros('like', afterIds);

for i = 1:numel(afterIds)
    uid = afterIds(i);
    mask = (afterTable == uid);
    beforeIds(i) = unique(beforeTable(mask)); % the unique assumption is important here!
end

% we need to perform the deletes in a descending fashion
[beforeIds, argsort] = sort(beforeIds, 'desc');
afterIds = afterIds(argsort);
end %fun

function isUndelete = checkIsUndelete(beforeTable, afterTable)
%CHECKISUNDELETE Given two spike tables, an undelete of one or more units
%implies the following:
%1. the maximum value of the after table is strictly greater than that of
%   the before table, as all units must be shifted up.
%2. the number of unique negative entries in the after table must be
%   strictly less than that of the before table, as all undeleted units
%   must be taken from negative values, possibly leaving gaps.
%3. each missing negative unit in the before table maps to exactly one
%   positive unit in the after table.
isUndelete = 0;

beforeUnits = unique(beforeTable);
afterUnits = unique(afterTable);

if max(afterUnits) <= max(beforeUnits)
    return;
end

negativeBefore = beforeUnits(beforeUnits < 0);
negativeAfter = afterUnits(afterUnits < 0);
if numel(negativeAfter) >= numel(negativeBefore)
    return;
end

missingNegatives = setdiff(negativeBefore, negativeAfter);
for i = 1:numel(missingNegatives)
    beforeEntry = missingNegatives(i);
    mask = (beforeTable == beforeEntry);

    afterEntry = unique(afterTable(mask));
    if numel(afterEntry) > 1 || afterEntry <= 0
        return;
    end
end

isUndelete = 1;
end %fun

function [beforeIds, afterIds] = getUndeletionMap(beforeTable, afterTable)
%GETUNDELETIONMAP Assuming the mapping between beforeTable and afterTable
%is an undelete, return the unique (negative) ids being undeleted along
%with their (positive) undeleted ids.
beforeUnits = unique(beforeTable);
negativeBefore = beforeUnits(beforeUnits < 0);

afterUnits = unique(afterTable);
negativeAfter = afterUnits(afterUnits < 0);

missingNegatives = setdiff(negativeBefore, negativeAfter);
beforeIds = zeros(numel(missingNegatives), 1);
afterIds = zeros(numel(missingNegatives), 1);

for i = 1:numel(missingNegatives)
    uid = missingNegatives(i);
    beforeIds(i) = uid;

    mask = (beforeTable == uid);    
    afterIds(i) = unique(afterTable(mask)); % the unique assumption is important here!
end
end %fun

function isMerge = checkIsMerge(beforeTable, afterTable)
%CHECKISMERGE Given two spike tables, a merge of two or more units implies
%the following:
%1. the maximum value of the after table is strictly less than that of the
%   before table, as all units must be shifted down.
%2. the minimum value of the after table is exactly equal to the minimum
%   value of the before table, as no units have been deleted.
%3. exactly one unit in the after table must map to more than one unit in
%   the before table, and that unit in the after table must be the minimum
%   of the corresponding units in the before table.
isMerge = 0;

afterIds = unique(afterTable);
if max(afterIds) >= max(beforeTable) || min(afterIds) ~= min(beforeTable)
    return;
end

afterIds = afterIds(afterIds > 0);
candidate = [];

for i = 1:numel(afterIds)
    aid = afterIds(i);

    mask = (afterTable == aid);
    beforeEntries = unique(beforeTable(mask));
    if numel(beforeEntries) > 1 && aid == min(beforeEntries) && isempty(candidate)
        candidate = aid;
    elseif numel(beforeEntries) > 1
        return;
    end
end

isMerge = 1;
end %fun

function [unitIds, partitioning] = getMergeMap(beforeTable, afterTable)
%GETMERGEMAP Assuming the mapping between beforeTable and afterTable
%is a merge, return the unit ids in the before table and their indices as a
%cell array representing the partitioning.
afterIds = unique(afterTable);
afterIds = afterIds(afterIds > 0);

for i = 1:numel(afterIds)
    aid = afterIds(i);

    mask = (afterTable == aid);
    unitIds = unique(beforeTable(mask));
    if numel(unitIds) > 1 && aid == min(unitIds)
        break;
    end
end

unitIds = sort(unitIds);
partitioning = arrayfun(@(iC) find(beforeTable == iC), unitIds, 'UniformOutput', 0);
end %fun

function isSplit = checkIsSplit(beforeTable, afterTable)
%CHECKISSPLIT Given two spike tables, a split of one unit into two or more
%units implies the following:
%1. the maximum value of the after table is strictly greater than that of
%   the before table, as all units must be shifted up.
%2. the minimum value of the after table is exactly equal to the minimum
%   value of the before table, as no units have been deleted.
%3. exactly one unit in the before table must map to more than one unit in
%   the after table, and that unit in the before table must be the minimum
%   of the corresponding units in the after table.
isSplit = 0;

beforeIds = unique(beforeTable);
if max(beforeIds) >= max(afterTable) || min(beforeIds) ~= min(afterTable)
    return;
end

beforeIds = beforeIds(beforeIds > 0);
candidate = [];

for i = 1:numel(beforeIds)
    bid = beforeIds(i);

    mask = (beforeTable == bid);
    afterEntries = unique(afterTable(mask));
    if numel(afterEntries) > 1 && bid == min(afterEntries) && isempty(candidate)
        candidate = bid;
    elseif numel(afterEntries) > 1
        return;
    end
end

isSplit = 1;
end %fun

function [unitIds, partitioning] = getSplitMap(beforeTable, afterTable)
%GETSPLITMAP Assuming the mapping between beforeTable and afterTable
%is a split, return the unit ids in the after table and their indices as a
%cell array representing the partitioning.
beforeIds = unique(beforeTable);
beforeIds = beforeIds(beforeIds > 0);

for i = 1:numel(beforeIds)
    bid = beforeIds(i);

    mask = (beforeTable == bid);
    unitIds = unique(afterTable(mask));
    if numel(unitIds) > 1 && bid == min(unitIds)
        break;
    end
end

unitIds = sort(unitIds);
partitioning = arrayfun(@(iC) find(afterTable == iC), unitIds, 'UniformOutput', 0);
end %fun