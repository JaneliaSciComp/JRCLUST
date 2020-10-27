function history = convertHistory(obj)
%CONVERTHISTORY Convert history from histFile construct to struct.
history = struct('optype', cell(1), 'message', cell(1), 'indices', cell(1));

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
fid = fopen(histFile, 'r');
fread(fid, 1, 'int32'); % discard checkInt

scBefore = fread(fid, nSpikes, 'int32');
for r = 2:nEdits
    fread(fid, 1, 'int32'); % discard checkInt
    scAfter = fread(fid, nSpikes, 'int32');
    diffMask = (scBefore ~= scAfter);
        
    before = scBefore(diffMask);
    after = scAfter(diffMask);

    deleteMask = (before > 0 & after < 0);
    if any(deleteMask) % a delete
        deletedUnits = unique(before(deleteMask));

        uid = deletedUnits(1);
        unitMask = deleteMask & (before == uid);
        deletedUid = unique(after(unitMask)); % assume 1 unique unit here

        history.optype{end+1} = 'delete';
        history.indices{end+1} = [uid, deletedUid];
        history.message{end+1} = sprintf('deleted %d', uid);

        shiftDownMask = (before > uid);
        before(shiftDownMask) = before(shiftDownMask) - 1;
        before(unitMask) = after(unitMask);
    end

    undeleteMask = (before < 0 & after > 0);
    if any(undeleteMask) % an undelete
        deletedUnits = unique(before(undeleteMask));

        deletedUid = deletedUnits(1);
        unitMask = undeleteMask & (before == deletedUid);
        uid = unique(after(unitMask)); % assume 1 unique unit here

        history.optype{end+1} = 'undelete';
        history.indices{end+1} = [deletedUid, uid];
        history.message{end+1} = sprintf('undeleted %d', uid);

        shiftUpMask = (before >= uid);
        before(shiftUpMask) = before(shiftUpMask) + 1;
        before(unitMask) = after(unitMask);
    end

    uniqueBefore = unique(before);
    maxBefore = max(uniqueBefore);

    uniqueAfter = unique(after);
    maxAfter = max(uniqueAfter);

    % a merge
    if maxBefore > maxAfter
        uid = uniqueAfter(1);
        unitMask = (after == uid);

        beforeUnits = unique(before(unitMask)); % these units merged into uid
        if numel(beforeUnits) > 1
            mergingUnits = [uid; beforeUnits(:)];
            partitioning = arrayfun(@(k) find(scBefore == k), mergingUnits, 'UniformOutput', 0);

            history.optype{end+1} = 'merge';
            history.indices{end+1} = {mergingUnits; partitioning};
            history.message{end+1} = sprintf('merged %s -> %d', ...
                jrclust.utils.field2str(mergingUnits), uid);

            before(unitMask) = uid;
            beforeUnits = sort(beforeUnits, 'descend');
            for j = 1:numel(beforeUnits)
                mask = before > beforeUnits(j);
                before(mask) = before(mask) - 1;
            end
        end
    end

    uniqueBefore = unique(before);
    maxBefore = max(uniqueBefore);

    uniqueAfter = unique(after);
    maxAfter = max(uniqueAfter);

    % a split
    if maxAfter > maxBefore
        nSplitOff = maxAfter - maxBefore;
        splitOffUnits = sort(uniqueAfter((1:nSplitOff)'));
        uid = splitOffUnits(1) - 1;
        splittingUnits = [uid; splitOffUnits];

        partitioning = arrayfun(@(i) find(scAfter == i), splittingUnits, 'UniformOutput', 0);

        history.optype{end+1} = 'split';
        history.indices{end+1} = {splittingUnits; partitioning};
        history.message{end+1} = sprintf('splitted %d -> %s', uid, ...
            jrclust.utils.field2str(splittingUnits));

        unitMask = ismember(after, splittingUnits);
        before(before > uid) = before(before > uid) + nSplitOff;
        before(unitMask) = after(unitMask);
    end

    diffMask = before ~= after;

    if all(diffMask)  % something else
        maxBefore = max(before);
        beforeAfter = arrayfun(@(iC) unique(after(before == iC)), ...
            1:maxBefore, 'UniformOutput', 0);

        if all(cellfun(@numel, beforeAfter) == 1) % a permutation
            history.optype{end+1} = 'reorder';
            history.indices{end+1} = [(1:maxBefore)' [beforeAfter{:}]'];
            history.message{end+1} = 'reorder units';
        else
            history.optype{end+1} = 'reassign';
            history.indices{end+1} = [scBefore(:) scAfter(:)];
            history.message{end+1} = 'reassign spikes';
        end
    end

    scBefore = scAfter;
end

fclose(fid);

% clean up
deleteFile('History has been converted. Delete the old history file?', obj.hCfg.histFile);

end % func
