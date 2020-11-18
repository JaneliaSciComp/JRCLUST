function success = splitCluster(obj, unitId, partitioning)
%SPLITCLUSTER Split a cluster according to a given partitioning.
%   The partitioning is expected in terms of *RELATIVE INDICES*.
success = 0;

if obj.isWorking
    jrclust.utils.qMsgBox('An operation is in progress.');
    return;
end

%% convert partitioning from relative to absolute indices
if nargin == 3 && iscell(partitioning)
    indices = find(obj.hClust.spikeClusters == unitId);

    if sum(cellfun(@numel, partitioning)) ~= numel(indices)
        warning('JRC:badPartitioning', 'Failed to split: incorrect number of spikes in partitioning.');
        jrclust.utils.qMsgBox('Operation failed.');
        return;
    end
else
    return;
end

try
    partitioning = cellfun(@(c) indices(c), partitioning, 'UniformOutput', 0);
catch ME
    warning('Failed to split: %s', ME.message);
    jrclust.utils.qMsgBox('Operation failed.');
    return;
end

unitIds = (unitId:(unitId + numel(partitioning) - 1))';

obj.isWorking = 1;

%% split cluster
try
    success = obj.hClust.splitSingle(unitIds, partitioning);
catch ME
    success = 0;
    warning('Failed to split: %s', ME.message);
    jrclust.utils.qMsgBox('Operation failed.');
end

obj.isWorking = 0;

if success
    success = obj.hClust.doRecompute();
end

if success
    % update showSubset
    showSubset = obj.showSubset(:);
    mask = showSubset > unitId;
    showSubset = [showSubset(~mask); unitId + (1:numel(partitioning))'; showSubset(mask) + numel(partitioning)];

    obj.showSubset = showSubset;

    obj.selected = [unitId, unitId + 1];

    % replot
    obj.updateHistMenu();
    obj.updateFigRD(); % centers changed, need replotting
    obj.replot();
end
end %fun