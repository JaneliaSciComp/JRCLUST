function success = deleteClusters(obj, unitIds)
%DELETECLUSTERS Delete clusters, either specified or selected.
success = 0;

if obj.isWorking
    jrclust.utils.qMsgBox('An operation is in progress.');
    return;
end

if nargin < 2 && numel(obj.selected) > 1
    return;
elseif nargin < 2
    unitIds = obj.selected(1);
end

obj.isWorking = 1;

%% delete clusters
try
    success = obj.hClust.deleteMultiple(unitIds);
catch ME
    success = 0;
    warning('Failed to delete: %s', ME.message);
    jrclust.utils.qMsgBox('Operation failed.');
end

obj.isWorking = 0;

if success
    % update showSubset
    unitIds = sort(unitIds, 'desc');

    showSubset = obj.showSubset;
    for i = 1:numel(unitIds)
        iCluster = unitIds(i);
        showSubset(showSubset == iCluster) = [];
        mask = (showSubset > iCluster);
        showSubset(mask) = showSubset(mask) - 1;
    end
    obj.showSubset = showSubset;

    % fix OOB error when deleting last cluster in updateFigSim
    obj.selected = min([max(obj.showSubset), obj.selected]);

    % replot
    obj.updateHistMenu();
    obj.updateFigRD(); % centers changed, need replotting
    obj.replot();
else
    warning('Failed to recompute derived values.');
    jrclust.utils.qMsgBox('Operation failed.');
end
end %fun
