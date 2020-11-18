function success = mergeSelected(obj)
%MERGESELECTED Merge a pair of clusters.
success = 0;

if obj.isWorking
    jrclust.utils.qMsgBox('An operation is in progress.');
    return;
end

if numel(unique(obj.selected)) < 2
    return;
end

obj.isWorking = 1;

%% merge clusters
try
    success = obj.hClust.mergeMultiple(obj.selected);
catch ME
    success = 0;
    warning('Failed to merge: %s', ME.message);
    jrclust.utils.qMsgBox('Operation failed.');
end

obj.isWorking = 0;

if success
    success = obj.hClust.doRecompute();
end

if success
    mergeTarget = min(obj.selected);
    mergingUnit = max(obj.selected);

    % update showSubset
    showSubset = obj.showSubset;
    showSubset(showSubset == mergingUnit) = [];
    mask = showSubset > mergingUnit;
    showSubset(mask) = showSubset(mask) - 1;
    obj.showSubset = showSubset;

    % select newly-merged unit
    obj.selected = mergeTarget;

    % replot
    obj.updateHistMenu();
    obj.updateFigRD(); % centers changed, need replotting
    obj.replot();
else
    warning('Failed to recompute derived values.');
    jrclust.utils.qMsgBox('Operation failed.');
end
end %fun