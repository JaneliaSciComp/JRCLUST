function mergeSelected(obj)
    %MERGESELECTED Merge a pair of clusters
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    if numel(obj.selected) < 2
        return;
    end

    iCluster = min(obj.selected);
    jCluster = max(obj.selected);

    obj.isWorking = 1;

    % speculatively merge clusters
    res = obj.hClust.mergeUnits(obj.hClust.spikeClusters, iCluster, jCluster);
    % operation found to be inconsistent
    if isempty(res.metadata)
        warning('failed to merge units %d and %d', iCluster, jCluster);
        obj.isWorking = 0;

        return;
    end

    msg = sprintf('merge %d and %d', iCluster, jCluster);
    try
        obj.hClust.commit(res(end).spikeClusters, res(end).metadata, msg);
    catch ME
        warning('Failed to merge: %s', ME.message);
        jrclust.utils.qMsgBox('Operation failed.');
    end

    obj.isWorking = 0; % in case updateSelect needs to zoom

    obj.selected = iCluster; % fix OOB error

    % replot
    obj.updateFigWav();
    obj.updateFigRD(); % centers changed, need replotting
    obj.updateFigSim();
    obj.updateSelect(iCluster);
end