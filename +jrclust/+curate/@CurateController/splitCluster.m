function splitCluster(obj, iCluster, unitPart)
    %SPLITCLUSTER Split off a cluster given retained spikes
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end
    
    if isempty(unitPart) || ~iscell(unitPart) || any(cellfun(@isempty, unitPart))
        return;
    end

    obj.isWorking = 1;

    showSubset = obj.showSubset;
    mask = showSubset > iCluster;
    try
        showSubset = [showSubset(~mask); iCluster+(1:numel(unitPart))'; showSubset(mask) + numel(unitPart)];
    catch
        showSubset = [showSubset(~mask) iCluster+(1:numel(unitPart)) showSubset(mask) + numel(unitPart)];        
    end

    try
        res = obj.hClust.splitUnit(obj.hClust.spikeClusters, iCluster, unitPart);
        if ~isempty(res.metadata)
            msg = sprintf('split %d -> %s', iCluster, strjoin(arrayfun(@num2str, iCluster+(0:numel(unitPart)), 'UniformOutput', 0), ', '));
            obj.hClust.commit(res.spikeClusters, res.metadata, msg);
            obj.showSubset = showSubset;
        end
    catch ME
        warning('Failed to split: %s', ME.message);
        obj.isWorking = 0;
        jrclust.utils.qMsgBox('Operation failed.');
        return;
    end
    obj.isWorking = 0;

    % replot
    obj.updateFigWav();
    obj.updateFigSim();
    obj.updateSelect([iCluster, iCluster + 1]);
end