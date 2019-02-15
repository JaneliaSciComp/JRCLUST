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
    try
        success = obj.hClust.mergeClusterPair(iCluster, jCluster);
        if success
            % save the new clustering
            commitMsg = sprintf('%s;merge;%d;%d', datestr(now, 31), iCluster, jCluster);
            obj.hClust.commit(commitMsg);

            % replot
            obj.updateFigWav();
            %obj.updateFigRD(); % centers changed, need replotting
            obj.updateFigSim();
            obj.updateSelect(min(obj.selected));
        else
            jrclust.utils.qMsgBox('Operation failed.');
        end
    catch ME
        warning('Failed to merge: %s', ME.message);
        jrclust.utils.qMsgBox('Operation failed.');
    end

    obj.isWorking = 0;
end