function mergeSelected(obj)
    %MERGESELECTED Merge a pair of clusters
    if obj.isWorking
        return;
    end
    obj.isWorking = 1;

    if numel(obj.selected) < 2
        return;
    end
    iCluster = min(obj.selected);
    jCluster = max(obj.selected);

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
    end

    obj.isWorking = 0;
end