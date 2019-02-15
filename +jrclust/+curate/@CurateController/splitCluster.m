function splitCluster(obj, iCluster, retainedSpikes)
    %SPLITCLUSTER Split off a cluster given retained spikes
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    iSpikes = obj.hClust.spikesByCluster{iCluster};

    obj.isWorking = 1;
    try
        [success, retainedSpikes] = obj.hClust.splitCluster(iCluster, iSpikes(retainedSpikes));
        if success
            % save the new clustering
            retainedSpikes = strjoin(arrayfun(@num2str, retainedSpikes, 'UniformOutput', 0), ',');
            commitMsg = sprintf('%s;split;%d;%s', datestr(now, 31), iCluster, retainedSpikes);
            obj.hClust.commit(commitMsg);

            % replot
            obj.updateFigWav();
            obj.updateFigSim();
            obj.updateSelect([iCluster, iCluster + 1]);
        else
            jrclust.utils.qMsgBox('Operation failed.');
        end
    catch ME
        warning('Failed to split: %s', ME.message);
        jrclust.utils.qMsgBox('Operation failed.');
    end

    obj.isWorking = 0;
end