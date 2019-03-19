function splitCluster(obj, iCluster, partition)
    %SPLITCLUSTER Split off a cluster given retained spikes
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    if ~iscell(partition)
        partition = {partition};
    end

    iSpikes = obj.hClust.spikesByCluster{iCluster};

    obj.isWorking = 1;
    try
        if numel(partition) == 1
            [success, retained] = obj.hClust.splitCluster(iCluster, iSpikes(partition{1}));

            if success
                % save the new clustering
                retained = strjoin(arrayfun(@num2str, retained, 'UniformOutput', 0), ',');
                commitMsg = sprintf('%s;split;%d;%s', datestr(now, 31), iCluster, retained);
                obj.hClust.commit(commitMsg);
            else
                obj.isWorking = 0;
                jrclust.utils.qMsgBox('Operation failed.');
                return;
            end
        else
            success = obj.hClust.partitionCluster(iCluster, partition);
            if success
                % save the new clustering
                partition = jrclust.utils.field2str(partition);
                commitMsg = sprintf('%s;partition;%d;%s', datestr(now, 31), iCluster, partition);
                obj.hClust.commit(commitMsg);
            else
                obj.isWorking = 0;
                jrclust.utils.qMsgBox('Operation failed.');
                return;
            end
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