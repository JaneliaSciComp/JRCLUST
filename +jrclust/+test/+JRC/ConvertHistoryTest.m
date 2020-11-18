classdef ConvertHistoryTest < jrclust.test.DensityPeakClustering.DensityPeakClusteringTestCase
%CONVERTHISTORYTEST Test conversion of old-style history to new-style.

%% FIRST-CLASS PROPS
properties
    hJRC;
    history;
end

%% HELPER METHODS
methods
    function success = deleteSingle(obj, unitId)
        %DELETESINGLE Delete a single unit in hClust.
        success = obj.hClust.deleteSingle(unitId);
    end

    function success = undeleteSingle(obj, deletedId, newId)
        %UNDELETESINGLE Undelete a single unit in hClust.
        success = obj.hClust.undeleteSingle(deletedId, newId);
    end

    function success = mergeMultiple(obj, unitIds)
        %MERGEMULTIPLE Merge several units in hClust.
        success = obj.hClust.mergeMultiple(unitIds);
    end

    function success = splitSingle(obj, unitIds, partitioning)
        %SPLITSINGLE Split a single unit in hClust.
        success = obj.hClust.splitSingle(unitIds, partitioning);
    end

    function success = reassignAll(obj, beforeTable, afterTable)
        %REASSIGNALL Reassign all spikes in the spike table.
        success = obj.hClust.reassignAll(beforeTable, afterTable);
    end

    function success = revertLast(obj, n)
        %REVERTLAST Revert the last `n` operations in hClust.
        success = obj.hClust.revertLast(n);
    end
end

%% SETUP METHODS
methods (TestClassSetup)
    function setupHistory(obj)
        %SETUPHISTORY Perform a bunch of random merges, splits,
        %deletes, and undeletes, saving them to our history file.
        fid = fopen(obj.histFile, 'w');
        fwrite(fid, int32(1), 'int32');
        fwrite(fid, int32(obj.hClust.spikeClusters), 'int32');

        rng(now);
        for i = 2:20
            % 1: delete; 2: undelete; 3: merge; 4: split; 5: reassign all
            op = randi(5);

            % if we manage to have deleted everything, revert the last
            % thing and finish up.
            if obj.hClust.nClusters == 0
                obj.revertLast(1);
                fwrite(fid, int32(i+1), 'int32');
                fwrite(fid, int32(obj.hClust.spikeClusters), 'int32');
                break;
            end

            switch op
                case 1 % delete
                    unitId = randi(obj.hClust.nClusters);
                    obj.assertEqual(obj.deleteSingle(unitId), 1);
                case 2 % undelete
                    % get a unit to undelete (check it hasn't already
                    % been undeleted!
                    undeletes = find(strcmp('undelete', obj.hClust.history.optype));
                    if isempty(undeletes)
                        lastUndelete = 0;
                    else
                        lastUndelete = max(undeletes);
                    end

                    deletes = find(strcmp('delete', obj.hClust.history.optype));
                    deletes(deletes < lastUndelete) = [];
                    if isempty(deletes)
                        continue;
                    end

                    if numel(deletes) == 1
                        idx = deletes;
                    else
                        idx = randsample(deletes, 1);
                    end

                    indices = obj.hClust.history.indices{idx};
                    newId = indices(1); deletedId = indices(2);
                    obj.assertEqual(obj.undeleteSingle(deletedId, newId), 1);
                case 3 % merge
                    nMerge = min(randi(7) + 1, obj.hClust.nClusters - 1); % merge anywhere from 2 to 8 units
                    unitIds = randsample(obj.hClust.nClusters, nMerge);
                    obj.assertEqual(obj.mergeMultiple(unitIds), 1);
                case 4 % split
                    nSplit = randi(4) + 1;
                    unitIds = randi(obj.hClust.nClusters) + (0:nSplit-1);

                    % create the partitioning
                    unitIndices = find(obj.hClust.spikeClusters == unitIds(1));
                    nSpikes = numel(unitIndices);

                    unitIndices = unitIndices(randperm(nSpikes));
                    partitioning = cell(nSplit, 1);
                    for j = 1:nSplit-1
                        start = floor(nSpikes * (j-1)/nSplit) + 1;
                        stop = floor(nSpikes * j/nSplit);
                        partitioning{j} = sort(unitIndices(start:stop));
                    end
                    start = floor(nSpikes * (nSplit - 1)/nSplit) + 1;
                    partitioning{end} = sort(unitIndices(start:end));

                    obj.assertEqual(obj.splitSingle(unitIds, partitioning), 1);
                case 5 % reassign all
                    beforeTable = obj.hClust.spikeClusters;
                    afterTable = randsample(beforeTable, obj.nSpikes, 1);

                    obj.assertEqual(obj.reassignAll(beforeTable, afterTable), 1);
            end

            fwrite(fid, int32(i), 'int32');
            fwrite(fid, int32(obj.hClust.spikeClusters), 'int32');
        end

        fclose(fid);

        obj.hClust.doRecompute();
        obj.history = obj.hClust.history;
    end

    function saveRes(obj)
        res_ = struct();
        skipFields = {'spikesRaw', 'spikesRawVolt', 'spikesFilt', 'spikesFiltVolt', ...
            'spikesFilt2', 'spikesFilt3', 'spikeFeatures', 'hClust', 'hRecs'};

        % load dRes fields
        fieldNames = fieldnames(obj.dRes);
        for i = 1:numel(fieldNames)
            fn = fieldNames{i};
            if ismember(fn, skipFields)
                continue;
            end

            res_.(fn) = obj.dRes.(fn);
        end

        % load sRes fields
        fieldNames = fieldnames(obj.sRes);
        for i = 1:numel(fieldNames)
            fn = fieldNames{i};
            if ismember(fn, skipFields)
                continue;
            end

            res_.(fn) = obj.sRes.(fn);
        end

        % load hClust fields
        fieldNames = obj.hClust.getSaveFields(); % obj.hClust -> obj.res.hClust
        for i = 1:numel(fieldNames)
            fn = fieldNames{i};
            res_.(fn) = obj.hClust.(fn);
        end

        jrclust.utils.saveStruct(res_, obj.hCfg.resFile);
    end

    function setupJRC(obj)
        obj.hJRC = JRC(obj.hCfg);
        obj.hJRC.loadRes();
    end
end

%% TEST METHODS
methods (Test)
    function ok(obj)
        obj.hJRC.convertHistory();
        for i = 1:numel(obj.history.optype)
            obj.assertEqual(obj.history.optype{i}, obj.hClust.history.optype{i});
            obj.assertEqual(obj.history.message{i}, obj.hClust.history.message{i});
            obj.assertEqual(obj.history.indices{i}, obj.hClust.history.indices{i});
        end
    end
end
end % classdef
