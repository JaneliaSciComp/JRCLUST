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

            for i = 2:20
                % 1: delete; 2: undelete; 3: merge; 4: split
                op = randi(4);

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
                end

                fwrite(fid, int32(i), 'int32');
                fwrite(fid, int32(obj.hClust.spikeClusters), 'int32');
            end

            fclose(fid);
            obj.history = obj.hClust.history;
        end

        function setupResFile(obj)
            spikeClusters = obj.hClust.spikeClusters;
            save(obj.hCfg.resFile, 'spikeClusters');
        end

        function setupJRC(obj)
            obj.hJRC = JRC(obj.hCfg);
        end
    end

    %% TEST METHODS
    methods (Test)
        function ok(obj)
            hist = obj.hJRC.convertHistory();
            for i = 1:numel(obj.history.optype)
                obj.assertEqual(obj.history.optype{i}, hist.optype{i});
                obj.assertEqual(obj.history.message{i}, hist.message{i});
                obj.assertEqual(obj.history.indices{i}, hist.indices{i});
            end
        end
    end
end

