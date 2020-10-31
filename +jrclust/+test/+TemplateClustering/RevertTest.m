classdef RevertTest < jrclust.test.TemplateClustering.TemplateClusteringTestCase
    %REVERTTEST Test reverting of curation operations.

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
        function setSeed(~)
            rng(now);
        end
    end

    %% TEST METHODS
    methods (Test)
        function revertDeleteIsUndelete(obj)
            %REVERTDELETEISUNDELETE Ensure that a reversion of a delete
            %operation is an undelete operation.
            unitId = randi(obj.nClusters);

            % perform the delete to revert
            obj.assertEqual(obj.deleteSingle(unitId), 1);
            obj.assertEqual(obj.hClust.nEdits, 1);

            % revert the delete
            obj.assertEqual(obj.revertLast(1), 1);
            obj.assertEqual(obj.hClust.nEdits, 0);

            % spikeClusters restored to original state
            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);

            % after revert, recompute is empty
            obj.assertEmpty(obj.hClust.recompute);
            
            % all fields are consistent
            obj.assertEmpty(obj.hClust.inconsistentFields());
        end

        function revertUndeleteIsDelete(obj)
            %REVERTUNDELETEISDELETE Ensure that a reversion of an undelete
            %operation is a delete operation.
            unitId = randi(obj.nClusters);

            % perform a delete to undelete later
            obj.assertEqual(obj.deleteSingle(unitId), 1);
            obj.assertEqual(obj.hClust.nEdits, 1);

            % revert the delete
            obj.assertEqual(obj.revertLast(1), 1);
            obj.assertEqual(obj.hClust.nEdits, 0);

            % spikeClusters restored to original state
            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);

            % after revert, recompute is empty
            obj.assertEmpty(obj.hClust.recompute);
            
            % all fields are consistent
            obj.assertEmpty(obj.hClust.inconsistentFields());
        end

        function revertMergeIsSplit(obj)
            %REVERTMERGEISSPLIT Ensure that a reversion of a merge
            %operation is a split operation.
            unitIds = randperm(obj.nClusters, 3)'; % merge 3 random units

            % perform the merge to revert
            obj.assertEqual(obj.mergeMultiple(unitIds), 1);
            obj.assertEqual(obj.hClust.nEdits, 1);

            % revert the merge
            obj.assertEqual(obj.revertLast(1), 1);
            obj.assertEqual(obj.hClust.nEdits, 0);

            % spikeClusters restored to original state
            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);

            % after revert, recompute is empty
            obj.assertEmpty(obj.hClust.recompute);
            
            % all fields are consistent
            obj.assertEmpty(obj.hClust.inconsistentFields());
        end

        function revertSplitIsMerge(obj)
            %REVERTSPLITISMERGE Ensure that a reversion of a split
            %operation is a merge operation.
            unitIds = randi(obj.nClusters - 2) + [0; 1; 2]; % split a unit into 3

            % create a new partitioning
            unitIndices = find(obj.spikeClusters == unitIds(1));
            nSpikes = numel(unitIndices);

            unitIndices = unitIndices(randperm(nSpikes));
            partitioning = {sort(unitIndices(1:floor(0.33*nSpikes))); ...
                sort(unitIndices(floor(0.33*nSpikes)+1:floor(0.67*nSpikes))); ...
                sort(unitIndices(floor(0.67*nSpikes)+1:end))};

            % perform the split to revert
            obj.assertEqual(obj.splitSingle(unitIds, partitioning), 1);
            obj.assertEqual(obj.hClust.nEdits, 1);

            % revert the split
            obj.assertEqual(obj.revertLast(1), 1);
            obj.assertEqual(obj.hClust.nEdits, 0);

            % spikeClusters restored to original state
            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);

            % after revert, recompute is empty
            obj.assertEmpty(obj.hClust.recompute);
            
            % all fields are consistent
            obj.assertEmpty(obj.hClust.inconsistentFields());
        end

        function revertSeveral(obj)
            %REVERTSEVERAL Perform one of each operation, then revert back
            %to the beginning.
            % first a delete
            unitId = randi(obj.hClust.nClusters);
            obj.assertEqual(obj.deleteSingle(unitId), 1);
            obj.assertEqual(obj.hClust.nEdits, 1);

            % now a second delete, to undelete immediately after
            unitId = randi(obj.hClust.nClusters);
            obj.assertEqual(obj.deleteSingle(unitId), 1);
            obj.assertEqual(obj.hClust.nEdits, 2);

            % recover deletedId (should be -2)
            indices = obj.hClust.history.indices{end};
            deletedId = indices(2);
            obj.assertEqual(deletedId, -2);

            % perform the undelete
            obj.assertEqual(obj.undeleteSingle(deletedId, unitId), 1);
            obj.assertEqual(obj.hClust.nEdits, 3);

            % now do a merge
            unitIds = sort(randperm(obj.hClust.nClusters, 4));
            obj.assertEqual(obj.mergeMultiple(unitIds), 1);
            obj.assertEqual(obj.hClust.nEdits, 4);

            % finally do a split
            unitIds = randi(obj.hClust.nClusters-4) + [0 1 2 3 4];

            % create the partitioning
            unitIndices = find(obj.hClust.spikeClusters == unitIds(1));
            nSpikes = numel(unitIndices);

            unitIndices = unitIndices(randperm(nSpikes));
            partitioning = { sort(unitIndices(1:floor(0.2*nSpikes))), ...
                             sort(unitIndices(floor(0.2*nSpikes)+1:floor(0.4*nSpikes))), ...
                             sort(unitIndices(floor(0.4*nSpikes)+1:floor(0.6*nSpikes))), ...
                             sort(unitIndices(floor(0.6*nSpikes)+1:floor(0.8*nSpikes))), ...
                             sort(unitIndices(floor(0.8*nSpikes)+1:end))
                           };

            obj.assertEqual(obj.splitSingle(unitIds, partitioning), 1);
            obj.assertEqual(obj.hClust.nEdits, 5);

            % now revert everything
            obj.assertEqual(obj.revertLast(5), 1);
            obj.assertEqual(obj.hClust.nEdits, 0);

            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);
            
            % after revert, recompute is empty
            obj.assertEmpty(obj.hClust.recompute);
            
            % all fields are consistent
            obj.assertEmpty(obj.hClust.inconsistentFields());
        end

        function crazyTown(obj)
            %CRAZYTOWN Perform a bunch of random merges, splits, deletes,
            %and undeletes. Revert them randomly. Ensure that our spike
            %table always matches up.
            for i = 1:20
                % 1: delete; 2: undelete; 3: merge; 4: split; 5: revert
                op = randi(5);
                
                % if we manage to have deleted everything, make sure we
                % revert on this iteration
                if obj.hClust.nClusters == 0
                    op = 5;
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
                    case 5 % revert
                        if obj.hClust.nEdits == 0
                            continue;
                        end

                        n = randi(obj.hClust.nEdits);
                        obj.assertEqual(obj.revertLast(n), 1);
                end
            end

            if obj.hClust.nEdits > 0
                obj.assertEqual(obj.revertLast(obj.hClust.nEdits), 1);
            end
            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);
            
            % after revert, recompute is empty
            obj.assertEmpty(obj.hClust.recompute);
            
            % all fields are consistent
            obj.assertEmpty(obj.hClust.inconsistentFields());
        end
    end
end

