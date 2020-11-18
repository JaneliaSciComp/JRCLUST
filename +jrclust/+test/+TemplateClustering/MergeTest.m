classdef MergeTest < jrclust.test.TemplateClustering.TemplateClusteringTestCase
    %MERGETEST Test merging of multiple clusters.

    %% HELPER METHODS
    methods
        function success = mergeMultiple(obj, unitIds)
            %MERGEMULTIPLE Merge several units in hClust.
            success = obj.hClust.mergeMultiple(unitIds);
        end

        function newTable = getSpikeTableAfterMergeMultiple(obj, unitIds)
            %GETSPIKETABLEAFTERMERGEMULTIPLE Get an updated spike table
            %when several units have been merged.
            newTable = obj.spikeClusters; % obj.spikeClusters is immutable
            
            % do this the dumb (but obviously correct) way
            unitIds = sort(unique(unitIds), 'descend');
            for i = 1:numel(unitIds)-1
                uid = unitIds(i);
                newTable(newTable == uid) = unitIds(end);
                newTable(newTable > uid) = newTable(newTable > uid) - 1;
            end
        end
    end

    %% TEST METHODS
    methods (Test)
        function mergeNonexistentUnitsErrors(obj)
            %MERGENONEXISTENTUNITCHANGESERRORS Ensure that an attempt to
            %merge a collection of units, all of which don't exist, will
            %throw an error.
            % no units exist
            obj.assertError(@() obj.hClust.mergeMultiple(obj.nClusters + [1, 2, 3]), ...
                ?MException);
            % some units don't exist
            obj.assertError(@() obj.hClust.mergeMultiple(obj.nClusters + [0 1 2]), ?MException);
            % noise/deleted units in merge
            obj.assertError(@() obj.hClust.mergeMultiple([0 1 2]), ?MException);
        end

        function mergeOneOrNoneErrors(obj)
            %MERGEONEORNONEERRORS Ensure that an error is thrown if
            %mergeMultiple is called with 0 or 1 units.
            
            % fails with no units
            obj.assertError(@() obj.hClust.mergeMultiple([]), ?MException);
            % fails with just one unit
            obj.assertError(@() obj.hClust.mergeMultiple(1), ?MException);
        end

        function mergeMultiplePropsChanged(obj)
            %MERGEMULTIPLEPROPSCHANGED Ensure that the number of edits, the
            %cluster count, and the history fields change after a merge.
            unitIds = [2; 4; 7; 8];
            obj.assertEqual(obj.mergeMultiple(unitIds), 1);

            % number of edits, number of clusters
            obj.assertEqual(obj.hClust.nEdits, 1);
            obj.assertEqual(obj.hClust.nClusters, obj.nClusters - numel(unitIds) + 1);

            % unit 2 needs its metadata recomputed
            obj.assertEqual(obj.hClust.recompute, unitIds(1));

            % history is updated to reflect units 2, 4, 7, and 8 having
            % been undeleted
            histEntry = obj.hClust.history;
            obj.assertEqual(histEntry.message{end}, 'merged [2; 4; 7; 8] -> 2');
            obj.assertEqual(histEntry.indices{end}, ...
                { unitIds; ...
                    { find(obj.spikeClusters == 2); ...
                      find(obj.spikeClusters == 4); ...
                      find(obj.spikeClusters == 7); ...
                      find(obj.spikeClusters == 8) } ...
                });
        end

        function mergeMultipleShiftDown(obj)
            %MERGEMULTIPLESHIFTDOWN Ensure that units that come after one
            %or more merged units in the spike table are shifted down by
            %the appropriate amount.
            unitIds = [2 4 7 8];
            newTable = obj.getSpikeTableAfterMergeMultiple(unitIds);
            obj.assertEqual(obj.mergeMultiple(unitIds), 1);

            obj.assertEqual(obj.hClust.spikeClusters, newTable);
        end

        function mergeMultipleVectorFieldsTruncated(obj)
            %MERGEMULTIPLEVECTORFIELDSTRUNCATED Ensure that vector fields
            %are properly truncated after a merge.
            unitIds = [1, 64, 10, 2]; % test an edge case
            obj.assertEqual(obj.mergeMultiple(unitIds), 1);

            % size is truncated
            obj.assertEqual(numel(obj.hClust.clusterNotes), obj.nClusters - 3);

            % proper values have been removed
            obj.assertTrue(ismember('1', obj.hClust.clusterNotes));
            obj.assertFalse(ismember('2', obj.hClust.clusterNotes));
            obj.assertFalse(ismember('10', obj.hClust.clusterNotes));
            obj.assertFalse(ismember('64', obj.hClust.clusterNotes));
        end

        function mergeMultipleMatrixFieldsTruncated(obj)
            %MERGEMULTIPLEMATRIXFIELDSTRUNCATED Ensure that the correct
            %entries in a matrix field have been truncated.
            unitIds = [64, 63, 14]; % test an edge case

            clusterCentroidsAfter = obj.hClust.clusterCentroids;
            clusterCentroidsAfter([63 64], :) = [];

            obj.assertEqual(obj.mergeMultiple(unitIds), 1);

            % size is truncated.
            obj.assertEqual(size(obj.hClust.clusterCentroids, 1), obj.nClusters - 2);
            % the correct entries are missing
            obj.assertEqual(obj.hClust.clusterCentroids, clusterCentroidsAfter);
        end
    end
end
