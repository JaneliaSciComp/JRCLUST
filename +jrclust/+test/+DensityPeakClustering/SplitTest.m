classdef SplitTest < jrclust.test.DensityPeakClustering.DensityPeakClusteringTestCase
    %SPLITTEST Test splitting of a single cluster.

    %% HELPER METHODS
    methods
        function success = splitSingle(obj, unitIds, partitioning)
            %SPLITSINGLE Split a single unit in hClust.
            success = obj.hClust.splitSingle(unitIds, partitioning);
        end
    end

    %% TEST METHODS
    methods (Test)
        function splitNonUniqueUnitsErrors(obj)
            %SPLITNONUNIQUEUNITSERRORS Ensure that an attempt to split a
            %unit into a collection of units which is not unique will throw
            %an error.
            obj.assertError(@() obj.splitSingle([1 1 2 3], cell(4)), ?MException);
        end

        function splitWithoutTargetUnitsErrors(obj)
            %SPLITWITHOUTTARGETUNITSERRORS Ensure that an attempt to split
            %a unit without another unit to split into will throw an error.
            % with 0 units should fail
            obj.assertError(@() obj.splitSingle([], cell(0)), ?MException);
            % with 1 unit should also fail
            obj.assertError(@() obj.splitSingle(1, cell(1)), ?MException);
        end

        function splitNonexistentUnitsErrors(obj)
            %SPLITNONEXISTENTUNITSERRORS Ensure that an attempt to split a
            %unit which is not found in the spike table will throw an
            %error.
            % splitting unit does not exist
            obj.assertError(@() obj.hClust.splitSingle(obj.nClusters + [1 2 3], cell(3)), ?MException);
            % include noise or deleted units in merge
            obj.assertError(@() obj.hClust.splitSingle([0 1 2], cell(3)), ?MException);
        end

        function partitioningSizeDoesNotMatchUnitsSize(obj)
            %PARTITIONINGSIZEDOESNOTMATCHUNITSSIZE Ensure that an attempt
            %to split a unit into a partitioning with a different number of
            %partitions than splitted units will throw an error.
            obj.assertError(@() obj.hClust.splitSingle([1 2], cell(1)), ?MException);
            obj.assertError(@() obj.hClust.splitSingle([1 2], cell(3)), ?MException);
        end
        
        function partitioningDoesNotAccountForAllSpikes(obj)
            %PARTITIONGDOESNOTACCOUNTFORALLSPIKES Ensure that an attempt to
            %split a unit into a partitioning which does not contain all
            %spikes, or contains too many spikes, or contains spikes
            %outside of the splitting unit, will throw an error.
            unitIds = [1 2];
            unit1Spikes = find(obj.hClust.spikeClusters == unitIds(1));
            unit2Spikes = find(obj.hClust.spikeClusters == unitIds(2));
            nSpikes = numel(unit1Spikes);

            tooSmallPartition = {unit1Spikes(1:floor(0.5*nSpikes)), ...
                                 unit1Spikes(floor(0.5*nSpikes)+1:floor(0.5*nSpikes) + 2)};
            tooLargePartition = {unit1Spikes, unit2Spikes};
            wrongUnitPartition = {unit1Spikes(1:floor(0.5*nSpikes)) + 1, ...
                                  unit1Spikes(floor(0.5*nSpikes):end) + 1};
            
            obj.assertError(@() obj.hClust.splitSingle(unitIds, tooSmallPartition), ?MException);
            obj.assertError(@() obj.hClust.splitSingle(unitIds, tooLargePartition), ?MException);
            obj.assertError(@() obj.hClust.splitSingle(unitIds, wrongUnitPartition), ?MException);
        end

        function splitSinglePropsChanged(obj)
            %SPLITSINGLEPROPSCHANGED Ensure that the number of edits, the
            %cluster count, and the history fields change after a split.
            unitIds = [2; 3; 5; 7];

            unitIndices = find(obj.spikeClusters == 2);
            nSpikes = numel(unitIndices);

            unitIndices = unitIndices(randperm(nSpikes));
            partitioning = {unitIndices(1:floor(0.25*nSpikes)), ...
                unitIndices(floor(0.25*nSpikes)+1:floor(0.5*nSpikes)), ...
                unitIndices(floor(0.5*nSpikes)+1:floor(0.75*nSpikes)), ...
                unitIndices(floor(0.75*nSpikes)+1:end)};
            
            obj.assertEqual(obj.splitSingle(unitIds, partitioning), 1);

            % number of edits, number of clusters
            obj.assertEqual(obj.hClust.nEdits, 1);
            obj.assertEqual(obj.hClust.nClusters, obj.nClusters + numel(unitIds) - 1);

            % all units should have metadata recomputed
            obj.assertEqual(obj.hClust.recompute, unitIds);

            % history is updated to reflect units 2, 4, 7, and 8 having
            % been undeleted
            histEntry = obj.hClust.history;
            obj.assertEqual(histEntry.message{end}, 'splitted 2 -> [2; 3; 5; 7]');
            obj.assertEqual(histEntry.indices{end}, {unitIds; partitioning});
        end

        function splitSingleShiftUp(obj)
            %SPLITSINGLESHIFTUP Ensure that the spike table is updated
            %correctly;
            unitIds = [2; 3; 5; 7];
            
            unitIndices = find(obj.spikeClusters == 2);
            nSpikes = numel(unitIndices);

            unitIndices = unitIndices(randperm(nSpikes));
            partitioning = {unitIndices(1:floor(0.25*nSpikes)), ...
                unitIndices(floor(0.25*nSpikes)+1:floor(0.5*nSpikes)), ...
                unitIndices(floor(0.5*nSpikes)+1:floor(0.75*nSpikes)), ...
                unitIndices(floor(0.75*nSpikes)+1:end)};
            
            obj.assertEqual(obj.splitSingle(unitIds, partitioning), 1); % successfully split
            
            % ensure unit 2 was properly split up
            obj.assertEqual(find(obj.hClust.spikeClusters == 2), sort(partitioning{1}));
            obj.assertEqual(find(obj.hClust.spikeClusters == 3), sort(partitioning{2}));
            obj.assertEqual(find(obj.hClust.spikeClusters == 5), sort(partitioning{3}));
            obj.assertEqual(find(obj.hClust.spikeClusters == 7), sort(partitioning{4}));
            
            % ensure everything else was shifted correctly
            mask = (obj.spikeClusters == 3); % 3->4
            obj.assertTrue(all(obj.hClust.spikeClusters(mask) == 4));
            mask = (obj.spikeClusters == 4); % 4->5->6
            obj.assertTrue(all(obj.hClust.spikeClusters(mask) == 6));
            mask = (obj.spikeClusters == 5); % 5->6->7->8
            obj.assertTrue(all(obj.hClust.spikeClusters(mask) == 8));
            % and so forth...
            for uid = 6:obj.nClusters
                mask = (obj.spikeClusters == uid);
                obj.assertTrue(all(obj.hClust.spikeClusters(mask) == uid + 3));
            end
        end

        function splitSingleVectorFieldsAugmented(obj)
            %SPLITSINGLEVECTORFIELDSAUGMENTED Ensure that vector fields
            %are properly augmented after a split.
            unitIds = [2; 3; 5; 7];
            
            unitIndices = find(obj.spikeClusters == 2);
            nSpikes = numel(unitIndices);

            unitIndices = unitIndices(randperm(nSpikes));
            partitioning = {unitIndices(1:floor(0.25*nSpikes)), ...
                unitIndices(floor(0.25*nSpikes)+1:floor(0.5*nSpikes)), ...
                unitIndices(floor(0.5*nSpikes)+1:floor(0.75*nSpikes)), ...
                unitIndices(floor(0.75*nSpikes)+1:end)};
            
            obj.assertEqual(obj.splitSingle(unitIds, partitioning), 1); % successfully split

            % size is augmented
            obj.assertEqual(numel(obj.hClust.clusterNotes), obj.nClusters + numel(unitIds) - 1);

            % proper values have been added
            obj.assertEqual(obj.hClust.clusterNotes{2}, '2');
            obj.assertEqual(obj.hClust.clusterNotes{3}, '');
            obj.assertEqual(obj.hClust.clusterNotes{4}, '3');
            obj.assertEqual(obj.hClust.clusterNotes{5}, '');
            obj.assertEqual(obj.hClust.clusterNotes{6}, '4');
            obj.assertEqual(obj.hClust.clusterNotes{7}, '');
            obj.assertEqual(obj.hClust.clusterNotes{8}, '5');
        end

        function splitSingleMatrixFieldsAugmented(obj)
            %SPLITSINGLEMATRIXFIELDSAUGMENTED Ensure that the correct
            %entries in a matrix field have been augmented.
            unitIds = [64; 65; 66]; % test an edge case

            clusterCentroidsBefore = obj.hClust.clusterCentroids;

            unitIndices = find(obj.spikeClusters == 64);
            nSpikes = numel(unitIndices);

            unitIndices = unitIndices(randperm(nSpikes));
            partitioning = {unitIndices(1:floor(0.33*nSpikes)), ...
                unitIndices(floor(0.33*nSpikes)+1:floor(0.67*nSpikes)), ...
                unitIndices(floor(0.67*nSpikes)+1:end)};
            
            obj.assertEqual(obj.splitSingle(unitIds, partitioning), 1); % successfully split

            % size is augmented.
            obj.assertEqual(size(obj.hClust.clusterCentroids, 1), obj.nClusters + 2);
            % everything that came before has been preserved
            obj.assertEqual(obj.hClust.clusterCentroids(1:obj.nClusters, :), clusterCentroidsBefore);
            % the new entries are zeros
            obj.assertTrue(all(all(obj.hClust.clusterCentroids(65:66, :) == 0)));
        end
    end
end
