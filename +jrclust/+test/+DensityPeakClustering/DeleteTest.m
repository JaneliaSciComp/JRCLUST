classdef DeleteTest < jrclust.test.DensityPeakClustering.DensityPeakClusteringTestCase
    %DELETETEST Test deletion of one or many clusters.

    %% HELPER METHODS
    methods
        function success = deleteSingle(obj, unitId)
            %DELETESINGLE Delete a single unit in hClust.
            success = obj.hClust.deleteSingle(unitId);
        end

        function newTable = getSpikeTableAfterDeleteSingle(obj, unitId)
            %GETSPIKETABLEAFTERDELETESINGLE Get an updated spike table when
            %a unit has been deleted.
            newTable = obj.spikeClusters; % obj.spikeClusters is immutable
            newTable(newTable == unitId) = -1; % set to a negative value

            % units that come after the deleted unit in the spike table
            mask = (newTable > unitId);
            newTable(mask) = newTable(mask) - 1;
        end
    end

    %% TEST METHODS
    methods (Test)
        function deleteNonexistentUnitChangesNothing(obj)
            %DELETENONEXISTENTUNITCHANGESNOTHING Ensure that an attempt to
            %delete a unit that doesn't exist (or is considered a 'noise'
            %unit) will return early.
            obj.assertEqual(obj.deleteSingle(0), 1); % try to delete the noise unit

            obj.assertEqual(obj.hClust.nEdits, 0);
            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);
            
            obj.assertEqual(obj.deleteSingle(obj.nClusters + 1), 1); % try to delete a nonexistent cluster

            obj.assertEqual(obj.hClust.nEdits, 0);
            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);
        end

        function deleteNoneOrMultipleErrors(obj)
            %DELETENONEORMULTIPLEERRORS Ensure that an error is thrown if
            %deleteSingle is called with 0 or > 1 units.
            
            % fails with no units
            obj.assertError(@() obj.hClust.deleteSingle([]), ?MException);
            % fails with more than one unit
            obj.assertError(@() obj.hClust.deleteSingle([1, 2]), ?MException);
        end

        function deleteSinglePropsChanged(obj)
            %DELETESINGLPROPSCHANGED Ensure that the number of edits, the
            %cluster count, and the history fields change after a delete.
            unitId = 6;
            obj.assertEqual(obj.deleteSingle(unitId), 1);

            obj.assertEqual(obj.hClust.nEdits, 1);
            obj.assertEqual(obj.hClust.nClusters, obj.nClusters - 1);

            histEntry = obj.hClust.history;
            obj.assertEqual(histEntry.message{end}, 'deleted 6');
            obj.assertEqual(histEntry.indices{end}, [6, -1]);
        end

        function deleteSingleShiftDown(obj)
            %DELETESINGLESHIFTDOWN Ensure that units that come after a
            %deleted unit in the spike table are shifted down by 1.
            unitId = 9;

            newTable = obj.getSpikeTableAfterDeleteSingle(unitId);
            obj.assertEqual(obj.deleteSingle(unitId), 1);

            obj.assertEqual(obj.hClust.spikeClusters, newTable);
        end

        function deleteSingleVectorFieldsTruncated(obj)
            %DELETESINGLEVECTORFIELDSTRUNCATED Ensure that vector fields
            %are properly truncated after a delete.
            unitId = 1; % test an edge case
            obj.assertEqual(obj.deleteSingle(unitId), 1);

            % size is truncated
            obj.assertEqual(numel(obj.hClust.clusterNotes), obj.nClusters - 1);
            % proper value has been removed
            obj.assertFalse(ismember(num2str(unitId), obj.hClust.clusterNotes));
        end

        function deleteSingleMatrixFieldsTruncated(obj)
            %DELETESINGLEMATRIXFIELDSTRUNCATED Ensure that the correct
            %entries in a matrix field have been truncated.
            unitId = obj.nClusters; % test an edge case

            clusterCentroidsBefore = obj.hClust.clusterCentroids;
            obj.assertEqual(obj.deleteSingle(unitId), 1);

            % size is truncated.
            obj.assertEqual(size(obj.hClust.clusterCentroids, 1), obj.nClusters - 1);
            % the correct entries are missing
            obj.assertEqual(obj.hClust.clusterCentroids, ...
                            [clusterCentroidsBefore(1:unitId-1, :); clusterCentroidsBefore(unitId+1:end, :)]);
        end
    end
end

