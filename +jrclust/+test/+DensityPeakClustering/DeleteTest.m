classdef DeleteTest < jrclust.test.DensityPeakClustering.DensityPeakClusteringTestCase
    %DELETETEST Test deletion of one or many clusters.

    %% HELPER METHODS
    methods
        function deleteSingle(obj, unitId)
            %DELETESINGLE Delete a single unit in hClust.
            obj.hClust.deleteSingle(unitId);
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
        function deleteSingleEditsChanged(obj)
            %DELETESINGLEDITSCHANGED Ensure that both the number of edits
            %and the edit position changes after a delete.
            unitId = 6;

            obj.assertEqual(unique(obj.hClust.spikeClusters), (1:obj.nClusters)');
            obj.deleteSingle(unitId);

            obj.assertEqual(obj.hClust.nEdits, 2);
            obj.assertEqual(obj.hClust.editPos, 2);

            obj.resetClustering(); % clean up
        end

        function deleteSingleShiftDown(obj)
            %DELETESINGLESHIFTDOWN Ensure that units that come after a
            %deleted unit in the spike table are shifted down by 1.
            unitId = 9;

            obj.assertEqual(unique(obj.hClust.spikeClusters), (1:obj.nClusters)');
            newTable = obj.getSpikeTableAfterDeleteSingle(unitId);

            obj.deleteSingle(unitId);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);

            obj.resetClustering(); % clean up
        end

        function deleteSingleVectorFieldsTruncated(obj)
            %DELETESINGLEVECTORFIELDSTRUNCATED Ensure that vector fields
            %are properly truncated after a delete.
            unitId = 6;

            obj.assertEqual(unique(obj.hClust.spikeClusters), (1:obj.nClusters)');
            obj.deleteSingle(unitId);

            % size is truncated
            obj.assertEqual(numel(obj.hClust.clusterNotes), obj.nClusters - 1);
            % proper value has been removed
            obj.assertFalse(ismember(num2str(unitId), obj.hClust.clusterNotes));

            obj.resetClustering(); % clean up
        end

        function deleteSingleMatrixFieldsTruncated(obj)
            unitId = 9;

            obj.assertEqual(unique(obj.hClust.spikeClusters), (1:obj.nClusters)');
            clusterCentroidsBefore = obj.hClust.clusterCentroids;
            obj.deleteSingle(unitId);

            obj.assertEqual(obj.hClust.clusterCentroids, ...
                            [clusterCentroidsBefore(1:unitId-1, :); clusterCentroidsBefore(unitId+1:end, :)]);

            obj.resetClustering(); % clean up
        end
% 
%             % check that tensor fields are appropriately shrunk
%             obj.assertEqual(size(obj.hClust.meanWfGlobal, 3), obj.nClusters - 1);
    end
end

