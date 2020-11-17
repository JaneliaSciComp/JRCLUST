classdef ReorderTest < jrclust.test.TemplateClustering.TemplateClusteringTestCase
    %REORDERTEST Test reordering of clusters.

    %% HELPER METHODS
    methods
        function success = reorderMultiple(obj, beforeIds, afterIds)
            %REORDERMULTIPLE Reorder units in hClust.
            success = obj.hClust.reorderMultiple(beforeIds, afterIds);
        end

        function newTable = getSpikeTableAfterReorderMultiple(obj, beforeIds, afterIds)
            %GETSPIKETABLEAFTERREORDERMULTIPLE Get an updated spike table
            %when units have been reordered.
            newTable = obj.spikeClusters; % obj.spikeClusters is immutable
            
            for i = 1:numel(beforeIds)
                newTable(obj.spikeClusters == beforeIds(i)) = afterIds(i);
            end
        end
    end

    %% TEST METHODS
    methods (Test)
        function swapTwoUnitsUpdatesSpikeTable(obj)
            %SWAPTWOUNITSUPDATESSPIKETABLE Test that swapping out two units
            %results in the expected spike table.
            beforeIds = [1, 2];
            afterIds = [2, 1];
            
            newTable = obj.getSpikeTableAfterReorderMultiple(beforeIds, afterIds);
            
            obj.assertEqual(obj.reorderMultiple(beforeIds, afterIds), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
        end

        function oneCycleUnitsUpdatesSpikeTable(obj)
            %ONECYCLEUNITSUPDATESSPIKETABLE Test that cycling the table
            %forward by one results in the expected spike table.
            beforeIds = 1:obj.nClusters;
            afterIds = [2:obj.nClusters 1];
            
            newTable = obj.getSpikeTableAfterReorderMultiple(beforeIds, afterIds);
            
            obj.assertEqual(obj.reorderMultiple(beforeIds, afterIds), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
        end
% 
%         function deleteNoneOrMultipleErrors(obj)
%             %DELETENONEORMULTIPLEERRORS Ensure that an error is thrown if
%             %deleteSingle is called with 0 or > 1 units.
%             
%             % fails with no units
%             obj.assertError(@() obj.hClust.deleteSingle([]), ?MException);
%             % fails with more than one unit
%             obj.assertError(@() obj.hClust.deleteSingle([1, 2]), ?MException);
%         end
% 
%         function deleteSinglePropsChanged(obj)
%             %DELETESINGLPROPSCHANGED Ensure that the number of edits, the
%             %cluster count, and the history fields change after a delete.
%             unitId = 6;
%             obj.assertEqual(obj.deleteSingle(unitId), 1);
% 
%             obj.assertEqual(obj.hClust.nEdits, 1);
%             obj.assertEqual(obj.hClust.nClusters, obj.nClusters - 1);
% 
%             histEntry = obj.hClust.history;
%             obj.assertEqual(histEntry.message{end}, 'deleted 6');
%             obj.assertEqual(histEntry.indices{end}, [6, -1]);
%         end
% 
%         function deleteSingleShiftDown(obj)
%             %DELETESINGLESHIFTDOWN Ensure that units that come after a
%             %deleted unit in the spike table are shifted down by 1.
%             unitId = 9;
% 
%             newTable = obj.getSpikeTableAfterDeleteSingle(unitId);
%             obj.assertEqual(obj.deleteSingle(unitId), 1);
% 
%             obj.assertEqual(obj.hClust.spikeClusters, newTable);
%         end
% 
%         function deleteSingleVectorFieldsTruncated(obj)
%             %DELETESINGLEVECTORFIELDSTRUNCATED Ensure that vector fields
%             %are properly truncated after a delete.
%             unitId = 1; % test an edge case
%             obj.assertEqual(obj.deleteSingle(unitId), 1);
% 
%             % size is truncated
%             obj.assertEqual(numel(obj.hClust.clusterNotes), obj.nClusters - 1);
%             % proper value has been removed
%             obj.assertFalse(ismember(num2str(unitId), obj.hClust.clusterNotes));
%         end
% 
%         function deleteSingleMatrixFieldsTruncated(obj)
%             %DELETESINGLEMATRIXFIELDSTRUNCATED Ensure that the correct
%             %entries in a matrix field have been truncated.
%             unitId = obj.nClusters; % test an edge case
% 
%             clusterCentroidsBefore = obj.hClust.clusterCentroids;
%             obj.assertEqual(obj.deleteSingle(unitId), 1);
% 
%             % size is truncated.
%             obj.assertEqual(size(obj.hClust.clusterCentroids, 1), obj.nClusters - 1);
%             % the correct entries are missing
%             obj.assertEqual(obj.hClust.clusterCentroids, ...
%                             [clusterCentroidsBefore(1:unitId-1, :); clusterCentroidsBefore(unitId+1:end, :)]);
%         end
    end
end

