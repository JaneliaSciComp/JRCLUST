classdef ReorderTest < jrclust.test.DensityPeakClustering.DensityPeakClusteringTestCase
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
            spikesByClusterBefore = obj.hClust.spikesByCluster;
            
            obj.assertEqual(obj.reorderMultiple(beforeIds, afterIds), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
            obj.assertEqual(obj.hClust.spikesByCluster{1}, spikesByClusterBefore{2});
            obj.assertEqual(obj.hClust.spikesByCluster{2}, spikesByClusterBefore{1});
        end

        function oneCycleUnitsUpdatesSpikeTable(obj)
            %ONECYCLEUNITSUPDATESSPIKETABLE Test that cycling the table
            %forward by one results in the expected spike table.
            beforeIds = 1:obj.nClusters;
            afterIds = [2:obj.nClusters 1];
            
            newTable = obj.getSpikeTableAfterReorderMultiple(beforeIds, afterIds);
            spikesByClusterBefore = obj.hClust.spikesByCluster;
            
            obj.assertEqual(obj.reorderMultiple(beforeIds, afterIds), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
            obj.assertEqual(obj.hClust.spikesByCluster, spikesByClusterBefore(afterIds));
        end
    end
end

