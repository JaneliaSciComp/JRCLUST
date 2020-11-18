classdef ReorderClustersTest < jrclust.test.TemplateClustering.CurateController.CurateControllerTestCase
    %REORDERCLUSTERSTEST Test of cluster reordering from curation UI.
    
    %% HELPER METHODS
    methods
        function success = reorderClusters(obj, by)
            %REORDERCLUSTERS Reorder clusters by some criterion.
            if nargin < 2
                success = obj.hCurate.reorderClusters();
            else
                success = obj.hCurate.reorderClusters(by);
            end
        end

        function [newTable, afterIds] = getSpikeTableAfterReorderClusters(obj, by)
            %GETSPIKETABLEAFTERREORDERCLUSTERS Return a snapshot of what
            %the spike table should look like after reordering clusters.
            if nargin < 2
                by = 'clusterSites';
            end

            beforeIds = 1:obj.nClusters;
            if strcmp(by, 'Y + X')
                [~, afterIds] = sort(sum(obj.hClust.clusterCentroids, 2));
            else
                [~, afterIds] = sort(obj.hClust.(by));
            end

            newTable = obj.spikeClusters;
            for i = 1:numel(afterIds)
                newTable(obj.spikeClusters == beforeIds(i)) = afterIds(i);
            end
        end
    end

    %% TEST METHODS
    methods (Test)
        function reorderClustersWithNoArgsReordersBySite(obj)
            %REORDERCLUSTERSWITHNOARGSREORDERSBYSITE Test that reordering
            %clusters without specifying a criterion defaults to
            %'clusterSites'.
            [newTable, argsort] = obj.getSpikeTableAfterReorderClusters('clusterSites');
            spikesByClusterBefore = obj.hClust.spikesByCluster;

            obj.assertEqual(obj.reorderClusters(), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
            obj.assertEqual(obj.hClust.spikesByCluster, spikesByClusterBefore(argsort));
        end

        function reorderClustersBySite(obj)
            %REORDERCLUSTERSBYSITE Test that reordering clusters by site
            %gives the expected result.
            [newTable, argsort] = obj.getSpikeTableAfterReorderClusters('clusterSites');
            spikesByClusterBefore = obj.hClust.spikesByCluster;

            obj.assertEqual(obj.reorderClusters('clusterSites'), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
            obj.assertEqual(obj.hClust.spikesByCluster, spikesByClusterBefore(argsort));
        end

        function reorderClustersByCentroid(obj)
            %REORDERCLUSTERSBYCENTROID Test that reordering clusters by
            %centroid gives the expected result.
            [newTable, argsort] = obj.getSpikeTableAfterReorderClusters('Y + X');
            spikesByClusterBefore = obj.hClust.spikesByCluster;

            obj.assertEqual(obj.reorderClusters('Y + X'), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
            obj.assertEqual(obj.hClust.spikesByCluster, spikesByClusterBefore(argsort));
        end
    end
end

