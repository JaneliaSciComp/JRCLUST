classdef ReorderClustersTest < jrclust.test.DensityPeakClustering.CurateController.CurateControllerTestCase
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

        function newTable = getSpikeTableAfterReorderClusters(obj, by)
            %GETSPIKETABLEAFTERREORDERCLUSTERS Return a snapshot of what
            %the spike table should look like after reordering clusters.
            if nargin < 2
                by = 'clusterSites';
            end

            beforeIds = 1:obj.nClusters;
            switch by
                case 'clusterSites'
                    [~, afterIds] = sort(obj.hClust.clusterSites);
                case 'Y + X'
                    [~, afterIds] = sort(sum(obj.hClust.clusterCentroids, 2));
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
            newTable = obj.getSpikeTableAfterReorderClusters('clusterSites');

            obj.assertEqual(obj.reorderClusters(), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
        end

        function reorderClustersBySite(obj)
            %REORDERCLUSTERSBYSITE Test that reordering clusters by site
            %gives the expected result.
            newTable = obj.getSpikeTableAfterReorderClusters('clusterSites');

            obj.assertEqual(obj.reorderClusters('clusterSites'), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
        end

        function reorderClustersByCentroid(obj)
            %REORDERCLUSTERSBYCENTROID Test that reordering clusters by
            %centroid gives the expected result.
            newTable = obj.getSpikeTableAfterReorderClusters('Y + X');

            obj.assertEqual(obj.reorderClusters('Y + X'), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
        end
    end
end

