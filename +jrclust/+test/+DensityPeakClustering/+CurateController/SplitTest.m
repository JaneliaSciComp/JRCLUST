classdef SplitTest < jrclust.test.DensityPeakClustering.CurateController.CurateControllerTestCase
    %SPLITTEST Test splitting of a cluster from the curation UI.
    
    %% HELPER METHODS
    methods
        function success = splitCluster(obj, unitId, partitioning)
            %SPLITCLUSTER Split a cluster.
            success = obj.hCurate.splitCluster(unitId, partitioning);
        end

        function newTable = getSpikeTableAfterSplitCluster(obj, unitId, partitioning)
            %GETSPIKETABLEAFTERSPLITCLUSTER Return a snapshot of what the
            %spike table should look like after splitting unitId.
            newTable = obj.spikeClusters; % obj.spikeClusters is immutable

            mask = (newTable > unitId);
            newTable(mask) = newTable(mask) + numel(partitioning) - 1;

            indices = find(newTable == unitId);
            for i = 1:numel(partitioning) - 1
                p = partitioning{i + 1};
                newTable(indices(p)) = unitId + i;
            end
        end

        function partitioning = randomPartitioning(obj, unitId, nPartitions)
            %RANDOMPARTITIONING Randomly partition a cluster's relative
            %indices into nPartitions partitions, making them as even as
            %possible.
            nSpikes = sum(obj.hClust.spikeClusters == unitId);
            indices = randsample((1:nSpikes)', nSpikes, 0); % permute spike indices
            partitioning = cell(nPartitions, 1);

            for i = 1:(nPartitions - 1)
                start = floor((i - 1) * nSpikes / nPartitions) + 1;
                stop = floor(i * nSpikes / nPartitions);
                partitioning{i} = indices(start:stop);
            end
            partitioning{end} = indices((stop+1):end);
        end
    end
    
    %% TEST METHODS
    methods (Test)
        function splitWithInsufficientArgumentsUnsuccessful(obj)
            %SPLITWITHINSUFFICIENTARGUMENTSUNSUCCESSFUL Test that calling
            %hCurate.splitCluster with 0 or 1 arguments returns with a
            %success flag of false.
            obj.assertEqual(obj.hCurate.splitCluster(), 0);
            obj.assertEqual(obj.hCurate.splitCluster(1), 0);
        end

        function splitWithNonCellPartitioningUnsuccessful(obj)
            %SPLITWITHNONCELLPARTITIONINGUNSUCCESSFUL Test that calling
            %hCurate.splitCluster with a partitioning that is not a cell
            %returns a false success value.
            obj.assertEqual(obj.splitCluster(1, []), 0);
        end

        function splitWithIncorrectSpikeCountUnsuccessful(obj)
            unitId = 1;
            partitioning = obj.randomPartitioning(unitId, 2);

            % chop off the last element from the last partition
            lastPartition = partitioning{end};
            partitioning{end} = lastPartition(1:end-1);

            success = obj.assertWarning(@() obj.splitCluster(unitId, partitioning), 'JRC:badPartitioning');
            obj.assertEqual(success, 0);
        end

        function splitFromUISplitsInClustering(obj)
            %SPLITFROMUISPLITSINCLUSTERING Test that performing a split in
            %the UI reflects that split in the underlying Clustering.
            unitId = 10;
            nPartitions = 3;
            partitioning = obj.randomPartitioning(unitId, nPartitions);
            obj.assertEqual(obj.splitCluster(unitId, partitioning), 1);

            newTable = obj.getSpikeTableAfterSplitCluster(unitId, partitioning);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);
        end

        function splitSelectsFirstTwoUnitIds(obj)
            %SPLITSELECTSFIRSTTWOUNITIDS Test that, after a split is
            %performed, the CurateController selects the originally
            %splitted unit and its immediate successor.
            unitId = obj.nClusters;
            nPartitions = 5;
            partitioning = obj.randomPartitioning(unitId, nPartitions);

            obj.assertEqual(obj.splitCluster(unitId, partitioning), 1);
            obj.assertEqual(obj.selected, [unitId, unitId + 1]);
        end

        function splitDoesRecompute(obj)
            %SPLITDOESRECOMPUTE Test that, after a split is
            %performed, doRecompute is called and recompute is empty.
            unitId = floor(obj.nClusters / 2);
            nPartitions = 2;
            partitioning = obj.randomPartitioning(unitId, nPartitions);

            obj.assertEqual(obj.splitCluster(unitId, partitioning), 1);
            obj.assertEmpty(obj.hClust.recompute);
        end

        function splitUpdatesHistMenu(obj)
            %SPLITUPDATESHISTMENU Test that a split updates the history
            %menu, and what the history menu looks like afterward.
            hHistMenu = obj.hCurate.hMenus('HistMenu');

            % Clustering is in its initial state
            obj.assertEqual(numel(hHistMenu.Children), 1);

            initialCommit = hHistMenu.Children(1);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertTrue(initialCommit.Checked);

            % split the unit and observe its effect on the history menu
            unitId = floor(obj.nClusters / 2);
            nPartitions = 2;
            partitioning = obj.randomPartitioning(unitId, nPartitions);
            obj.assertEqual(obj.splitCluster(unitId, partitioning), 1);

            % investigate properties of new menu entry
            obj.assertEqual(numel(hHistMenu.Children), 2);

            splittedEntry = hHistMenu.Children(1);
            obj.assertEqual(splittedEntry.Text, ...
                sprintf('splitted %d -> [%s]', unitId, strjoin(arrayfun(@num2str, unitId:(unitId + nPartitions - 1), 'UniformOutput', 0), '; ')));
            obj.assertTrue(splittedEntry.Checked);

            % effect on Initial commit entry is to push it down and set
            % Checked to 'off'
            initialCommit = hHistMenu.Children(2);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertFalse(initialCommit.Checked);
        end
    end
end

