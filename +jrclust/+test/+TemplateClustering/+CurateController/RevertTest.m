classdef RevertTest < jrclust.test.TemplateClustering.CurateController.CurateControllerTestCase
    %REVERTTEST Test reversion of one or more operations from the UI.

    %% HELPER METHODS
    methods
        function success = deleteClusters(obj, unitIds)
            %DELETECLUSTERS Delete a single unit from the UI.
            success = obj.hCurate.deleteClusters(unitIds);
        end

        function success = mergeTwo(obj, unitIds)
            %MERGETWO Select, then merge, two units in the UI.
            obj.hCurate.updateSelect(unitIds);
            success = obj.hCurate.mergeSelected();
        end

        function success = splitCluster(obj, unitId, partitioning)
            %SPLITCLUSTER Split a unit in the UI
            success = obj.hCurate.splitCluster(unitId, partitioning);
        end

        function success = revertLast(obj, n)
            %REVERTLAST Revert the last `n` operations in the UI.
            success = obj.hCurate.revertLast(n);
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
        function revertFirstOperationGivesInitialCommit(obj)
            %REVERTFIRSTOPERATIONGIVESINITIALCOMMIT Test that reverting the
            %very first operation gives only the 'Initial commit' menu
            %entry.
            hHistMenu = obj.hCurate.hMenus('HistMenu');

            % Clustering is in its initial state
            obj.assertEqual(numel(hHistMenu.Children), 1);

            initialCommit = hHistMenu.Children(1);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertTrue(initialCommit.Checked);

            % delete a cluster
            unitId = 1;
            obj.assertEqual(obj.deleteClusters(unitId), 1);
            obj.assertNotEqual(obj.hClust.spikeClusters, obj.spikeClusters);

            % investigate properties of new menu entry
            obj.assertEqual(numel(hHistMenu.Children), 2);
            deletedEntry = hHistMenu.Children(1);
            obj.assertEqual(deletedEntry.Text, sprintf('deleted %d', unitId));
            obj.assertTrue(deletedEntry.Checked);

            % revert
            obj.assertEqual(obj.revertLast(1), 1);
            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);

            % Clustering is back in its initial state
            obj.assertEqual(numel(hHistMenu.Children), 1);

            initialCommit = hHistMenu.Children(1);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertTrue(initialCommit.Checked);
        end

        function revertFirstOperationFromHistMenuGivesInitialCommit(obj)
            %REVERTFIRSTOPERATIONGIVESINITIALCOMMIT Test that reverting the
            %very first operation, this time from the history menu, gives
            %only the 'Initial commit' menu entry.
            hHistMenu = obj.hCurate.hMenus('HistMenu');

            % Clustering is in its initial state
            obj.assertEqual(numel(hHistMenu.Children), 1);

            initialCommit = hHistMenu.Children(1);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertTrue(initialCommit.Checked);

            % delete a cluster
            unitId = 1;
            obj.assertEqual(obj.deleteClusters(unitId), 1);
            obj.assertNotEqual(obj.hClust.spikeClusters, obj.spikeClusters);

            % investigate properties of new menu entry
            obj.assertEqual(numel(hHistMenu.Children), 2);
            deletedEntry = hHistMenu.Children(1);
            obj.assertEqual(deletedEntry.Text, sprintf('deleted %d', unitId));
            obj.assertTrue(deletedEntry.Checked);

            % revert from the history menu
            initialCommit = hHistMenu.Children(2);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertFalse(initialCommit.Checked);

            % invoke this menu entry's callback
            obj.assertEqual(initialCommit.MenuSelectedFcn(), 1);
            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);

            % Clustering is back in its initial state
            obj.assertEqual(numel(hHistMenu.Children), 1);

            initialCommit = hHistMenu.Children(1);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertTrue(initialCommit.Checked);
        end

        function revertSixOfTenOperations(obj)
            %REVERTSIXOFTENOPERATIONS Test that reverting six out of ten
            %operations gives the correct spike table.
            hHistMenu = obj.hCurate.hMenus('HistMenu');

            % Clustering is in its initial state
            obj.assertEqual(numel(hHistMenu.Children), 1);

            initialCommit = hHistMenu.Children(1);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertTrue(initialCommit.Checked);

            table4 = struct();

            rng(now);
            for i = 1:10
                % 1: delete; 2: merge; 3: split
                op = randi(3);

                switch op
                    case 1 % delete
                        unitId = randi(obj.hClust.nClusters);
                        obj.assertEqual(obj.deleteClusters(unitId), 1);

                        msg = sprintf('deleted %d', unitId);
                    case 2 % merge
                        unitIds = randsample(obj.hClust.nClusters, 2, 0);
                        obj.assertEqual(obj.mergeTwo(unitIds), 1);

                        msg = sprintf('merged [%d; %d] -> %d', min(unitIds), max(unitIds), min(unitIds));
                    case 3 % split
                        unitId = randi(obj.hClust.nClusters);
                        nPartitions = randi(4) + 1; % from 2 to 5 partitions
                        partitioning = obj.randomPartitioning(unitId, nPartitions);
                        obj.assertEqual(obj.splitCluster(unitId, partitioning), 1);

                        msg = sprintf('splitted %d -> [%s]', unitId, strjoin(arrayfun(@num2str, unitId:(unitId + nPartitions - 1), 'UniformOutput', 0), '; '));
                end

                % investigate properties of new menu entry
                obj.assertEqual(numel(hHistMenu.Children), min(5, i + 1));
                newestEntry = hHistMenu.Children(1);
                obj.assertEqual(newestEntry.Text, msg);
                obj.assertTrue(newestEntry.Checked);

                if i == 4
                    table4.spikeClusters = obj.hClust.spikeClusters;
                    table4.msg = msg;
                end
            end

            % perform the revert
            obj.assertEqual(obj.revertLast(6), 1);

            % after the revert, the spike table should be as we remembered
            % it
            obj.assertEqual(obj.hClust.spikeClusters, table4.spikeClusters);

            % inspect the state of the menu after the revert
            obj.assertEqual(numel(hHistMenu.Children), 5);
            entry = hHistMenu.Children(1);
            obj.assertEqual(entry.Text, table4.msg);
            obj.assertTrue(entry.Checked);

            % the very first (i.e., last) entry in the menu should be our
            % initial commit again
            obj.assertEqual(hHistMenu.Children(end).Text, 'Initial commit');
        end

        function revertFourOfTenOperationsFromHistMenu(obj)
            %REVERTFOUROFTENOPERATIONSFROMHISTMENU Test that reverting four
            %out of ten operations from the history menu gives the correct
            %spike table and history entries.
            hHistMenu = obj.hCurate.hMenus('HistMenu');

            % Clustering is in its initial state
            obj.assertEqual(numel(hHistMenu.Children), 1);

            initialCommit = hHistMenu.Children(1);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertTrue(initialCommit.Checked);

            table6 = struct();
            messages = cell(5, 1);

            rng(now);
            for i = 1:10
                % 1: delete; 2: merge; 3: split
                op = randi(3);

                switch op
                    case 1 % delete
                        unitId = randi(obj.hClust.nClusters);
                        obj.assertEqual(obj.deleteClusters(unitId), 1);

                        msg = sprintf('deleted %d', unitId);
                    case 2 % merge
                        unitIds = randsample(obj.hClust.nClusters, 2, 0);
                        obj.assertEqual(obj.mergeTwo(unitIds), 1);

                        msg = sprintf('merged [%d; %d] -> %d', min(unitIds), max(unitIds), min(unitIds));
                    case 3 % split
                        unitId = randi(obj.hClust.nClusters);
                        nPartitions = randi(4) + 1; % from 2 to 5 partitions
                        partitioning = obj.randomPartitioning(unitId, nPartitions);
                        obj.assertEqual(obj.splitCluster(unitId, partitioning), 1);

                        msg = sprintf('splitted %d -> [%s]', unitId, strjoin(arrayfun(@num2str, unitId:(unitId + nPartitions - 1), 'UniformOutput', 0), '; '));
                end

                % investigate properties of new menu entry
                obj.assertEqual(numel(hHistMenu.Children), min(5, i + 1));
                newestEntry = hHistMenu.Children(1);
                obj.assertEqual(newestEntry.Text, msg);
                obj.assertTrue(newestEntry.Checked);

                if i >= 2 && i <= 6
                    messages{i - 1} = msg;
                end
                if i == 6
                    table6.spikeClusters = obj.hClust.spikeClusters;
                    table6.msg = msg;
                end
            end

            % perform the revert from the history menu
            lastEntry = hHistMenu.Children(end);
            obj.assertEqual(lastEntry.MenuSelectedFcn(), 1);

            % after the revert, the spike table should be as we remembered
            % it
            obj.assertEqual(obj.hClust.spikeClusters, table6.spikeClusters);

            % inspect the state of the menu after the revert
            obj.assertEqual(numel(hHistMenu.Children), 5);
            entry = hHistMenu.Children(1);
            obj.assertEqual(entry.Text, table6.msg);
            obj.assertTrue(entry.Checked);

            % the history menu should show us entries corresponding to
            % Clustering history entries 2-6
            for i = 1:5
                child = hHistMenu.Children(6 - i); % i=1 => entry 5, i=2 => entry 4, &c.
                obj.assertEqual(child.Text, messages{i});
            end
        end
    end
end

