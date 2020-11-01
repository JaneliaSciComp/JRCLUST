classdef DeleteTest < jrclust.test.DensityPeakClustering.CurateController.CurateControllerTestCase
    %DELETETEST Test deletion of one or many clusters from the curation UI.
    
    %% HELPER METHODS
    methods
        function success = deleteClusters(obj, unitIds)
            %DELETECLUSTERS Delete one or more clusters.
            if nargin < 2
                success = obj.hCurate.deleteClusters();
            else
                success = obj.hCurate.deleteClusters(unitIds);
            end
        end

        function newTable = getSpikeTableAfterDeleteClusters(obj, unitIds)
            %GETSPIKETABLEAFTERDELETECLUSTERS Return a snapshot of what the
            %spike table should look like after deleting unitIds.
            unitIds = sort(unitIds, 'desc');
            newTable = obj.spikeClusters; % obj.spikeClusters is immutable
            for i = 1:numel(unitIds)
                iCluster = unitIds(i);
                newTable(newTable == iCluster) = -i;
                newTable(newTable > iCluster) = newTable(newTable > iCluster) - 1;
            end
        end
    end

    %% TEST METHODS
    methods (Test)
        function deleteWithNoArgumentsDeletesSelected(obj)
            %DELETEWITHNOARGUMENTSDELETESSELECTED Test whether a call to
            %deleteClusters with no argument deletes the currently selected
            %cluster.
            unitId = 7;
            newTable = obj.getSpikeTableAfterDeleteClusters(unitId);
            obj.updateSelect(unitId); % select the unit to be deleted

            obj.assertEqual(obj.deleteClusters(), 1);
            obj.assertEqual(obj.hCurate.hClust.spikeClusters, newTable);
        end

        function deleteSingleFromUIDeletesFromClustering(obj)
            %DELETESINGLEFROMUIDELETSFROMCLUSTERING Test whether deleting a
            %single unit from the UI deletes that same unit from the
            %Clustering.
            unitId = 1;
            newTable = obj.getSpikeTableAfterDeleteClusters(unitId);

            obj.assertEqual(obj.deleteClusters(unitId), 1);
            obj.assertEqual(obj.hCurate.hClust.spikeClusters, newTable);
        end

        function deleteUnselectedClusterDoesNotAffectSelectedCluster(obj)
            %DELETEUNSELECTEDCLUSTERDOESNOTAFFECTSELECTEDCLUSTER Test
            %whether deleting a cluster which is not selected changes the
            %selected cluster in any way.
            selectedUnitId = 4;
            deletedUnitId = 11;

            % select this unit; it should remain unchanged after deletion
            obj.updateSelect(selectedUnitId);

            obj.assertEqual(obj.deleteClusters(deletedUnitId), 1);
            obj.assertEqual(obj.selected(1), selectedUnitId);
        end

        function deleteSelectedLastClusterSelectsNextLastCluster(obj)
            %DELETESELECTEDLASTCLUSTERSELECTSNEXTLASTCLUSTER Test that,
            %having selected the last cluster, deleting that cluster
            %updates the selection to the formerly next-to-last cluster,
            %which becomes the new last cluster.
            unitId = obj.nClusters;

            obj.updateSelect(unitId);

            % deletes selected cluster, i.e., the last one
            obj.assertEqual(obj.deleteClusters(), 1)
            obj.assertEqual(obj.selected(1), obj.nClusters - 1);
        end

        function deleteDoesNotRequireRecompute(obj)
            %DELETEDOESNOTEREQUIRERECOMPUTE Test that on a successful
            %delete, recompute is empty.
            unitIds = [obj.nClusters - 1, obj.nClusters - 2];

            obj.assertEqual(obj.deleteClusters(unitIds), 1);
            obj.assertEmpty(obj.hClust.recompute);
        end

        function deleteUpdatesHistMenu(obj)
            %DELETEUPDATESHISTMENU Test that a deletion updates the
            %history menu, and what the history menu looks like afterwards.
            unitId = 10;

            hHistMenu = obj.hCurate.hMenus('HistMenu');

            % Clustering is in its initial state
            obj.assertEqual(numel(hHistMenu.Children), 1);

            initialCommit = hHistMenu.Children(1);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertTrue(initialCommit.Checked);

            % delete the unit and observe its effect on the history menu
            obj.assertEqual(obj.deleteClusters(unitId), 1);

            % investigate properties of new menu entry
            obj.assertEqual(numel(hHistMenu.Children), 2);
            deletedEntry = hHistMenu.Children(1);
            obj.assertEqual(deletedEntry.Text, sprintf('deleted %d', unitId));
            obj.assertTrue(deletedEntry.Checked);

            % effect on Initial commit entry is to push it down and set
            % Checked to 'off'
            initialCommit = hHistMenu.Children(2);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertFalse(initialCommit.Checked);
        end
    end
end

