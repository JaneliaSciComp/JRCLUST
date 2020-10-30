classdef MergeTest < jrclust.test.CurateController.CurateControllerTestCase
    %MERGETEST Test merging of two clusters from the curation UI.
    
    %% HELPER METHODS
    methods
        function success = mergeSelected(obj)
            %MERGESELECTED Merge the two selected clusters.
            success = obj.hCurate.mergeSelected();
        end

        function newTable = getSpikeTableAfterMergeUnits(obj, unitIds)
            %GETSPIKETABLEAFTERMERGEUNITS Return a snapshot of what the
            %spike table should look like after merging unitIds.
            newTable = obj.spikeClusters; % obj.spikeClusters is immutable

            if numel(unitIds) == 2
                mergeTarget = min(unitIds);
                mergingUnit = max(unitIds);

                newTable(newTable == mergingUnit) = mergeTarget;
                newTable(newTable > mergingUnit) = newTable(newTable > mergingUnit) - 1;
            end
        end
    end

    %% TEST METHODS
    methods (Test)
        function mergeWithOneSelectedFails(obj)
            %MERGEWITHONESELECTEDFAILS Test that attempting to merge with
            %only one unit selected results in failure.
            obj.updateSelect(1);

            obj.assertEqual(obj.mergeSelected(), 0);
        end

        function mergeFromUIMergesInClustering(obj)
            %MERGEFROMUIMERGESINCLUSTERING Test that merging two units from
            %the UI merges the same two units in the Clustering.
            unitIds = [1, 2];
            obj.updateSelect(unitIds);

            obj.assertEqual(obj.mergeSelected(), 1);

            newTable = obj.getSpikeTableAfterMergeUnits(unitIds);
            obj.assertEqual(obj.hCurate.hClust.spikeClusters, newTable);
        end

        function mergeFromUISelectsMergeTarget(obj)
            %MERGEFROMUISELECTSMERGETARGET Test that merging two units from
            %the UI results in the merge target (the smaller unit) being
            %updated.
            unitIds = [obj.nClusters - 1, obj.nClusters - 2];
            obj.updateSelect(unitIds);

            obj.assertEqual(obj.mergeSelected(), 1);
            obj.assertEqual(obj.selected, min(unitIds));
        end

        function mergeDoesRecompute(obj)
            %MERGEDOESRECOMPUTE Test that on a successful merge,
            %doRecompute is called.
            unitIds = [obj.nClusters - 1, obj.nClusters - 2];
            obj.updateSelect(unitIds);

            obj.assertEqual(obj.mergeSelected(), 1);
            obj.assertEmpty(obj.hClust.recompute);
        end

        function mergeUpdatesHistMenu(obj)
            %MERGEUPDATESHISTMENU Test that a merge updates the history
            %menu, and what the history menu looks like afterwards.
            hHistMenu = obj.hCurate.hMenus('HistMenu');

            % Clustering is in its initial state
            obj.assertEqual(numel(hHistMenu.Children), 1);

            initialCommit = hHistMenu.Children(1);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertTrue(initialCommit.Checked);

            % merge the units and observe its effect on the history menu
            unitIds = [obj.nClusters - 1, obj.nClusters - 2];
            obj.updateSelect(unitIds);

            obj.assertEqual(obj.mergeSelected(), 1);

            % investigate properties of new menu entry
            obj.assertEqual(numel(hHistMenu.Children), 2);
            mergedEntry = hHistMenu.Children(1);
            obj.assertEqual(mergedEntry.Text, ...
                sprintf('merged [%d; %d] -> %d', min(unitIds), max(unitIds), min(unitIds)));
            obj.assertTrue(mergedEntry.Checked);

            % effect on Initial commit entry is to push it down and set
            % Checked to 'off'
            initialCommit = hHistMenu.Children(2);
            obj.assertEqual(initialCommit.Text, 'Initial commit');
            obj.assertFalse(initialCommit.Checked);
        end
    end
end

