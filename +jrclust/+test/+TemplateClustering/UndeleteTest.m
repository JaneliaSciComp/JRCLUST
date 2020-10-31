classdef UndeleteTest < jrclust.test.TemplateClustering.TemplateClusteringTestCase
    %UNDELETETEST Test undeletion of one cluster.
    
    %% HELPER METHODS
    methods
        function success = deleteSingle(obj, unitId)
            %DELETESINGLE Delete a single unit in hClust.
            success = obj.hClust.deleteSingle(unitId);
        end

        function success = undeleteSingle(obj, deletedId, newId)
            %UNDELETESINGLE Undelete a single unit in hClust.
            success = obj.hClust.undeleteSingle(deletedId, newId);
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
        function undeleteNonexistentErrors(obj)
            %UNDELETENONEXISTENTUNITCHANGESNOTHING Ensure that an attempt
            %to undelete a unit that already exists (or is considered a
            %'noise' unit, or is not found among the deleted entries) will
            %throw an error.
            obj.assertError(@() obj.hClust.undeleteSingle(0, 1), ?MException); % try to undelete the noise unit
        end

        function undeleteNoneOrMultipleErrors(obj)
            %DELETENONEORMULTIPLEERRORS Ensure that an error is thrown if
            %undeleteSingle is called with 0 or > 1 units.
            
            % fails with no units
            obj.assertError(@() obj.hClust.undeleteSingle([], 1), ?MException);
            obj.assertError(@() obj.hClust.undeleteSingle(-1, []), ?MException);
            % fails with more than one unit
            obj.assertError(@() obj.hClust.undeleteSingle([-1, -2], 1), ?MException);
            obj.assertError(@() obj.hClust.undeleteSingle(-1, [1, 2]), ?MException);
        end

        function undeleteSinglePropsChange(obj)
            %UNDELETESINGLWPROPSCHANGE Ensure that the number of edits, the
            %cluster count, and the history fields change after an undelete.
            unitId = 6;
            obj.assertEqual(obj.deleteSingle(unitId), 1);

            obj.assertEqual(obj.hClust.nEdits, 1);
            obj.assertEqual(obj.hClust.nClusters, obj.nClusters - 1);

            % now undelete
            obj.assertEqual(obj.undeleteSingle(min(obj.hClust.spikeClusters), unitId), 1);

            % number of edits, number of clusters
            obj.assertEqual(obj.hClust.nEdits, 2);
            obj.assertEqual(obj.hClust.nClusters, obj.nClusters);

            % unit 6 needs its metadata recomputed
            obj.assertEqual(obj.hClust.recompute, unitId);

            % history is updated to reflect unit 6 having been undeleted
            histEntry = obj.hClust.history;
            obj.assertEqual(histEntry.message{end}, 'undeleted 6');
            obj.assertEqual(histEntry.indices{end}, [-1, 6]);
        end

        function undeleteSingleShiftUp(obj)
            %UNDELETESINGLESHIFTUP Ensure that units that come after an
            %undeleted unit in the spike table are shifted up by 1.
            unitId = 9;

            newTable = obj.getSpikeTableAfterDeleteSingle(unitId);
            obj.assertEqual(obj.deleteSingle(unitId), 1);
            obj.assertEqual(obj.hClust.spikeClusters, newTable);

            % now undelete
            obj.assertEqual(obj.undeleteSingle(min(obj.hClust.spikeClusters), unitId), 1);

            % everything is as it was before
            obj.assertEqual(obj.hClust.spikeClusters, obj.spikeClusters);
        end

        function undeleteAddsNewUnitToRecompute(obj)
            %UNDELETEADDSNEWUNITTORECOMPUTE Ensure that a unit which is
            %undeleted is added to the array of values to recompute.
            unitId = 12;
            obj.assertEqual(obj.deleteSingle(unitId), 1);

            % now undelete
            obj.assertEqual(obj.undeleteSingle(min(obj.hClust.spikeClusters), unitId), 1);
        end

        function undeleteSingleVectorFieldsAugmented(obj)
            %UNDELETESINGLEVECTORFIELDSAUGMENTED Ensure that vector fields
            %are properly augmented after an delete.
            unitId = 1; % test an edge case
            unitISIRatioBefore = obj.hClust.unitISIRatio;

            obj.assertEqual(obj.deleteSingle(unitId), 1);

            % size is truncated
            obj.assertEqual(numel(obj.hClust.unitISIRatio), obj.nClusters - 1);
            % the correct value is missing
            obj.assertEqual(obj.hClust.unitISIRatio, [unitISIRatioBefore(1:unitId-1); unitISIRatioBefore(unitId+1:end)])

            % now undelete
            obj.assertEqual(obj.undeleteSingle(min(obj.hClust.spikeClusters), unitId), 1);
    
            % size is back to normal
            obj.assertEqual(numel(obj.hClust.unitISIRatio), obj.nClusters);

            % everything that wasn't deleted is back in place
            obj.assertEqual(obj.hClust.unitISIRatio(1:unitId-1), unitISIRatioBefore(1:unitId-1));
            obj.assertEqual(obj.hClust.unitISIRatio(unitId+1:end), unitISIRatioBefore(unitId+1:end));
        end

        function undeleteSingleMatrixFieldsAugmented(obj)
            %UNDELETESINGLEMATRIXFIELDSAUGMENTED Ensure that the correct
            %entries in a matrix field have been truncated.
            unitId = obj.nClusters; % test an edge case

            clusterCentroidsBefore = obj.hClust.clusterCentroids;
            obj.assertEqual(obj.deleteSingle(unitId), 1);

            % size is truncated.
            obj.assertEqual(size(obj.hClust.clusterCentroids, 1), obj.nClusters - 1);
            % the correct entries are missing
            obj.assertEqual(obj.hClust.clusterCentroids, ...
                            [clusterCentroidsBefore(1:unitId-1, :); clusterCentroidsBefore(unitId+1:end, :)]);

            % now undelete
            obj.assertEqual(obj.undeleteSingle(min(obj.hClust.spikeClusters), unitId), 1);
    
            % size is back to normal
            obj.assertEqual(size(obj.hClust.clusterCentroids, 1), obj.nClusters);

            % everything that wasn't deleted is back in place
            obj.assertEqual(obj.hClust.clusterCentroids(1:unitId-1, :), clusterCentroidsBefore(1:unitId-1, :));
            obj.assertEqual(obj.hClust.clusterCentroids(unitId+1:end, :), clusterCentroidsBefore(unitId+1:end, :));
        end
    end
end

