classdef RevertTest < jrclust.test.DensityPeakClustering.DensityPeakClusteringTestCase
    %REVERTTEST Test reverting of curation operations.

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
        function initialCommitOk(obj)
            fid = fopen(obj.histFile);
            checkInt = fread(fid, 1, 'int32');
            obj.assertEqual(checkInt, 1);

            committedClusters = fread(fid, obj.nSpikes, 'int32');
            obj.assertEqual(committedClusters, obj.sRes.spikeClusters);

            fclose(fid);
        end
    end
end

