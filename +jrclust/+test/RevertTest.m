classdef RevertTest < jrclust.test.Clustering.ClusteringTestCase
    %REVERTTEST Test reverting of curation operations.

    %% SETUP METHODS
    methods (TestClassSetup)
        
        function setupProps(obj)
            setupProps@jrclust.test.Clustering.ClusteringTestCase(obj);            
        end
    end

    %% TEARDOWN METHODS
    methods (TestClassTeardown)
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

