classdef AutoDeleteTest < jrclust.test.DensityPeakClustering.CurateController.CurateControllerTestCase
    %AUTODELETETEST Test automatic deletion of clusters from the curation
    %UI.
    
    %% SETUP METHODS
    methods (TestMethodSetup)
        function updateUnitSNR(obj)
            %UPDATEUNITSNR Set unitSNR to allow roughly 1/3 of units to be
            %deleted.
            obj.hClust.unitSNR = 15 * rand(size(obj.hClust.unitSNR));
        end
    end

    %% HELPER METHODS
    methods
        function success = autoDelete(obj)
            %AUTODELETE Automatically delete clusters.
            success = obj.hCurate.autoDelete();
        end
    end
    
    %% TEST METHODS
    methods (Test)
        function autoDeleteOK(obj)
            obj.assertEqual(obj.autoDelete(), 1);
            obj.verifyLessThan(obj.hClust.nClusters, obj.nClusters);
            obj.assertEmpty(obj.hClust.recompute);
        end
    end
end

