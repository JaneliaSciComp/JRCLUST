classdef AutoDeleteTest < jrclust.test.TemplateClustering.CurateController.CurateControllerTestCase
    %AUTODELETETEST Test automatic deletion of clusters from the curation
    %UI.

    %% SETUP METHODS
    methods (TestMethodSetup)
        function updateUnitVpp(obj)
            %UPDATEUNITVPP Set unitVpp to allow roughly 1/3 of units to be
            %deleted.
            obj.hClust.unitVpp = 60 * rand(size(obj.hClust.unitVpp));
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

