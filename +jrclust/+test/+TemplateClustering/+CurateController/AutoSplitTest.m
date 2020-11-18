classdef AutoSplitTest < jrclust.test.TemplateClustering.CurateController.CurateControllerTestCase
    %AUTOSPLITTEST Test automatic splitting of clusters from the curation
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
        function success = autoSplit(obj, multisite)
            %AUTOSPLIT Automatically split clusters.
            success = obj.hCurate.autoSplit(multisite);
        end
    end
    
    %% TEST METHODS
    methods (Test)
        function autoSplitMultiChanOK(obj)
            obj.assertEqual(obj.autoSplit(1), 1);
            obj.verifyGreaterThan(obj.hClust.nClusters, obj.nClusters);
            obj.assertEmpty(obj.hClust.recompute);
        end
    end
end

