classdef AutoMergeTest < jrclust.test.TemplateClustering.CurateController.CurateControllerTestCase
    %AUTOMERGETEST Test automatic merging of clusters from the curation UI.

    %% SETUP METHODS
    methods (TestMethodSetup)
        function setupWaveformSim(obj)
            %SETUPWAVEFORMSIM Give our units higher similarity scores so
            %merges will happen.
            sim = rand(obj.nClusters);
            sim = min(sim + sim', 1.09) - 0.1;
            obj.hClust.waveformSim = sim;
        end
    end

    %% HELPER METHODS
    methods
        function success = autoMerge(obj, maxUnitSim)
            %AUTOMERGE Automatically merge clusters.
            if nargin < 2
                success = obj.hCurate.autoMerge();
            else
                success = obj.hCurate.autoMerge(maxUnitSim);
            end
        end
    end
    
    %% TEST METHODS
    methods (Test)
        function autoMergeDefaultMaxUnitSimOK(obj)
            %AUTOMERGEDEFAULTMAXUNITSIMOK Test that, using the default
            %maxUnitSim, autoMerge succeeds from the curation UI.
            obj.assertEqual(obj.autoMerge(), 1);
            obj.verifyLessThan(obj.hClust.nClusters, obj.nClusters);
            obj.assertEmpty(obj.hClust.recompute);
        end
    end
end

