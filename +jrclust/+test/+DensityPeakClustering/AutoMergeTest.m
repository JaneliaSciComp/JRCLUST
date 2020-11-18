classdef AutoMergeTest < jrclust.test.DensityPeakClustering.DensityPeakClusteringTestCase
    %AUTOMERGETEST Test automatic merge of clusters.

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
        function success = autoMerge(obj)
            %DELETESINGLE Delete a single unit in hClust.
            success = obj.hClust.autoMerge();
        end
    end

    %% TEST METHODS
    methods (Test)
        function autoMergeOK(obj)
            obj.assertEqual(obj.autoMerge(), 1);
            obj.verifyLessThan(obj.hClust.nClusters, obj.nClusters);
        end
    end
end
