classdef ManualTest < matlab.unittest.TestCase
    %MANUALTEST Test manual UI
    properties
        hJRC;
    end

    properties (Dependent)
        hCfg;
        hClust;
        hCurate;
    end

    %% SETUP METHODS
    methods (TestClassSetup)
        function setupProps(obj)
            hCfg_ = jrclust.Config(fullfile(getenv('JRCTESTDATA'), 'large', 'test_large.prm'));
            hCfg_.outlierThresh = 0; % disable `rmOutlierSpikes`
            obj.hJRC = JRC(hCfg_);
            obj.hCfg.testRun = 1;
        end

        function doCurate(obj)
            obj.hJRC.curate();
        end
    end

    %% TEARDOWN METHODS
    methods (TestClassTeardown)
        function endSesh(obj)
            cRes = obj.hCurate.endSession();
            close all;
        end
    end
    
    methods (Test)
        function selectOkay(obj)
            obj.hCurate.updateSelect(obj.hClust.nClusters); % select last cluster
            lastCluster = obj.hCurate.selected;
            obj.hCurate.updateSelect(obj.hClust.nClusters + 1); % try to select OOB
            obj.assertEqual(obj.hCurate.selected, lastCluster);

            obj.hCurate.updateSelect(-10.5); % try to select OOB, noninteger
            obj.assertEqual(obj.hCurate.selected, 1);
        end

        function autoMergeOkay(obj)
            spikeClusters = obj.hClust.spikeClusters;

            % null automerge (cluster assignments shouldn't change)
            obj.hCurate.autoMerge(1);
            obj.assertEqual(obj.hClust.spikeClusters, spikeClusters);

            % revert (still no change expected)
            obj.hCurate.restoreHistory(1);
            obj.assertEqual(obj.hClust.spikeClusters, spikeClusters);
        end
    end

    %% GETTERS/SETTERS
    methods
        function hCfg = get.hCfg(obj)
            hCfg = obj.hJRC.hCfg;
        end
        function set.hCfg(obj, hCfg)
            obj.hJRC.hCfg = hCfg;
        end

        function hClust = get.hClust(obj)
            hClust = obj.hCurate.hClust;
        end

        function hCurate = get.hCurate(obj)
            hCurate = obj.hJRC.hCurate;
        end
    end
end

