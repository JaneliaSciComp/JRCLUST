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
            cRes = obj.hCurate.forceQuit();
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

            waveformSimBeforeMerge = obj.hClust.waveformSim;
            waveformSimBeforeMerge(isnan(waveformSimBeforeMerge)) = 0;

            % null automerge (cluster assignments shouldn't change)
            obj.hCurate.autoMerge(1);
            obj.assertEqual(obj.hClust.spikeClusters, spikeClusters);

            % revert (still no change expected)
            obj.hCurate.restoreHistory(1);
            obj.assertEqual(obj.hClust.spikeClusters, spikeClusters);

            obj.hCurate.autoMerge(0.95);

            waveformSim = obj.hClust.waveformSim;
            % clear nans and set diagonal to 0
            waveformSim = jrclust.utils.setDiag(waveformSim, zeros(size(waveformSim, 1)));
            waveformSim(isnan(obj.hClust.waveformSim)) = 0;
            obj.assertLessThanOrEqual(waveformSim, 0.95);

            % restore and check our waveformSim is restored as well
            obj.hCurate.restoreHistory(1);
            obj.assertEqual(obj.hClust.spikeClusters, spikeClusters);
            waveformSim = obj.hClust.waveformSim;
            waveformSim(isnan(obj.hClust.waveformSim)) = 0;
            obj.assertTrue(any(waveformSim(:) >= 0.95));

            obj.assertLessThan(norm(waveformSim - waveformSimBeforeMerge), 1);
        end

        function opMonteCarlo(obj)
            opsGroup = {'merge', 'split', 'delete'};
            operations = {};
            opClustersBefore = {};

            for i = 1:6
                op = randsample(opsGroup, 1);
                op = op{1}; % cell array
                operations{end+1} = op;
                spikeClustersAfter = obj.hClust.spikeClusters;
                opClustersBefore{end+1} = spikeClustersAfter;

                switch op
                    case 'merge'
                        iCluster = randsample(1:(obj.hClust.nClusters-1), 1);
                        jCluster = iCluster + 1;
                        obj.hCurate.updateSelect([iCluster, jCluster]);
                        obj.hCurate.mergeSelected();

                    case 'split'
                        iCluster = randsample(1:obj.hClust.nClusters, 1);
                        iSpikes = obj.hClust.spikesByCluster{iCluster};
                        retainedSpikes = randsample(iSpikes, floor(numel(iSpikes)/2));
                        obj.hCurate.splitCluster(iCluster, ismember(iSpikes, retainedSpikes));

                    case 'delete'
                        iCluster = randsample(1:obj.hClust.nClusters, 2); % delete 2 at a time
                        obj.hCurate.deleteClusters(iCluster);
                end
            end

            % still here? peel off operations and check our clusters match
            for i = 1:6
                obj.hCurate.restoreHistory(obj.hClust.nEdits);

                spikeClustersBefore = opClustersBefore{end};
                opClustersBefore = opClustersBefore(1:end-1);
                
                spikeClustersAfter = obj.hClust.spikeClusters;

                % if a cluster <= 0 it doesn't matter if it matches exactly
                goodClustersBefore = spikeClustersBefore > 0;
                goodClustersAfter = spikeClustersAfter > 0;
                obj.assertEqual(goodClustersBefore, goodClustersAfter);

                % but good clusters must match
                spikeClustersAfter = spikeClustersAfter(goodClustersAfter);
                spikeClustersBefore = spikeClustersBefore(goodClustersBefore);
                obj.assertEqual(spikeClustersAfter, spikeClustersBefore);
            end
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

