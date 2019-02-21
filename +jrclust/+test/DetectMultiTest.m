classdef DetectMultiTest < matlab.unittest.TestCase
    %DETECTMULTITEST Test spike detection from multiple files
    properties
        hJRC;
        hJRC0;
        hJRC1;
        nSpikes0 = 9891;
        nSpikes1 = 9000;
        nSpikes = 9891 + 9000;
    end

    properties (Dependent)
        res;
        res0;
        res1;
        hCfg;
        hCfg0;
        hCfg1;
    end

    methods (TestClassSetup)
        function setupProps(obj)
            hCfg_ = jrclust.Config(fullfile(getenv('JRCTESTDATA'), 'multi', 'test.prm'));
            obj.hJRC = JRC(hCfg_);
            obj.hCfg.testRun = 1;
            obj.hCfg.extractAfterDetect = 1;

            hCfg0_ = jrclust.Config(fullfile(getenv('JRCTESTDATA'), 'multi', 'test0.prm'));
            obj.hJRC0 = JRC(hCfg0_);
            obj.hCfg0.testRun = 1;
            obj.hCfg0.extractAfterDetect = 1;

            hCfg1_ = jrclust.Config(fullfile(getenv('JRCTESTDATA'), 'multi', 'test1.prm'));
            obj.hJRC1 = JRC(hCfg1_);
            obj.hCfg1.testRun = 1;
            obj.hCfg1.extractAfterDetect = 1;
        end

        function detectSingly(obj)
            obj.hJRC0.detect();
            obj.hJRC1.detect();
        end

        function detectTogether(obj)
            obj.hJRC.detect();
        end
    end

    methods (Test)
        function isSorted(obj)
            obj.assertTrue(issorted(obj.res.spikeTimes));
        end

        function consistentAcrossDetects(obj)
            % count of spikes in multiple detect should equal sum of counts spikes in single detects
            obj.assertEqual(numel(obj.res.spikeTimes), obj.nSpikes);
        end

        function allCountsEqual(obj)
            % vectors
            obj.assertEqual(numel(obj.res.spikeAmps), obj.nSpikes);
            obj.assertEqual(numel(obj.res.spikeSites), obj.nSpikes);
            obj.assertEqual(numel(obj.res.spikeSites2), obj.nSpikes);

            % matrices
            obj.assertEqual(size(obj.res.spikePositions, 1), obj.nSpikes);

            % tensors
            obj.assertEqual(size(obj.res.spikesRaw, 3), obj.nSpikes);
            obj.assertEqual(size(obj.res.spikesFilt, 3), obj.nSpikes);
            obj.assertEqual(size(obj.res.spikeFeatures, 3), obj.nSpikes);
        end

        function partitionOkay(obj)
            obj.assertEqual(sum(cellfun(@(c) numel(c), obj.res.spikesBySite)), obj.nSpikes);
            obj.assertEqual(sum(cellfun(@(c) numel(c), obj.res.spikesBySite2)), obj.nSpikes);

            nSites = numel(obj.res.spikesBySite);
            for i = 1:nSites
                iSpikes = obj.res.spikesBySite{i};
                iSpikes2 = obj.res.spikesBySite2{i};

                % no spike should have its secondary site equal to its primary site
                obj.assertEmpty(intersect(iSpikes, iSpikes2));
                for j = i+1:nSites
                    jSpikes = obj.res.spikesBySite{j};
                    jSpikes2 = obj.res.spikesBySite2{j};

                    obj.assertEmpty(intersect(iSpikes, jSpikes));
                    obj.assertEmpty(intersect(iSpikes2, jSpikes2));
                end
            end
        end

        function detIndependence(obj)
            % detect spikes at the same times in the first file
            firstHalfSame = all(obj.res.spikeTimes(1:obj.nSpikes0) == obj.res0.spikeTimes);
            obj.assertTrue(firstHalfSame);

            % detect spikes at the same times in the second file
            secondHalfOffset = unique(obj.res.spikeTimes(end-obj.nSpikes1+1:end) - obj.res1.spikeTimes);
            obj.assertEqual(numel(secondHalfOffset), 1);

            % detect the same amplitudes in each half
            firstHalfSame = all(obj.res.spikeAmps(1:obj.nSpikes0) == obj.res0.spikeAmps);
            secondHalfSame = all(obj.res.spikeAmps(end-obj.nSpikes1+1:end) == obj.res1.spikeAmps);
            obj.assertTrue(firstHalfSame);
            obj.assertTrue(secondHalfSame);

            % detect the same primary/secondary sites in each half
            firstHalfSame = all(obj.res.spikeSites(1:obj.nSpikes0) == obj.res0.spikeSites);
            secondHalfSame = all(obj.res.spikeSites(end-obj.nSpikes1+1:end) == obj.res1.spikeSites);
            obj.assertTrue(firstHalfSame);
            obj.assertTrue(secondHalfSame);

            firstHalfSame = all(obj.res.spikeSites2(1:obj.nSpikes0) == obj.res0.spikeSites2);
            secondHalfSame = all(obj.res.spikeSites2(end-obj.nSpikes1+1:end) == obj.res1.spikeSites2);
            obj.assertTrue(firstHalfSame);
            obj.assertTrue(secondHalfSame);

            % compute the same spike positions
            firstHalfNorm = norm(obj.res.spikePositions(1:obj.nSpikes0, :) - obj.res0.spikePositions);
            secondHalfNorm = norm(obj.res.spikePositions(end-obj.nSpikes1+1:end, :) - obj.res1.spikePositions);
            obj.assertEqual(firstHalfNorm, single(0));
            obj.assertEqual(secondHalfNorm, single(0));

            % extract the same raw/filtered spikes, features
            firstHalfSame = all(all(all(obj.res.spikesRaw(:, :, 1:obj.nSpikes0) - obj.res0.spikesRaw == 0)));
            secondHalfSame = all(all(all(obj.res.spikesRaw(:, :, end-obj.nSpikes1+1:end) - obj.res1.spikesRaw == 0)));
            obj.assertTrue(firstHalfSame);
            obj.assertTrue(secondHalfSame);

            firstHalfSame = all(all(all(obj.res.spikesFilt(:, :, 1:obj.nSpikes0) - obj.res0.spikesFilt == 0)));
            secondHalfSame = all(all(all(obj.res.spikesFilt(:, :, end-obj.nSpikes1+1:end) - obj.res1.spikesFilt == 0)));
            obj.assertTrue(firstHalfSame);
            obj.assertTrue(secondHalfSame);

            firstHalfSame = all(all(all(obj.res.spikeFeatures(:, :, 1:obj.nSpikes0) - obj.res0.spikeFeatures < 1e-8)));
            secondHalfSame = all(all(all(obj.res.spikeFeatures(:, :, end-obj.nSpikes1+1:end) - obj.res1.spikeFeatures < 1e-8)));
            obj.assertTrue(firstHalfSame);
            obj.assertTrue(secondHalfSame);
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

        function hCfg0 = get.hCfg0(obj)
            hCfg0 = obj.hJRC0.hCfg;
        end
        function set.hCfg0(obj, hCfg0)
            obj.hJRC0.hCfg = hCfg0;
        end

        function hCfg1 = get.hCfg1(obj)
            hCfg1 = obj.hJRC1.hCfg;
        end
        function set.hCfg1(obj, hCfg1)
            obj.hJRC1.hCfg = hCfg1;
        end

        function res = get.res(obj)
            res = obj.hJRC.res;
        end
        function res0 = get.res0(obj)
            res0 = obj.hJRC0.res;
        end
        function res1 = get.res1(obj)
            res1 = obj.hJRC1.res;
        end
    end
end

