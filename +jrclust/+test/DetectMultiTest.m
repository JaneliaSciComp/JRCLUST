classdef DetectMultiTest < matlab.unittest.TestCase
    %DETECTMULTITEST

    properties
        dRes0;
        dRes1;
        dRes;
        testConfig;
        nSpikes0 = 11208;
        nSpikes1 = 9635;
        nSpikes = 11208 + 9635;
    end

    methods (TestClassSetup)
        function detectSingly(testCase)
            tc0 = fullfile(getenv('JRCTESTDATA'), 'multi', 'test0.prm');
            hJRC = jrc('detect', tc0);
            testCase.dRes0 = hJRC.dRes;

            tc1 = fullfile(getenv('JRCTESTDATA'), 'multi', 'test1.prm');
            hJRC = jrc('detect', tc1);
            testCase.dRes1 = hJRC.dRes;
        end

        function detectTogether(testCase)
            testCase.testConfig = fullfile(getenv('JRCTESTDATA'), 'multi', 'test.prm');
            hJRC = jrc('detect', testCase.testConfig);
            testCase.dRes = hJRC.dRes;
        end
    end

    methods (Test)
        function isSorted(testCase)
            testCase.assertTrue(issorted(testCase.dRes.spikeTimes));
        end

        function consistentAcrossDetects(testCase)
            % count of spikes in multiple detect should equal sum of counts spikes in single detects
            testCase.assertEqual(numel(testCase.dRes.spikeTimes), testCase.nSpikes);
        end

        function allCountsEqual(testCase)
            % vectors
            testCase.assertEqual(numel(testCase.dRes.spikeAmps), testCase.nSpikes);
            testCase.assertEqual(numel(testCase.dRes.spikeSites), testCase.nSpikes);
            testCase.assertEqual(numel(testCase.dRes.spikeSites2), testCase.nSpikes);

            % matrices
            testCase.assertEqual(size(testCase.dRes.spikePositions, 1), testCase.nSpikes);

            % tensors
            testCase.assertEqual(size(testCase.dRes.spikesRaw, 3), testCase.nSpikes);
            testCase.assertEqual(size(testCase.dRes.spikesFilt, 3), testCase.nSpikes);
            testCase.assertEqual(size(testCase.dRes.spikeFeatures, 3), testCase.nSpikes);
        end

        function partitionOkay(testCase)
            testCase.assertEqual(sum(cellfun(@(c) numel(c), testCase.dRes.spikesBySite)), testCase.nSpikes);
            testCase.assertEqual(sum(cellfun(@(c) numel(c), testCase.dRes.spikesBySite2)), testCase.nSpikes);

            nSites = numel(testCase.dRes.spikesBySite);
            for i = 1:nSites
                iSpikes = testCase.dRes.spikesBySite{i};
                iSpikes2 = testCase.dRes.spikesBySite2{i};

                % no spike should have its secondary site equal to its primary site
                testCase.assertEmpty(intersect(iSpikes, iSpikes2));
                for j = i+1:nSites
                    jSpikes = testCase.dRes.spikesBySite{j};
                    jSpikes2 = testCase.dRes.spikesBySite2{j};

                    testCase.assertEmpty(intersect(iSpikes, jSpikes));
                    testCase.assertEmpty(intersect(iSpikes2, jSpikes2));
                end
            end
        end

        function detIndependence(testCase)
            % detect spikes at the same times in the first file
            firstHalfSame = all(testCase.dRes.spikeTimes(1:testCase.nSpikes0) == testCase.dRes0.spikeTimes);
            testCase.assertTrue(firstHalfSame);

            % detect spikes at the same times in the second file
            secondHalfOffset = unique(testCase.dRes.spikeTimes(end-testCase.nSpikes1+1:end) - testCase.dRes1.spikeTimes);
            testCase.assertEqual(numel(secondHalfOffset), 1);

            % detect the same amplitudes in each half
            firstHalfSame = all(testCase.dRes.spikeAmps(1:testCase.nSpikes0) == testCase.dRes0.spikeAmps);
            secondHalfSame = all(testCase.dRes.spikeAmps(end-testCase.nSpikes1+1:end) == testCase.dRes1.spikeAmps);
            testCase.assertTrue(firstHalfSame);
            testCase.assertTrue(secondHalfSame);

            % detect the same primary/secondary sites in each half
            firstHalfSame = all(testCase.dRes.spikeSites(1:testCase.nSpikes0) == testCase.dRes0.spikeSites);
            secondHalfSame = all(testCase.dRes.spikeSites(end-testCase.nSpikes1+1:end) == testCase.dRes1.spikeSites);
            testCase.assertTrue(firstHalfSame);
            testCase.assertTrue(secondHalfSame);

            firstHalfSame = all(testCase.dRes.spikeSites2(1:testCase.nSpikes0) == testCase.dRes0.spikeSites2);
            secondHalfSame = all(testCase.dRes.spikeSites2(end-testCase.nSpikes1+1:end) == testCase.dRes1.spikeSites2);
            testCase.assertTrue(firstHalfSame);
            testCase.assertTrue(secondHalfSame);

            % compute the same spike positions
            firstHalfNorm = norm(testCase.dRes.spikePositions(1:testCase.nSpikes0, :) - testCase.dRes0.spikePositions);
            secondHalfNorm = norm(testCase.dRes.spikePositions(end-testCase.nSpikes1+1:end, :) - testCase.dRes1.spikePositions);
            testCase.assertEqual(firstHalfNorm, single(0));
            testCase.assertEqual(secondHalfNorm, single(0));

            % extract the same raw/filtered spikes, features
            firstHalfSame = all(all(all(testCase.dRes.spikesRaw(:, :, 1:testCase.nSpikes0) - testCase.dRes0.spikesRaw == 0)));
            secondHalfSame = all(all(all(testCase.dRes.spikesRaw(:, :, end-testCase.nSpikes1+1:end) - testCase.dRes1.spikesRaw == 0)));
            testCase.assertTrue(firstHalfSame);
            testCase.assertTrue(secondHalfSame);

            firstHalfSame = all(all(all(testCase.dRes.spikesFilt(:, :, 1:testCase.nSpikes0) - testCase.dRes0.spikesFilt == 0)));
            secondHalfSame = all(all(all(testCase.dRes.spikesFilt(:, :, end-testCase.nSpikes1+1:end) - testCase.dRes1.spikesFilt == 0)));
            testCase.assertTrue(firstHalfSame);
            testCase.assertTrue(secondHalfSame);

            firstHalfSame = all(all(all(testCase.dRes.spikeFeatures(:, :, 1:testCase.nSpikes0) - testCase.dRes0.spikeFeatures < 1e-8)));
            secondHalfSame = all(all(all(testCase.dRes.spikeFeatures(:, :, end-testCase.nSpikes1+1:end) - testCase.dRes1.spikeFeatures < 1e-8)));
            testCase.assertTrue(firstHalfSame);
            testCase.assertTrue(secondHalfSame);
        end
    end
end

