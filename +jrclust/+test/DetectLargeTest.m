classdef DetectLargeTest < matlab.unittest.TestCase
    %DETECTLARGETEST

    properties
        dRes;
        siteLoc;
        testConfig;
        nSpikes = 1433802;
    end

    methods (TestClassSetup)
        function doDetect(testCase)
            testCase.testConfig = fullfile(getenv('JRCTESTDATA'), 'large', 'test.prm');
            hJRC = jrc('detect', testCase.testConfig);

            testCase.dRes = hJRC.dRes;
            testCase.siteLoc = hJRC.hCfg.siteLoc;
        end
    end

    methods (Test)
        function isSorted(testCase)
            testCase.assertTrue(issorted(testCase.dRes.spikeTimes));
        end

        function consistentAcrossDetects(testCase)
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

        function sameTimes(testCase) % this fails pretty hard
            timesAreEqual = find(diff(testCase.dRes.spikeTimes) == 0);
            for i = 1:numel(timesAreEqual)
                t = timesAreEqual(i); % spikeTimes(t) == spikeTimes(t+1)
                sites1 = [testCase.dRes.spikeSites(t) testCase.dRes.spikeSites2(t)];
                sites2 = [testCase.dRes.spikeSites(t+1) testCase.dRes.spikeSites2(t+1)];
                testCase.assertEmpty(intersect(sites1, sites2));
            end
        end
    end
end

