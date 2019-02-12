classdef DetectTest < matlab.unittest.TestCase
    %DETECTTEST Test spike detection
    properties
        hJRC;
        nSpikes = 9891;
    end

    properties (Dependent)
        res;
        hCfg;
    end

    methods (TestClassSetup)
        function setupProps(obj)
            hCfg_ = jrclust.Config(fullfile(getenv('JRCTESTDATA'), 'single', 'test.prm'));
            obj.hJRC = jrclust.JRC(hCfg_);
            obj.hCfg.testRun = 1;
        end

        function doDetect(obj)
            obj.hJRC.detect();
        end
    end

    methods (Test)
        function isSorted(obj)
            %ISSORTED Assert all spike times are in order
            obj.assertTrue(issorted(obj.res.spikeTimes));
        end

        function consistentAcrossDetects(obj)
            %CONSISTENTACROSSDETECTS Assert spike count stays the same from
            %detect to detect
            obj.assertEqual(numel(obj.res.spikeTimes), obj.nSpikes);
        end

        function allCountsEqual(obj)
            %ALLCOUNTSEQUAL Assert spikeAmps, spikeSites, spikeSites2, etc.
            %all have the same number of elements
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
            %PARTITIONOKAY Assert spikesBySite is a partition of spikes
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

        function testMerge(obj)
            %TESTMERGE Assert that no two spiking events are too close
            %together
            timesNearby = find(abs(diff(obj.res.spikeTimes) < obj.hCfg.refracIntSamp));
            for i = 1:numel(timesNearby)
                t = timesNearby(i); % spikeTimes(t) == spikeTimes(t+1)
                sites1 = [obj.res.spikeSites(t) obj.res.spikeSites2(t)];
                sites2 = [obj.res.spikeSites(t+1) obj.res.spikeSites2(t+1)];
                siteIntersection = intersect(sites1, sites2);

                % get site locations for event 1 and event 2 and ensure
                % they're not too close together
                if ~isempty(siteIntersection)
                    site1XY = obj.hCfg.siteLoc(sites1(1), :);
                    site2XY = obj.hCfg.siteLoc(sites2(1), :);
                    obj.assertGreaterThan(pdist2(site1XY, site2XY), obj.hCfg.evtDetectRad);
                end
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

        function res = get.res(obj)
            res = obj.hJRC.res;
        end
    end
end


