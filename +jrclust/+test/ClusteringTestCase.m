classdef (Abstract) ClusteringTestCase < jrclust.test.MockConfigTestCase
    %CLUSTERINGTESTCASE Superclass of tests for all subclasses of Clusterin
    %Clustering, or objects with a Clustering member.

    %% FIRST-CLASS PROPS
    properties
        dRes = struct();                % detect results
        sRes = struct();                % sort results
        hClust;                         % Clustering object
    end

    %% DEPENDENT PROPS
    properties (Dependent)
        nClusters;      % number of clusters (one per site)
        spikeAmps;      % spike amplitudes
        spikeClusters;  % spike cluster labels
        spikeFeatures;  % spike clustering features
        spikeSites;     % detected sites for each spike
        spikesBySite;   % indices in spike table of spikes, grouped by site
        spikeTimes;     % spike times
    end

    %% SETUP METHODS
    methods (TestClassSetup)
        function setupProps(obj)
            %SETUPPROPS Create the necessary data for a mock Clustering.
            setupProps@jrclust.test.MockConfigTestCase(obj);

            rng('default'); rng(10191);
            obj.spikeAmps = randi([-256, 255], obj.nSpikes, 1);
            obj.spikeFeatures = rand(obj.hCfg.nSitesEvt, 1, obj.nSpikes);
            obj.spikeSites = repmat((1:obj.nSites)', obj.nSpikes/obj.nSites, 1);
            obj.spikesBySite = arrayfun(@(iS) find(obj.dRes.spikeSites == iS), 1:obj.nSites, 'UniformOutput', 0);
            obj.spikeTimes = (1:obj.nSpikes)';

            obj.spikeClusters = repmat((1:obj.nClusters)', obj.nSpikes/obj.nClusters, 1);
        end
    end

    %% TEARDOWN METHODS
    methods (TestMethodTeardown)
        function reset(obj)
            %RESET Restore the Clustering to its initial state.
            rng('default'); rng(10191);
            obj.hClust.spikeClusters = obj.spikeClusters;
            obj.hClust.clusterNotes = arrayfun(@(i) num2str(i), (1:obj.nClusters)', 'UniformOutput', 0);
            obj.hClust.clusterSites = (1:obj.nClusters)';
            obj.hClust.spikesByCluster = arrayfun(@(i) find(obj.spikeClusters == i), (1:obj.nClusters)', 'UniformOutput', 0);
            obj.hClust.spikesRaw = randi([-1024, 1024], diff(obj.hCfg.evtWindowRawSamp) + 1, 2*obj.hCfg.nSiteDir + 1, obj.nSpikes, 'int16');
            obj.hClust.spikesFilt = randi([-1024, 1024], diff(obj.hCfg.evtWindowSamp) + 1, 2*obj.hCfg.nSiteDir + 1, obj.nSpikes, 'int16');
            obj.hClust.computeQualityScores([]);

            obj.hClust.history = struct('optype', cell(1), 'message', cell(1), 'indices', cell(1));
            obj.hClust.recompute = [];
        end
    end

    %% GETTERS/SETTERS
    methods
        % nClusters
        function nc = get.nClusters(obj)
            nc = obj.nSites;
        end

        % spikeAmps
        function sa = get.spikeAmps(obj)
            sa = obj.dRes.spikeAmps;
        end
        function set.spikeAmps(obj, sa)
           obj.dRes.spikeAmps = sa; 
        end

        % spikeClusters
        function sc = get.spikeClusters(obj)
            sc = obj.sRes.spikeClusters;
        end
        function set.spikeClusters(obj, sc)
            obj.sRes.spikeClusters = sc;
        end

        % spikeFeatures
        function sf = get.spikeFeatures(obj)
            sf = obj.dRes.spikeFeatures;
        end
        function set.spikeFeatures(obj, sf)
            obj.dRes.spikeFeatures = sf;
        end

        % spikeSites
        function ss = get.spikeSites(obj)
            ss = obj.dRes.spikeSites;
        end

        function set.spikeSites(obj, ss)
            obj.dRes.spikeSites = ss;
        end

        % spikesBySite
        function ss = get.spikesBySite(obj)
            ss = obj.dRes.spikesBySite;
        end

        function set.spikesBySite(obj, ss)
            obj.dRes.spikesBySite = ss;
        end

        % spikeTimes
        function st = get.spikeTimes(obj)
            st = obj.dRes.spikeTimes;
        end
        function set.spikeTimes(obj, st)
            obj.dRes.spikeTimes = st;
        end
    end
end

