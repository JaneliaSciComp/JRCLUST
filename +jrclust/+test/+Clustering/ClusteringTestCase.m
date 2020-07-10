classdef ClusteringTestCase < matlab.mock.TestCase
    %CLUSTERINGTESTCASE Superclass of tests for Clustering objects.

    %% FIRST-CLASS PROPS
    properties
        dRes = struct();        % detect results
        sRes = struct();        % sort results
        hCfg;                   % mock jrclust.Config object
        hCfgBehavior;           % behavior object for mock Config
        hClust;                 % clustering object
        histFile = tempname();  % temporary history file
        nSpikes = 131072;       % 2^17 spikes
        nSites = 64;            % 64 sites
    end

    %% DEPENDENT PROPS
    properties (Dependent)
        nClusters;      % number of clusters (one per site)
        spikeAmps;      % spike amplitudes
        spikeClusters;  % spike cluster labels
        spikeFeatures;  % spike clustering features
        spikeTimes;     % spike times
    end

    %% SETUP METHODS
    methods (TestClassSetup)
        function setupConfig(obj)
            %SETUPCONFIG Create a mock jrclust.Config object with just the
            % necessary properties.
            import matlab.mock.actions.AssignOutputs;

            siteNeighbors = zeros(5, obj.nSites);
            for i = 1:obj.nSites
                siteNeighbors(:, i) = mod((i-1:i+3)', obj.nSites) + 1;
            end

            [obj.hCfg, obj.hCfgBehavior] = obj.createMock( ...
                'AddedProperties', ["bitScaling", ...
                                    "filterType", ...
                                    "histFile", ...
                                    "nSitesEvt", ...
                                    "nSitesExcl", ...
                                    "qqFactor", ...
                                    "sampleRate", ...
                                    "siteLoc", ...
                                    "siteNeighbors"], ...
                'AddedMethods', ["isa", ...
                                 "updateLog"] ...
                );
            obj.assignOutputsWhen(obj.hCfgBehavior.isa('jrclust.Config'), true)
            
            obj.assignOutputsWhen(get(obj.hCfgBehavior.bitScaling), 1/pi);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.histFile), 'ndiff');
            obj.assignOutputsWhen(get(obj.hCfgBehavior.histFile), obj.histFile);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.nSitesEvt), 4);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.nSitesExcl), 1);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.qqFactor), 5);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.sampleRate), 25000);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.siteLoc), rand(obj.nSites, 2));
            obj.assignOutputsWhen(get(obj.hCfgBehavior.siteNeighbors), siteNeighbors);
        end

        function setupProps(obj)
            %SETUPPROPS Create the necessary data for a mock clustering.
            rng('default'); rng(10191);
            obj.spikeAmps = randi([-256, 255], obj.nSpikes, 1);
            obj.spikeFeatures = rand(4, 1, obj.nSpikes);
            obj.spikeTimes = (1:obj.nSpikes)';
            obj.spikeClusters = repmat((1:obj.nClusters)', obj.nSpikes/obj.nClusters, 1);
        end
    end

    %% TEARDOWN METHODS
    methods (TestClassTeardown)
        function rmHistFile(obj)
            fclose all;
            delete(obj.histFile);
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

        % spikeTimes
        function st = get.spikeTimes(obj)
            st = obj.dRes.spikeTimes;
        end
        function set.spikeTimes(obj, st)
            obj.dRes.spikeTimes = st;
        end
    end
end

