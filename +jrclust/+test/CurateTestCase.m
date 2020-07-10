classdef CurateTestCase < matlab.mock.TestCase
    %CURATETESTCASE Superclass of tests for CurateController objects.

    %% FIRST-CLASS PROPS
    properties
        dRes = struct();        % detect results
        sRes = struct();        % sort results
        hCfg;                   % mock jrclust.Config object
        hCfgBehavior;           % behavior object for mock Config
        histFile = tempname();  % temporary history file
        hCurate;                % CurateController
        nSpikes = 131072;       % 2^17 spikes
        nSites = 64;            % 64 sites
    end

    %% DEPENDENT PROPS
    properties (Dependent)
        hClust;         % clustering object
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

            % touch histFile
            fclose(fopen(obj.histFile, 'w'));
        end

        function setupProps(obj)
            %SETUPPROPS Create the necessary data for a mock clustering.
            rng('default'); rng(10191);
            obj.spikeFeatures = rand(4, 1, obj.nSpikes);
            obj.spikeTimes = (1:obj.nSpikes)';
            obj.spikeClusters = repmat((1:obj.nSites)', 2^11, 1);
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

