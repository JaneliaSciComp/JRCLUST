classdef ClusteringTestCase < matlab.mock.TestCase
    %CLUSTERINGTESTCASE Superclass of tests for Clustering objects.

    %% FIRST-CLASS PROPS
    properties
        dRes = struct();                % detect results
        sRes = struct();                % sort results
        hCfg;                           % mock jrclust.Config object
        hCfgBehavior;                   % behavior object for mock Config
        hClust;                         % clustering object
        histFile = tempname();          % temporary history file
        resFile = [tempname() '.mat'];  % temporary res file
        nSpikes = 8192;                 % 2^13 spikes
        nSites = 64;                    % 64 sites
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
            import matlab.mock.actions.Invoke;

            params = jrclust.utils.getDefaultParams();
            defaultParamNames = fieldnames(params)';

            paramNames = [defaultParamNames ...
                {'histFile', 'resFile'}...
                {'evtWindowSamp', 'evtWindowRawSamp', 'nSites', 'nSitesEvt', 'siteNeighbors'}];

            [obj.hCfg, obj.hCfgBehavior] = obj.createMock( ...
                'AddedProperties', paramNames, ...
                'AddedMethods', ["getOr", ...
                                 "isa", ...
                                 "updateLog", ...
                                ] ...
                );

            % set default param values
            for i = 1:numel(defaultParamNames)
                paramName = defaultParamNames{i};
                param = params.(paramName);

                obj.assignOutputsWhen(get(obj.hCfgBehavior.(paramName)), param.default_value);
            end

            % set unavoidably user-defined values
            nSiteDir = 7;
            nSitesExcl = 2;
            
            obj.assignOutputsWhen(get(obj.hCfgBehavior.histFile), obj.histFile);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.nSiteDir), nSiteDir);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.nSitesExcl), nSitesExcl);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.resFile), obj.resFile);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.siteMap), (1:obj.nSites)');
            obj.assignOutputsWhen(get(obj.hCfgBehavior.siteLoc), rand(obj.nSites, 2));

            % set derived param values
            obj.assignOutputsWhen(get(obj.hCfgBehavior.evtWindowSamp), ...
                round(params.evtWindow.default_value * params.sampleRate.default_value / 1000));
            obj.assignOutputsWhen(get(obj.hCfgBehavior.evtWindowRawSamp), ...
                round(params.evtWindowRaw.default_value * params.sampleRate.default_value / 1000));
            obj.assignOutputsWhen(get(obj.hCfgBehavior.nSites), obj.nSites);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.nSitesEvt), ...
                1 + 2*nSiteDir - nSitesExcl);

            siteNeighbors = zeros(1 + 2*nSiteDir, obj.nSites);
            for i = 1:obj.nSites
                siteNeighbors(:, i) = mod((i-1:i + 2*nSiteDir - 1)', obj.nSites) + 1;
            end
            obj.assignOutputsWhen(get(obj.hCfgBehavior.siteNeighbors), siteNeighbors);

            % set method return values
            when(withAnyInputs(obj.hCfgBehavior.getOr), Invoke(@(varargin) obj.getOr(varargin)));
            obj.assignOutputsWhen(obj.hCfgBehavior.isa('jrclust.Config'), true); % for isa checking
        end

        function setupProps(obj)
            %SETUPPROPS Create the necessary data for a mock clustering.
            rng('default'); rng(10191);
            obj.spikeAmps = randi([-256, 255], obj.nSpikes, 1);
            obj.spikeFeatures = rand(obj.hCfg.nSitesEvt, 1, obj.nSpikes);
            obj.spikeTimes = (1:obj.nSpikes)';
            obj.spikeClusters = repmat((1:obj.nClusters)', obj.nSpikes/obj.nClusters, 1);

            % touch histFile
            fclose(fopen(obj.histFile, 'w'));

            % touch resFile
            fclose(fopen(obj.resFile, 'w'));
        end
    end

    %% TEARDOWN METHODS
    methods (TestClassTeardown)
        function rmHistFile(obj)
            fclose all;
            delete(obj.histFile);
        end
    end

    methods (TestMethodTeardown)
        function resetClustering(obj)
            %RESETCLUSTERING Restore the clustering to its initial state.
            obj.hClust.spikeClusters = obj.spikeClusters;
            obj.hClust.clusterNotes = arrayfun(@(i) num2str(i), (1:obj.nClusters)', 'UniformOutput', 0);
            obj.hClust.clusterSites = (1:obj.nClusters)';
            obj.hClust.spikesByCluster = arrayfun(@(i) find(obj.spikeClusters == i), (1:obj.nClusters)', 'UniformOutput', 0);

            obj.hClust.history = struct('optype', cell(1), 'message', cell(1), 'indices', cell(1));

            obj.hClust.recompute = [];
        end
    end
    
    %% UTILITY METHODS
    methods (Access = protected)
        function val = getOr(varargin)
            val = [];
            args = varargin{2};
            if numel(args) == 3
                val = args{3};
            end
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

