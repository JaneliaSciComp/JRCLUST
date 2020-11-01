classdef (Abstract) TemplateClusteringTestCase < jrclust.test.ClusteringTestCase
    %TEMPLATECLUSTERINGTESTCASE Superclass of tests for TemplateClustering,
    %or objects with a TemplateClustering member.

    %% DEPENDENT PROPS
    properties (Dependent)
        amplitudes;         % extracted amplitude for each spike
        templateFeatures;   % matrix giving the magnitude of the projection of each spike onto nTempFeatures other features
        templateFeatureInd; % for each template, the channels with largest amplitudes are indexed in order
        pcFeatures;         % projections of each detected spike onto the top templates most similar to the spike's assigned template
        pcFeatureInd;       % for each template, the other templates with largest similarity are indexed in order
        simScore;           % template-template correlation matrix
        spikeTemplates;     % template assignments for each spike
    end

    %% SETUP METHODS
    methods (TestClassSetup)
        function setupProps(obj)
            setupProps@jrclust.test.ClusteringTestCase(obj);

            % set some fields
            obj.amplitudes = rand(obj.nSpikes, 1);
            obj.spikeTemplates = obj.spikeClusters;
            obj.templateFeatures = rand(3, obj.nSpikes);
            obj.templateFeatureInd = randi(obj.nSites, 3, obj.nSpikes);
            obj.pcFeatures = [];
            obj.pcFeatureInd = [];

            % make a new TemplateClustering
            obj.hClust = jrclust.sort.TemplateClustering(obj.hCfg, obj.sRes, obj.dRes);

            % fill in other appropriate fields
            obj.reset();
            obj.hClust.doRecompute();
        end
    end

    %% TEARDOWN METHODS
    methods (TestMethodTeardown)
        function reset(obj)
            reset@jrclust.test.ClusteringTestCase(obj);

            obj.hClust.simScore = rand(obj.nClusters);
            obj.hClust.spikeTemplates = obj.spikeTemplates;

            obj.hClust.templatesByCluster = arrayfun(@(iC) unique(obj.spikeTemplates(obj.spikeClusters == iC)), ...
                1:obj.nClusters, 'UniformOutput', 0);
            obj.hClust.templateSim = obj.hClust.waveformSim;
        end
    end

    %% GETTERS/SETTERS
    methods
        % amplitudes
        function a = get.amplitudes(obj)
            a = obj.sRes.amplitudes;
        end
        
        function set.amplitudes(obj, a)
            obj.sRes.amplitudes = a;
        end

        % templateFeatures
        function tf = get.templateFeatures(obj)
            tf = obj.sRes.templateFeatures;
        end

        function set.templateFeatures(obj, tf)
            obj.sRes.templateFeatures = tf;
        end

        % templateFeatureInd
        function tfi = get.templateFeatureInd(obj)
            tfi = obj.sRes.templateFeatureInd;
        end

        function set.templateFeatureInd(obj, tfi)
            obj.sRes.templateFeatureInd = tfi;
        end

        % pcFeatures
        function pcf = get.pcFeatures(obj)
            pcf = obj.sRes.pcFeatures;
        end

        function set.pcFeatures(obj, pcf)
            obj.sRes.pcFeatures = pcf;
        end

        % pcFeatureInd
        function pcfi = get.pcFeatureInd(obj)
            pcfi = obj.sRes.pcFeatureInd;
        end

        function set.pcFeatureInd(obj, pcfi)
            obj.sRes.pcFeatureInd = pcfi;
        end

        % simScore
        function ss = get.simScore(obj)
            ss = obj.sRes.simScore;
        end

        function set.simScore(obj, ss)
            obj.sRes.simScore = ss;
        end

        % spikeTemplates
        function st = get.spikeTemplates(obj)
            st = obj.sRes.spikeTemplates;
        end

        function set.spikeTemplates(obj, st)
            obj.sRes.spikeTemplates = st;
        end
    end
end

