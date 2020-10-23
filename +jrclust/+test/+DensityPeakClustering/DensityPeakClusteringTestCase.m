classdef DensityPeakClusteringTestCase < jrclust.test.Clustering.ClusteringTestCase

    %% DEPENDENT PROPS
    properties (Dependent)
        meanSiteThresh; % mean detection threshold per site
    end

    %% SETUP METHODS
    methods (TestClassSetup)
        function setupProps(obj)
            setupProps@jrclust.test.Clustering.ClusteringTestCase(obj);

            % set some fields
            obj.dRes.meanSiteThresh = rand(obj.nSites, 1);
            obj.dRes.spikeSites = repmat((1:obj.nSites)', obj.nSpikes/obj.nSites, 1);
            obj.dRes.spikesBySite = arrayfun(@(iS) find(obj.dRes.spikeSites == iS), 1:obj.nSites, 'UniformOutput', 0);
            obj.dRes.spikesBySite2 = cell(obj.nSites, 1);

            % make a new DensityPeakClustering
            obj.hClust = jrclust.sort.DensityPeakClustering(obj.hCfg, obj.sRes, obj.dRes);

            % fill in the appropriate fields
            obj.resetClustering();
        end
    end

    %% TEARDOWN METHODS
    methods (TestMethodTeardown)
        function resetClustering(obj)
            %RESETCLUSTERING Restore the clustering to its initial state.
            resetClustering@jrclust.test.Clustering.ClusteringTestCase(obj);

            obj.hClust.computeCentroids([]);
            obj.hClust.computeQualityScores([]);
        end
    end

    %% GETTERS/SETTERS
    methods
        function mst = get.meanSiteThresh(obj)
            mst = obj.dRes.meanSiteThresh;
        end
    end
end

