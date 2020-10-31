classdef (Abstract) DensityPeakClusteringTestCase < jrclust.test.ClusteringTestCase
    %DENSITYPEAKCLUSTERINGTESTCASE Superclass of tests for
    %DensityPeakClustering, or objects with a DensityPeakClustering member.

    %% DEPENDENT PROPS
    properties (Dependent)
        meanSiteThresh; % mean detection threshold per site
        spikesBySite2;  % indices in spike table of spikes, grouped by secondary site
        spikeRho;       % rho values for each spike
        spikeDelta;     % delta values for each spike
    end

    %% SETUP METHODS
    methods (TestClassSetup)
        function setupProps(obj)
            setupProps@jrclust.test.ClusteringTestCase(obj);

            % set some fields
            obj.meanSiteThresh = rand(obj.nSites, 1);
            obj.spikesBySite2 = cell(obj.nSites, 1);
            obj.spikeRho = rand(obj.nSpikes, 1);
            obj.spikeDelta = rand(obj.nSpikes, 1);

            % make a new DensityPeakClustering
            obj.hClust = jrclust.sort.DensityPeakClustering(obj.hCfg, obj.sRes, obj.dRes);

            % fill in the appropriate fields
            obj.reset();
        end
    end

    %% TEARDOWN METHODS
    methods (TestMethodTeardown)
        function reset(obj)
            reset@jrclust.test.ClusteringTestCase(obj);

            obj.hClust.assignClusterCenters();
        end
    end

    %% GETTERS/SETTERS
    methods
        % meanSiteThresh
        function mst = get.meanSiteThresh(obj)
            mst = obj.dRes.meanSiteThresh;
        end
        
        function set.meanSiteThresh(obj, mst)
            obj.dRes.meanSiteThresh = mst;
        end

        % spikesBySite2
        function ss = get.spikesBySite2(obj)
            ss = obj.dRes.spikesBySite2;
        end

        function set.spikesBySite2(obj, ss)
            obj.dRes.spikesBySite2 = ss;
        end

        % spikeRho
        function sr = get.spikeRho(obj)
            sr = obj.sRes.spikeRho;
        end

        function set.spikeRho(obj, sr)
            obj.sRes.spikeRho = sr;
        end

        % spikeDelta
        function sd = get.spikeDelta(obj)
            sd = obj.sRes.spikeDelta;
        end

        function set.spikeDelta(obj, sd)
            obj.sRes.spikeDelta = sd;
        end
    end
end

