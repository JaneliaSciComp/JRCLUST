classdef Clustering < handle
    %CLUSTERING Model representing clustering of spike data

    %% OLD-STYLE PROPERTIES, publicly gettable (will be deprecated after a grace period)
    properties (SetObservable, Dependent, Hidden, Transient)
        cviSpk_clu;         % => spikesByCluster
        vnSpk_clu;          % => clusterCounts
    end

    %% NEW-STYLE PROPERTIES
    properties (Dependent, Transient)
        nClusters;          % number of clusters
        clusterCounts;      % number of spikes per cluster
        clusterSites;       % site on which spikes in this cluster most often occur
    end

    properties (SetAccess=private)
        initialClustering;  % initial assignment of spikes to cluster
        spikeClusters;      % individual spike assignments
        spikesByCluster;    % cell array of spike indices per cluster
    end

    %% LIFECYCLE
    methods
        function obj = Clustering(spikeClusters)
            obj.initialClustering = int32(spikeClusters);
            obj.spikeClusters = obj.initialClustering;
            %cviSpk_clu
            obj.spikesByCluster = arrayfun(@(iC) find(obj.spikeClusters == iC), 1:nClusters, 'UniformOutput', 0);
            %vnSpk_clu
            %viSite_clu
            obj.clusterSites = double(arrayfun(@(iC) mode(spikeSites(obj.spikesByCluster{iC})), 1:nClusters));
        end
    end

    %% GETTERS/SETTERS
    methods
        % clusterCounts/vnSpk_clu
        function cc = get.clusterCounts(obj)
            cc = cellfun(@numel, obj.spikesByCluster);
        end
        function cc = get.vnSpk_clu(obj)
            cc = obj.clusterCounts;
        end

        % clusterSites/viSite_clu
        function cs = get.clusterSites(obj)
            cs = double(arrayfun(@(iC) mode(spikeSites(obj.spikesByCluster{iC})), 1:obj.nClusters));
        end

        % nClusters
        function nc = get.nClusters(obj)
            nc = max(double(obj.spikeClusters));
        end
    end
end

