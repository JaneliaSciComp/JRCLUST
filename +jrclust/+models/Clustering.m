classdef Clustering < handle
    %CLUSTERING Model representing clustering of spike data

    %% OLD-STYLE PROPERTIES, publicly gettable (will be deprecated after a grace period)
    properties (SetObservable, Dependent, Hidden, Transient)
        cviSpk_clu;         % => spikesByCluster
        viSite_clu;         % => clusterSites
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

            obj.refresh(true);
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        function refresh(obj, doRemoveEmpty)
            %REFRESH Recount and store spikes by cluster, optionally removing empty clusters
            obj.spikesByCluster = arrayfun(@(iC) find(obj.spikeClusters == iC), 1:obj.nClusters, 'UniformOutput', false);

            if doRemoveEmpty
                obj.removeEmptyClusters();
            end
        end

        function removeEmptyClusters(obj)
            %REMOVEEMPTYCLUSTERS Find and remove empty clusters
            keepClusters = obj.clusterCounts > 0;
            if all(keepClusters)
                return;
            end

            % subset all fields indexed by cluster
            % obj.subsetFields(obj, keepClusters);

            if min(obj.spikeClusters) < 1 % noise or garbage clusters
                obj.spikeClusters(obj.spikeClusters < 1) = 0;

                % renumber clusters 1:numel(unique(spikeClusters))
                [~, ~, obj.spikeClusters] = unique(obj.spikeClusters + 1);
                obj.spikeClusters = obj.spikeClusters - 1;
            else
                % renumber clusters 1:numel(unique(spikeClusters))
                [~, ~, obj.spikeClusters] = unique(obj.spikeClusters);
            end

            obj.spikeClusters = int32(obj.spikeClusters);
        end

        function subsetFields(obj, keepClusters) % WIP
            %SELECT Subset all data fields, taking only those we're interested in
            % excl vnSpk_clu, viSite_clu, vrPosX_clu, vrPosY_clu

            csNames = fieldnames(obj);
            if isempty(csNames)
                return;
            end

            % subset vector fields
            viMatch_v = cellfun(@(vi) ~isempty(vi), cellfun(@(cs)regexp(cs, '^v\w*_clu$'), csNames, 'UniformOutput', false));
            obj = struct_select_(obj, csNames(viMatch_v), keepClusters);

            % subset tensor fields
            viMatch_t = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^t\w*_clu$'), csNames, 'UniformOutput', false));
            obj = struct_select_(obj, csNames(viMatch_t), keepClusters, 3);

            % subset cell fields
            viMatch_c = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^c\w*_clu$'), csNames, 'UniformOutput', false));
            obj = struct_select_(obj, csNames(viMatch_c), keepClusters);

            % remap mrWavCor
            if isprop(obj, 'mrWavCor')
        %         obj.mrWavCor = S_clu_wavcor_remap_(obj, viKeep_clu);
                obj.mrWavCor = obj.mrWavCor(keepClusters, keepClusters);
            end

            % remap mrSim_clu
            if isprop(obj, 'mrSim_clu')
                obj.mrSim_clu = obj.mrSim_clu(keepClusters, keepClusters);
            end
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
        function cs = get.viSite_clu(obj)
            cs = obj.clusterSites;
        end

        % nClusters
        function nc = get.nClusters(obj)
            nc = double(max(obj.spikeClusters));
        end
    end
end
