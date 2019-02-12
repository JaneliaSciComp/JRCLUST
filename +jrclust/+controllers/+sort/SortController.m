classdef SortController < handle
    %SORTCONTROLLER Handle for sorting spikes into clusters using Rodriguez-Laio
    properties (Access=private, Transient)
        rhoCK;          % CUDA kernel for rho computation
        deltaCK;        % CUDA kernel for delta computation

        hCfg;           % Config object
    end

    properties(SetAccess=private, Transient)
        errMsg;         % error message, if any
        isError;        % flag indicating an error occurred in sorting
    end

    %% LIFECYCLE
    methods
        function obj = SortController(hCfg)
            %SORTCONTROLLER Constructor
            obj.hCfg = hCfg;
            obj.isError = 0;
        end
    end

    %% USER METHODS
    methods
        function res = sort(obj, dRes)
            %SORT Cluster the spikes given in dRes
            res = struct();
            t0 = tic();

            if ~isfield(dRes, 'spikeFeatures')
                obj.errMsg = 'cannot sort without features';
                obj.isError = 1;
                return;
            end

            nSpikes = numel(dRes.spikeTimes);

            res.spikeRho = zeros(nSpikes, 1, 'single');
            res.spikeDelta = zeros(nSpikes, 1, 'single');
            res.spikeNeigh = zeros(nSpikes, 1, 'uint32');

            % compute rho
            res = jrclust.cluster.densitypeaks.computeRho(dRes, res, obj.hCfg);

            % compute delta
            res = jrclust.cluster.densitypeaks.computeDelta(dRes, res, obj.hCfg);

            % assign clusters
            [~, res.ordRho] = sort(res.spikeRho, 'descend');

            res = jrclust.cluster.densitypeaks.assignClusters(dRes, res, obj.hCfg);
            hClust = jrclust.models.clustering.DensityPeakClustering(res, dRes, obj.hCfg);
            hClust.autoMerge();

            res.hClust = hClust;

            % if get_set_(P, 'fCorrect_overlap', 0) % correct waveforms and features after correcting clusters
            %     S_clu = sort_overlap_(S0, S_clu, P);
            % end

            % summarize
            res.sortTime = toc(t0);
            res.sortedOn = now();
        end
    end
end

