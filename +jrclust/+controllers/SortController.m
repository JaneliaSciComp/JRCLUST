classdef SortController < handle
    %SORTCONTROLLER Handle for sorting spikes into clusters using Rodriguez-Laio

    properties (SetAccess=private, SetObservable, Transient)
        rhoCK;          % CUDA kernel for rho computation
        deltaCK;        % CUDA kernel for delta computation

        chunkSize;
        fTwoStep;
        nC_max;         % maximum number of dimensions per event
        nSubsample;     % number of spikes to sample in estimating cutoff distance

        spikeFeatures;

        hCfg;           % Config object
        errMsg;         % error message, if any
        isError;        % flag indicating an error occurred in sorting
    end

    % LIFECYCLE
    methods
        function obj = SortController(hCfg)
            %SORTCONTROLLER Constructor
            obj.hCfg = hCfg;

            % set default clustering params
            obj.chunkSize = 16;
            obj.fTwoStep = false;
            obj.nC_max = 45;
            obj.nSubsample = 1000;

            obj.isError = false;
        end
    end

    % USER METHODS
    methods
        function res = sort(obj, dRes)
            %SORT Cluster the spikes given in dRes
            res = struct();
            t0 = tic();

            if isfield(dRes, 'spikeFeatures')
                obj.spikeFeatures = dRes.spikeFeatures;
            else
                obj.errMsg = 'cannot sort without features';
                obj.isError = true;
                return;
            end

            nSpikes = numel(dRes.spikeTimes);

            res.spikeRho = zeros(nSpikes, 1, 'single');
            res.spikeDelta = zeros(nSpikes, 1, 'single');
            res.spikeNeigh = zeros(nSpikes, 1, 'uint32');

            % if obj.hCfg.getOr('fDenoise_fet', 0) % not in default prm; deprecated?
            %     trFet_spk = denoise_fet_(trFet_spk, vlRedo_spk);
            % end

            % compute rho
            res = obj.computeRho(dRes, res);

            % compute delta
            res = obj.computeDelta(dRes, res);

            % assign clusters
            [~, res.ordRho] = sort(res.spikeRho, 'descend');
            res = obj.assignClusters(dRes, res);
            
            res.runtime = toc(t0);
        end
    end

    % RODRIGUEZ-LAIO METHODS
    methods (Access=protected, Hidden)
        function res = computeRho(obj, dRes, res)
            %COMPUTERHO Compute rho values for spike features
            nSites = numel(obj.hCfg.siteMap);

            res.rhoCutSite = zeros(nSites, 1, 'single');
            res.rhoCutGlobal = [];

            % compute rho cutoff globally
            if obj.hCfg.getOr('fDc_global', true)
                res.rhoCutGlobal = obj.estRhoCutGlobal(dRes);
            end

            if obj.hCfg.verbose
                fprintf('Computing rho\n\t');
                t1 = tic;
            end

            spikeData = struct('spikeTimes', dRes.spikeTimes);
            for iSite = 1:nSites
                if isfield(dRes, 'spikesBySite')
                    spikeData.spikes1 = dRes.spikesBySite{iSite};
                else
                    spikeData.spikes1 = dRes.spikeSites(dRes.spikeSites==iSite);
                end

                if isfield(dRes, 'spikesBySite2') && ~isempty(dRes.spikesBySite2)
                    spikeData.spikes2 = dRes.spikesBySite2{iSite};
                end

                if isfield(dRes, 'spikesBySite3') && ~isempty(dRes.spikesBySite3)
                    spikeData.spikes3 = dRes.spikesBySite3{iSite};
                end

                [siteFeatures, ~, n1, n2, spikeOrder] = jrclust.features.getSiteFeatures(obj.spikeFeatures, iSite, spikeData, obj.hCfg);

                if isempty(siteFeatures)
                    continue;
                end

                siteFeatures = jrclust.utils.tryGpuArray(siteFeatures);
                spikeOrder = jrclust.utils.tryGpuArray(spikeOrder);

                if isempty(res.rhoCutGlobal) % estimate rhoCut in CPU
                    siteCut = obj.estRhoCutSite(siteFeatures, spikeOrder, n1, n2);
                else
                    siteCut = res.rhoCutGlobal.^2;
                end

                siteRho = obj.computeRhoSite(siteFeatures, spikeOrder, n1, n2, siteCut);

%                 if ~isempty(vlRedo_spk), viSpk_site_ = viSpk_site_(vlRedo_spk(viSpk_site_)); end
                res.spikeRho(spikeData.spikes1) = jrclust.utils.tryGather(siteRho);
                res.rhoCutSite(iSite) = jrclust.utils.tryGather(siteCut);
                clear siteFeatures spikeOrder siteRho;
                fprintf('.');
            end

            if obj.hCfg.verbose
                fprintf('\n\ttook %0.1fs\n', toc(t1));
            end
        end

        function rho = computeRhoSite(obj, siteFeatures, spikeOrder, n1, n2, rhoCut)
            %COMPUTERHOSITE Compute site-wise rho for spike features
            [nC, n12] = size(siteFeatures); % nc is constant with the loop
            dn_max = int32(round((n1 + n2) / obj.hCfg.nTime_clu));

            rhoCut = single(rhoCut);

            if obj.hCfg.useGPU
                try
                    if isempty(obj.rhoCK) % create CUDA kernel
                        ptxFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_rho.ptx');
                        cuFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_rho.cu');
                        obj.rhoCK = parallel.gpu.CUDAKernel(ptxFile, cuFile);
                        obj.rhoCK.ThreadBlockSize = [obj.hCfg.nThreads, 1];
                        obj.rhoCK.SharedMemorySize = 4 * obj.chunkSize * (2 + obj.nC_max + 2 * obj.hCfg.nThreads); % @TODO: update the size
                    end

                    obj.rhoCK.GridSize = [ceil(n1 / obj.chunkSize / obj.chunkSize), obj.chunkSize]; %MaxGridSize: [2.1475e+09 65535 65535]
                    rho = zeros(1, n1, 'single', 'gpuArray');
                    consts = int32([n1, n12, nC, dn_max, obj.hCfg.getOr('fDc_spk', 0)]);
                    rho = feval(obj.rhoCK, rho, siteFeatures, spikeOrder, consts, rhoCut);

                    return;
                catch ME
                    disperr_('CUDA kernel failed. Re-trying in CPU.');
                    obj.rhoCK = [];
                end
            end

            % retry in CPU
            [siteFeatures, spikeOrder] = jrclust.utils.tryGather(siteFeatures, spikeOrder);
            spikeOrderN1 = spikeOrder(1:n1)';
            nearby = abs(bsxfun(@minus, spikeOrder, spikeOrderN1)) <= dn_max;
            rho = sum(nearby & pdist2(siteFeatures', siteFeatures(:,1:n1)').^2 < rhoCut) - 1; %do not include self
            rho = single(rho ./ sum(nearby)); % normalize
        end

        function rhoCut = estRhoCutGlobal(obj, dRes, vlRedo_spk)
            %ESTRHOCUTGLOBAL Estimate rho cutoff distance over all sites
            if nargin < 3
                vlRedo_spk = [];
            end

            spikeData = struct('spikeTimes', dRes.spikeTimes, ...
                               'vlRedo_spk', vlRedo_spk);

            if obj.hCfg.verbose
                fprintf('Estimating cutoff distance\n\t');
                t1 = tic;
            end

            nSites = numel(obj.hCfg.siteMap);

            % estimate rho cutoff for each site
            siteCuts = nan(1, nSites);
            for iSite = 1:nSites
                if isfield(dRes, 'spikesBySite')
                    spikeData.spikes1 = dRes.spikesBySite{iSite};
                else
                    spikeData.spikes1 = dRes.spikeSites(dRes.spikeSites==iSite);
                end

                if isfield(dRes, 'spikesBySite2') && ~isempty(dRes.spikesBySite2)
                    spikeData.spikes2 = dRes.spikesBySite2{iSite};
                end

                if isfield(dRes, 'spikesBySite3') && ~isempty(dRes.spikesBySite3)
                    spikeData.spikes3 = dRes.spikesBySite3{iSite};
                end

                [siteFeatures, ~, n1, n2, spikeOrder] = jrclust.features.getSiteFeatures(obj.spikeFeatures, iSite, spikeData, obj.hCfg);

                if isempty(siteFeatures)
                    continue;
                end

                siteCuts(iSite) = obj.estRhoCutSite(siteFeatures, spikeOrder, n1, n2);
            end

            rhoCut = sqrt(abs(quantile(siteCuts, .5)));

            if obj.hCfg.verbose
                fprintf('\n\ttook %0.1fs\n', toc(t1));
            end
        end

        function rhoCut = estRhoCutSite(obj, siteFeatures, spikeOrder, n1, n2)
            %ESTRHOCUTSITE Estimate site-wise rho cutoff
            if obj.hCfg.getOr('fDc_spk', 0)
                rhoCut = (obj.hCfg.dc_percent/100).^2;
                return;
            end

            [nPrimary, nSecondary] = deal(obj.nSubsample, 4*obj.nSubsample);

            ss = jrclust.utils.subsample(1:n1, nPrimary);
            spikeOrderPrimary = spikeOrder(ss);
            featuresPrimary = siteFeatures(:, ss);

            ss = jrclust.utils.subsample(1:n1, nSecondary);
            spikeOrder = spikeOrder(ss);
            siteFeatures = siteFeatures(:, ss);

            for iRetry = 1:2
                try
                    featureDists = pdist2(siteFeatures', featuresPrimary').^2;

                    fSubset = abs(bsxfun(@minus, spikeOrder, spikeOrderPrimary')) < (n1 + n2) / obj.hCfg.nTime_clu;
                    featureDists(~fSubset) = nan;
                catch ME
                    siteFeatures = jrclust.utils.tryGather(siteFeatures);
                end
            end

            if obj.hCfg.getOr('fDc_subsample_mode', false)
                featureDists(featureDists <= 0) = nan;
                rhoCut = quantile(featureDists(~isnan(featureDists)), obj.hCfg.dc_percent/100);
            else
                featureDists = jrclust.utils.tryGather(featureDists);
                featureDists(featureDists <=0) = nan;
                rhoCut = nanmedian(quantile(featureDists, obj.hCfg.dc_percent/100));

                if isnan(rhoCut) % featureDists was completely nan
                    rhoCut = quantile(featureDists(:), obj.hCfg.dc_percent/100);
                end
            end

            fprintf('.');
        end

        function res = computeDelta(obj, dRes, res)
            %COMPUTEDELTA Compute delta for spike features
            nSites = numel(obj.hCfg.siteMap);

            if obj.hCfg.verbose
                fprintf('Computing delta\n\t');
                t2 = tic;
            end

            spikeData = struct('spikeTimes', dRes.spikeTimes);
            for iSite = 1:nSites
                if isfield(dRes, 'spikesBySite')
                    spikeData.spikes1 = dRes.spikesBySite{iSite};
                else
                    spikeData.spikes1 = dRes.spikeSites(dRes.spikeSites==iSite);
                end

                if isfield(dRes, 'spikesBySite2') && ~isempty(dRes.spikesBySite2)
                    spikeData.spikes2 = dRes.spikesBySite2{iSite};
                end

                if isfield(dRes, 'spikesBySite3') && ~isempty(dRes.spikesBySite3)
                    spikeData.spikes3 = dRes.spikesBySite3{iSite};
                end

                [siteFeatures, spikes, n1, n2, spikeOrder] = jrclust.features.getSiteFeatures(obj.spikeFeatures, iSite, spikeData, obj.hCfg);

                if isempty(siteFeatures)
                    continue;
                end

                rhoOrder = jrclust.utils.rankorder(res.spikeRho(spikes), 'descend');
                siteFeatures = jrclust.utils.tryGpuArray(siteFeatures);
                rhoOrder = jrclust.utils.tryGpuArray(rhoOrder);
                spikeOrder = jrclust.utils.tryGpuArray(spikeOrder);

                try
                    [siteDelta, siteNN] = obj.computeDeltaSite(siteFeatures, spikeOrder, rhoOrder, n1, n2, res.rhoCutSite(iSite));
                    [siteDelta, siteNN] = jrclust.utils.tryGather(siteDelta, siteNN);
                catch ME
                    disperr_(sprintf('error at site# %d', iSite), ME);
                end

                % if ~isempty(vlRedo_spk), viSpk_site_ = viSpk_site_(vlRedo_spk(viSpk_site_)); end

                res.spikeDelta(spikeData.spikes1) = siteDelta;
                res.spikeNeigh(spikeData.spikes1) = spikes(siteNN);
                [siteFeatures, rhoOrder, spikeOrder] = jrclust.utils.tryGather(siteFeatures, rhoOrder, spikeOrder);
                clear('siteFeatures', 'rhoOrder', 'spikeOrder');
                fprintf('.');
            end

            % set delta for spikes with no nearest neighbor of higher
            % density to maximum distance to ensure they don't get
            % grouped into another cluster
            nanDelta = find(isnan(res.spikeDelta));
            if ~isempty(nanDelta)
                res.spikeDelta(nanDelta) = max(spikeDelta);
            end

            if obj.hCfg.verbose
                fprintf('\n\ttook %0.2fs\n', toc(t2));
            end
        end

        function [delta, nNeigh] = computeDeltaSite(obj, siteFeatures, spikeOrder, rhoOrder, n1, n2, rhoCut)
            %COMPUTEDELTASITE Compute site-wise delta for spike features
            [nC, n12] = size(siteFeatures); % nc is constant with the loop
            dn_max = int32(round((n1 + n2) / obj.hCfg.nTime_clu));

            if obj.hCfg.fGpu
                try
                    if isempty(obj.deltaCK) % create cuda kernel
                        ptxFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_delta.ptx');
                        cuFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_delta.cu');
                        obj.deltaCK = parallel.gpu.CUDAKernel(ptxFile, cuFile);
                        obj.deltaCK.ThreadBlockSize = [obj.hCfg.nThreads, 1];
                        obj.deltaCK.SharedMemorySize = 4 * obj.chunkSize * (3 + obj.nC_max + 2*obj.hCfg.nThreads); % @TODO: update the size
                    end

                    % set every time function is called
                    obj.deltaCK.GridSize = [ceil(n1 / obj.chunkSize / obj.chunkSize), obj.chunkSize]; %MaxGridSize: [2.1475e+09 65535 65535]
                    delta = zeros([1, n1], 'single', 'gpuArray');
                    nNeigh = zeros([1, n1], 'uint32', 'gpuArray');
                    consts = int32([n1, n12, nC, dn_max, obj.hCfg.getOr('fDc_spk', 0)]);
                    [delta, nNeigh] = feval(obj.deltaCK, delta, nNeigh, siteFeatures, spikeOrder, rhoOrder, consts, rhoCut);

                    return;
                catch ME
                    disperr_('CUDA kernel failed. Re-trying in CPU.');
                    obj.deltaCK = [];
                end
            end

            dists = pdist2(siteFeatures', siteFeatures(:,1:n1)').^2;
            nearby = bsxfun(@lt, rhoOrder, rhoOrder(1:n1)') & abs(bsxfun(@minus, spikeOrder, spikeOrder(1:n1)')) <= dn_max;
            dists(~nearby) = nan;
            [delta, nNeigh] = min(dists);
        end
    end

    % CLUSTER ASSIGNMENT/MERGING METHODS
    methods (Access=protected, Hidden)
        function res = assignClusters(obj, dRes, res) % WIP
            if ~isfield(dRes, 'spikesBySite')
                spikesBySite = arrayfun(@(iSite) dRes.spikes(dRes.spikeSites==iSite), obj.hCfg.siteMap, 'UniformOutput', false);
            end

            if strcmp(obj.hCfg.rlDetrendMode, 'local')      % perform detrending site by site
                res.centers = jrclust.clustering.detrendRhoDelta(res, spikesBySite, true, obj.hCfg);
            elseif strcmp(obj.hCfg.rlDetrendMode, 'global') % detrend over all sites
                res.centers = jrclust.clustering.detrendRhoDelta(res, spikesBySite, false, obj.hCfg);
            elseif strcmp(obj.hCfg.rlDetrendMode, 'logz')   % identify centers with sufficiently high z-scores
                % res.centers = log_ztran_(res.spikeRho, res.spikeDelta, obj.hCfg.log10RhoCut, 4 + obj.hCfg.log10DeltaCut);
                x = log10(res.spikeRho(:));
                y = log10(res.spikeDelta(:));

                mask = isfinite(x) & isfinite(y);
                y(mask) = jrclust.utils.zscore(y(mask));

                % from postCluster_: 4+P.delta1_cut
                res.centers = find(x >= obj.hCfg.log10RhoCut & y >= 4 + obj.hCfg.log10DeltaCut);
            else % strcmp(obj.hCfg.rlDetrendMode, 'none')
                res.centers = find(res.spikeRho(:) > 10^(obj.hCfg.log10RhoCut) & res.spikeDelta(:) > 10^(obj.hCfg.log10DeltaCut));                
            end

            spikeClusters = [];

            % res = assign_clu_count_(res, obj.hCfg); % enforce min count algorithm
            nRepeat_max = 1000;
            if isempty(spikeClusters)
                nClu_pre = [];
            else
                nClu_pre = res.nClu;
            end
            nClu_rm = 0;
            fprintf('assigning clusters, nClu:%d\n', numel(res.centers));
            t1 = tic;

            % assign spikes to clusters
            for iRepeat = 1:nRepeat_max % repeat 1000 times max
                % [spikeClusters, res.centers] = assignCluster_(spikeClusters, res.ordRho, res.nneigh, res.centers);
                nSpikes = numel(res.ordRho);
                nClusters = numel(res.centers);

                if isempty(spikeClusters)
                    spikeClusters = zeros([nSpikes, 1], 'int32');
                    spikeClusters(res.centers) = 1:nClusters;
                end

                % one or no center, assign all spikes to one cluster
                if numel(res.centers) == 0 || numel(res.centers) == 1
                    spikeClusters = ones([nSpikes, 1], 'int32');
                    res.centers = res.ordRho(1);
                else
                    nneigh1 = res.spikeNeigh(res.ordRho);
                    for i = 1:10
                        vi = find(spikeClusters(res.ordRho) <=0);
                        if isempty(vi)
                            break;
                        end

                        vi = vi(:)';

                        for j = vi
                            spikeClusters(res.ordRho(j)) = spikeClusters(nneigh1(j));
                        end
                        nUnassigned = sum(spikeClusters <= 0);

                        if nUnassigned==0
                            break;
                        end

                        fprintf('i:%d, n0=%d, ', i, nUnassigned);
                    end
                    spikeClusters(spikeClusters <= 0) = 1; %background
                end

                obj.hCfg.minClusterSize = max(obj.hCfg.minClusterSize, 2*size(dRes.spikeFeatures, 1));
                % http://scikit-learn.org/stable/modules/lda_qda.html

                hClust = jrclust.models.Clustering(spikeClusters);

                % remove clusters unused
                viCluKill = find(res.vnSpk_clu <= obj.hCfg.minClusterSize);
                if isempty(viCluKill)
                    break;
                end

                res.centers(viCluKill) = [];
                spikeClusters = [];
                nClu_rm = nClu_rm + numel(viCluKill);

                if iRepeat == nRepeat_max
                    fprintf(2, 'assign_clu_count_: exceeded nRepeat_max=%d\n', nRepeat_max);
                end
            end

            res.hClust = hClust;

            fprintf('\n\ttook %0.1fs. Removed %d clusters having <%d spikes: %d->%d\n', ...
            toc(t1), nClu_rm, obj.hCfg.minClusterSize, nClu_pre, res.nClu);
        end
    end
end

