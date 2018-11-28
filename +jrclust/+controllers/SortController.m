classdef SortController < handle
    %SORTCONTROLLER Handle for sorting spikes into clusters using Rodriguez-Laio

    properties (SetAccess=private, SetObservable, Transient)
        chunkSize;      
        nC_max;         % maximum number of dimensions per event
        nSubsample;     % number of spikes to sample in estimating cutoff distance
        spikeFeatures;  % features to cluster
    end

    properties (Access=private, Transient)
        rhoCK;          % CUDA kernel for rho computation
        deltaCK;        % CUDA kernel for delta computation

        hCfg;           % Config object
        errMsg;         % error message, if any
        isError;        % flag indicating an error occurred in sorting
    end

    %% LIFECYCLE
    methods
        function obj = SortController(hCfg)
            %SORTCONTROLLER Constructor
            obj.hCfg = hCfg;

            % set default clustering params
            obj.chunkSize = 16;
            obj.nC_max = 45;
            obj.nSubsample = 1000;

            obj.isError = false;
        end
    end

    %% USER METHODS
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

            % compute rho
            res = obj.computeRho(dRes, res);

            % compute delta
            res = obj.computeDelta(dRes, res);

            % assign clusters
            [~, res.ordRho] = sort(res.spikeRho, 'descend');
            res = obj.assignClusters(dRes, res);

%             if obj.hCfg.repeatLower
%                 
%             end

            % summarize
            res.runtime = toc(t0);
            res.completedAt = now();
        end
    end

    %% RODRIGUEZ-LAIO METHODS
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

    %% CLUSTER ASSIGNMENT/MERGING METHODS
    methods (Access=protected, Hidden)
        function res = computeCenters(obj, dRes, res)
            %COMPUTECENTERS Find cluster centers
            if ~isfield(dRes, 'spikesBySite')
                dRes.spikesBySite = arrayfun(@(iSite) dRes.spikes(dRes.spikeSites==iSite), obj.hCfg.siteMap, 'UniformOutput', false);
            end

            if strcmp(obj.hCfg.rlDetrendMode, 'local')      % perform detrending site by site
                res.centers = jrclust.clustering.detrendRhoDelta(res, dRes.spikesBySite, true, obj.hCfg);
            elseif strcmp(obj.hCfg.rlDetrendMode, 'global') % detrend over all sites
                res.centers = jrclust.clustering.detrendRhoDelta(res, dRes.spikesBySite, false, obj.hCfg);
            elseif strcmp(obj.hCfg.rlDetrendMode, 'logz')   % identify centers with sufficiently high z-scores
                % res.centers = log_ztran_(res.spikeRho, res.spikeDelta, obj.hCfg.log10RhoCut, 4 + obj.hCfg.log10DeltaCut);
                x = log10(res.spikeRho(:));
                y = log10(res.spikeDelta(:));

                mask = isfinite(x) & isfinite(y);
                y(mask) = jrclust.utils.zscore(y(mask));

                % from postCluster_: 4+P.delta1_cut
                res.centers = find(x >= obj.hCfg.log10RhoCut & y >= 4 + obj.hCfg.log10DeltaCut);
            else                                            % don't detrend
                res.centers = find(res.spikeRho(:) > 10^(obj.hCfg.log10RhoCut) & res.spikeDelta(:) > 10^(obj.hCfg.log10DeltaCut));                
            end
        end

        function res = assignClusters(obj, dRes, res)
            %ASSIGNCLUSTERS Given rho-delta information, assign spikes to clusters
            res = obj.computeCenters(dRes, res);
            res.spikeClusters = [];

            % res = assign_clu_count_(res, obj.hCfg); % enforce min count algorithm
            nRepeat_max = 1000;
            if isempty(res.spikeClusters)
                nClustersPrev = [];
            else
                nClustersPrev = res.nClusters;
            end

            removedClusters = 0;
            if obj.hCfg.verbose
                fprintf('assigning clusters, nClusters:%d\n', numel(res.centers));
                t = tic;
            end

            % assign spikes to clusters
            for iRepeat = 1:nRepeat_max % repeat 1000 times max
                % [res.spikeClusters, res.centers] = assignCluster_(res.spikeClusters, res.ordRho, res.nneigh, res.centers);
                nSpikes = numel(res.ordRho);
                nClusters = numel(res.centers);

                if isempty(res.spikeClusters)
                    res.spikeClusters = zeros([nSpikes, 1], 'int32');
                    res.spikeClusters(res.centers) = 1:nClusters;
                end

                % one or no center, assign all spikes to one cluster
                if numel(res.centers) == 0 || numel(res.centers) == 1
                    res.spikeClusters = ones([nSpikes, 1], 'int32');
                    res.centers = res.ordRho(1);
                else
                    nNeigh = res.spikeNeigh(res.ordRho);

                    for i = 1:10
                        unassigned = find(res.spikeClusters(res.ordRho) <=0);
                        if isempty(unassigned)
                            break;
                        end

                        unassigned = unassigned(:)';

                        for j = unassigned
                            res.spikeClusters(res.ordRho(j)) = res.spikeClusters(nNeigh(j));
                        end
                        nUnassigned = sum(res.spikeClusters <= 0);

                        if nUnassigned==0
                            break;
                        end

                        fprintf('i:%d, n0=%d, ', i, nUnassigned);
                    end
                    res.spikeClusters(res.spikeClusters <= 0) = 1; %background
                end

                obj.hCfg.minClusterSize = max(obj.hCfg.minClusterSize, 2*size(obj.spikeFeatures, 1));
                % http://scikit-learn.org/stable/modules/lda_qda.html

                % count spikes in clusters
                %cviSpk_clu
                res.spikesByCluster = arrayfun(@(iC) find(res.spikeClusters == iC), 1:nClusters, 'UniformOutput', 0);
                %vnSpk_clu
                res.clusterCounts = cellfun(@numel, res.spikesByCluster);
                %viSite_clu
                res.clusterSites = double(arrayfun(@(iC) mode(dRes.spikeSites(res.spikesByCluster{iC})), 1:nClusters));

                % remove clusters unused
                smallClusters = find(res.clusterCounts <= obj.hCfg.minClusterSize);
                if isempty(smallClusters)
                    break;
                end

                res.centers(smallClusters) = [];
                res.spikeClusters = [];
                removedClusters = removedClusters + numel(smallClusters);

                if iRepeat == nRepeat_max
                    warning('assignClusters: exceeded nRepeat_max = %d\n', nRepeat_max);
                end
            end

            res.hClust = jrclust.models.Clustering(res.spikeClusters, obj.spikeFeatures, dRes.spikeSites);

            if obj.hCfg.verbose
                fprintf('\n\ttook %0.1fs. Removed %d clusters having <%d spikes: %d->%d\n', toc(t), removedClusters, obj.hCfg.minClusterSize, nClustersPrev, res.hClust.nClusters);
            end
        end

%         function res = S_clu_reclust_(obj, dRes, res, S0)
%             vcMode_divide = 'amp'; % {'amp', 'density', 'fet'}
% 
%             trFet_spk0 = obj.spikeFeatures;
%             nSites_fet = obj.hCfg.maxSite*2+1-obj.hCfg.nSites_ref;
%             nFetPerSite = size(trFet_spk,1) / nSites_fet;
% 
%             switch vcMode_divide
% %                 case 'nneigh'
% %                     % nearest neighbor averaging per same ecluster for feature enhancement
% %                     trFet_spk = nneigh_ave_(S_clu, obj.hCfg, trFet_spk);
% %                     P1 = setfield(obj.hCfg, 'nRepeat_fet', 1);
% %                     S_clu = postCluster_(cluster_spacetime_(S0, P1), obj.hCfg);
% %                     trFet_spk = trFet_spk0; % undo fet change (undo fet cleanup)
% % 
% %                 case 'fet'
% %                     % recompute pca and
% %                     vrSnr_clu = S_clu_snr_(S_clu);
% %                     vlRedo_clu = vrSnr_clu < quantile(vrSnr_clu, 1/nFetPerSite);
% %                     vlRedo_spk = ismember(S_clu.viClu, find(vlRedo_clu));
% %                     tnWav_spk = get_spkwav_(obj.hCfg, 0);
% %                     trWav2_spk = single(permute(tnWav_spk(:,:,vlRedo_spk), [1,3,2]));
% %                     trWav2_spk = spkwav_car_(trWav2_spk, obj.hCfg);
% %                     [mrPv, vrD1] = tnWav2pv_(trWav2_spk, obj.hCfg);
% %                     dimm1 = size(trWav2_spk);
% %                     mrWav_spk1 = reshape(trWav2_spk, dimm1(1), []);
% %                     trFet_spk_ = reshape(mrPv(:,1)' * mrWav_spk1, dimm1(2:3))';
% 
%                 case 'amp'
%                     vrSnr_clu = S_clu_snr_(res);
%                     try
%                         %         for iRepeat = (nFetPerSite-1):-1:1
%                         vlRedo_clu = vrSnr_clu < quantile(vrSnr_clu, 1/2);
%                         vlRedo_spk = ismember(res.viClu, find(vlRedo_clu));
% 
%                         % reproject the feature
%                         %         nSpk_ = sum(vlRedo_spk);
%                         %         nFets_spk_ = ceil(size(trFet_spk,1)/2);
%                         %         trFet_spk_ = pca(reshape(trFet_spk(:,:,vlRedo_spk), size(trFet_spk,1), []), 'NumComponents', nFets_spk_);
%                         %         trFet_spk_ = permute(reshape(trFet_spk_, [size(trFet_spk,2), nSpk_, nFets_spk_]), [3,1,2]);
%                         %         trFet_spk = trFet_spk(1:nFets_spk_,:,:);
%                         %         trFet_spk(:,:,vlRedo_spk) = trFet_spk_;
% 
%                         mlFet_ = false(nSites_fet, nFetPerSite);
%                         nSites_fet = ceil(nSites_fet*.5); %*.75
%                         mlFet_(1:nSites_fet, 1) = 1;
%                         %             mlFet_(1,:) = 1;
%                         %             mlFet_(:, 1) = 1;
%                         trFet_spk = trFet_spk0(find(mlFet_),:,:);
% 
%                         %             trFet_spk = trFet_spk0(1:1*nSites_fet,:,:);
%                         S_clu_B = postCluster_(cluster_spacetime_(S0, obj.hCfg, vlRedo_spk), obj.hCfg);
%                         res = S_clu_combine_(res, S_clu_B, vlRedo_clu, vlRedo_spk);
%                     catch
%                         disperr_();
%                     end
%                     trFet_spk = trFet_spk0; %restore
% 
%                 case 'density'
%                     vlRedo_clu = res.vnSpk_clu > quantile(res.vnSpk_clu, 1/2); %ilnear selection %2^(-iRepeat_clu+1)
%                     vlRedo_spk = ismember(res.viClu, find(vlRedo_clu));
%                     S_clu_A = postCluster_(cluster_spacetime_(S0, obj.hCfg, ~vlRedo_spk), obj.hCfg);
%                     S_clu_B = postCluster_(cluster_spacetime_(S0, obj.hCfg, vlRedo_spk), obj.hCfg);
%                     res.viClu(~vlRedo_spk) = S_clu_A.viClu;
%                     res.viClu(vlRedo_spk) = S_clu_B.viClu + max(S_clu_A.viClu);
% 
%                 otherwise
%                     return;
%             end
%         end
    end
end

