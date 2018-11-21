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

            nSites = numel(obj.hCfg.siteMap);
            nSpikes = numel(dRes.spikeTimes);

            res.spikeRho = zeros(nSpikes, 1, 'single');
            res.spikeDelta = zeros(nSpikes, 1, 'single');
            res.spikeNeigh = zeros(nSpikes, 1, 'uint32');
            res.rhoCutSite = zeros(nSites, 1, 'single');
            res.rhoCutGlobal = [];

            % clear memory
            % cuda_rho_(); % set nC_ to 0 and return
            % cuda_delta_(); % set nC_ to 0 and return
            % if get_set_(obj.hCfg, 'fDenoise_fet', 0) % not in default prm; deprecated?
            %     trFet_spk = denoise_fet_(trFet_spk, vlRedo_spk);
            % end

            % compute rho cutoff globally
            if get_set_(obj.hCfg, 'fDc_global', true)
                res.rhoCutGlobal = obj.estRhoCutGlobal(dRes);
            end

            % compute rho
            if obj.hCfg.verbose
                fprintf('Calculating Rho\n\t');
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

                siteRho = obj.cuda_rho_(siteFeatures, spikeOrder, n1, n2, siteCut);

%                 if ~isempty(vlRedo_spk), viSpk_site_ = viSpk_site_(vlRedo_spk(viSpk_site_)); end
                res.spikeRho(spikeData.spikes1) = jrclust.utils.tryGather(siteRho);
                res.rhoCutSite(iSite) = jrclust.utils.tryGather(siteCut);
                [siteFeatures, spikeOrder, siteRho] = deal([]);
                fprintf('.');
            end

            if obj.hCfg.verbose
                fprintf('\n\ttook %0.1fs\n', toc(t1));
            end

            res.runtime = toc(t0);
        end
    end

    % UTILITY METHODS
    methods (Access=protected, Hidden)
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
            if get_set_(obj.hCfg, 'fDc_spk', 0)
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
                catch
                    siteFeatures = jrclust.utils.tryGather(siteFeatures);
                end
            end

            if get_set_(obj.hCfg, 'fDc_subsample_mode', 0)
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

        function rho = cuda_rho_(obj, siteFeatures, spikeOrder, n1, n2, rhoCut)
%             persistent nC_
            if nargin == 0
                obj.rhoCK = [];
            end

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
                    consts = int32([n1, n12, nC, dn_max, get_set_(obj.hCfg, 'fDc_spk', 0)]);
                    rho = feval(obj.rhoCK, rho, siteFeatures, spikeOrder, consts, rhoCut);

                    return;
                catch ME
                    disperr_('CUDA kernel failed. Re-trying in CPU.');
                    obj.rhoCK = [];
                end
            end

            % retry in CPU
            viiSpk1_ord = spikeOrder(1:n1)';
            mlKeep12_ = abs(bsxfun(@minus, spikeOrder, viiSpk1_ord)) <= dn_max;
            if FLAG_FIXN == 0
                rho = sum(pdist2(siteFeatures', siteFeatures(:,1:n1')) < rhoCut & mlKeep12_) - 1; %do not include self
                rho = single(rho ./ sum(mlKeep12_));
            else
                mlKeep12_(1:min(2*dn_max+1,n12), viiSpk1_ord <= min(dn_max,n1)) = 1;
                mlKeep12_(max(n12-2*dn_max,1):n12,  viiSpk1_ord >= max(n1-dn_max,1)) = 1;
                rho = sum(pdist2(siteFeatures', siteFeatures(:,1:n1')) < rhoCut & mlKeep12_) - 1; %do not include self
                rho = single(rho ./ (2*dn_max));
            end
        end

        function [vrDelta1, viNneigh1] = cuda_delta_(obj, mrFet12, viiSpk12_ord, viiRho12_ord, n1, n2, dc2)
            % Ultimately use CUDA to do this distance computation
            if nargin==0, nC_ = 0; return; end
            if isempty(nC_), nC_ = 0; end
            [nC, n12] = size(mrFet12); %nc is constant with the loop
            dn_max = int32(round((n1+n2) / obj.hCfg.nTime_clu));
            if obj.hCfg.fGpu
                try
                    if (nC_ ~= nC) % create cuda kernel
                        nC_ = nC;
                        obj.deltaCK = parallel.gpu.CUDAKernel('jrc_cuda_delta.ptx','jrc_cuda_delta.cu');
                        obj.deltaCK.ThreadBlockSize = [obj.hCfg.nThreads, 1];
                        obj.deltaCK.SharedMemorySize = 4 * obj.chunkSize * (3 + obj.nC_max + 2*obj.hCfg.nThreads); % @TODO: update the size
                    end
                    obj.deltaCK.GridSize = [ceil(n1 / obj.chunkSize / obj.chunkSize), obj.chunkSize]; %MaxGridSize: [2.1475e+09 65535 65535]
                    vrDelta1 = zeros([1, n1], 'single', 'gpuArray');
                    viNneigh1 = zeros([1, n1], 'uint32', 'gpuArray');
                    vnConst = int32([n1, n12, nC, dn_max, get_set_(obj.hCfg, 'fDc_spk', 0)]);
                    [vrDelta1, viNneigh1] = feval(obj.deltaCK, vrDelta1, viNneigh1, mrFet12, viiSpk12_ord, viiRho12_ord, vnConst, dc2);
                    % [vrDelta1_, viNneigh1_] = deal(vrDelta1, viNneigh1);
                    return;
                catch
                    disperr_('CUDA kernel failed. Re-trying in CPU.');
                    nC_ = 0;
                end
            end

            mrDist12_ = pdist2(mrFet12', mrFet12(:,1:n1')); %not sqrt
            mlRemove12_ = bsxfun(@ge, viiRho12_ord, viiRho12_ord(1:n1)') ...
            | abs(bsxfun(@minus, viiSpk12_ord, viiSpk12_ord(1:n1)')) > dn_max;
            mrDist12_(mlRemove12_) = nan;
            [vrDelta1, viNneigh1] = min(mrDist12_);
        end
    end
end

