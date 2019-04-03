function res = computeDelta(dRes, res, hCfg)
    %COMPUTEDELTA Compute delta for spike features
    hCfg.updateLog('computeDelta', 'Computing delta', 1, 0);

    % create CUDA kernel
    chunkSize = 16;
    nC_max = 45;
    if hCfg.useGPU
        ptxFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_delta.ptx');
        cuFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_delta.cu');
        deltaCK = parallel.gpu.CUDAKernel(ptxFile, cuFile);
        deltaCK.ThreadBlockSize = [hCfg.nThreadsGPU, 1];
        deltaCK.SharedMemorySize = 4 * chunkSize * (3 + nC_max + 2*hCfg.nThreadsGPU);
    else
        deltaCK = [] ;
    end

    spikeData = struct('spikeTimes', dRes.spikeTimes);
    for iSite = 1:hCfg.nSites
        if isfield(dRes, 'spikesBySite')
            spikeData.spikes1 = dRes.spikesBySite{iSite};
        else
            spikeData.spikes1 = dRes.spikeSites(dRes.spikeSites==iSite);
        end

        if isfield(dRes, 'spikesBySite2') && ~isempty(dRes.spikesBySite2)
            spikeData.spikes2 = dRes.spikesBySite2{iSite};
        elseif isfield(dRes, 'spikeSites2')
            spikeData.spikes2 = dRes.spikeSites2(dRes.spikeSites2==iSite);
        end

        if isfield(dRes, 'spikesBySite3') && ~isempty(dRes.spikesBySite3)
            spikeData.spikes3 = dRes.spikesBySite3{iSite};
        elseif isfield(dRes, 'spikeSites3')
            spikeData.spikes3 = dRes.spikeSites3(dRes.spikeSites3==iSite);
        end

        [siteFeatures, spikes, n1, n2, spikeOrder] = jrclust.features.getSiteFeatures(dRes.spikeFeatures, iSite, spikeData, hCfg);
        if hCfg.useGPU
            deltaCK.GridSize = [ceil(n1/chunkSize^2), chunkSize]; % MaxGridSize: [2.1475e+09 65535 65535]
        end

        if isempty(siteFeatures)
            continue;
        end

        rhoOrder = jrclust.utils.rankorder(res.spikeRho(spikes), 'descend');
        siteFeatures = jrclust.utils.tryGpuArray(siteFeatures, hCfg.useGPU);
        rhoOrder = jrclust.utils.tryGpuArray(rhoOrder, hCfg.useGPU);
        spikeOrder = jrclust.utils.tryGpuArray(spikeOrder, hCfg.useGPU);

        try
            [siteDelta, siteNN] = computeDeltaSite(siteFeatures, spikeOrder, rhoOrder, n1, n2, res.rhoCutSite(iSite), deltaCK, hCfg);
            [siteDelta, siteNN] = jrclust.utils.tryGather(siteDelta, siteNN);
        catch ME % can't continue!
            error('Error at site %d: %s', iSite, ME.message);
        end

        res.spikeDelta(spikeData.spikes1) = siteDelta;
        res.spikeNeigh(spikeData.spikes1) = spikes(siteNN);
        [siteFeatures, rhoOrder, spikeOrder] = jrclust.utils.tryGather(siteFeatures, rhoOrder, spikeOrder); %#ok<ASGLU>
        hCfg.updateLog('deltaSite', sprintf('Site %d: median delta: %0.5f (%d spikes)', iSite, median(siteDelta), n1), 0, 0);
    end

    % set delta for spikes with no nearest neighbor of higher density to
    % maximum distance, ensure they don't get grouped into another cluster
    nanDelta = find(isnan(res.spikeDelta));
    if ~isempty(nanDelta)
        res.spikeDelta(nanDelta) = max(res.spikeDelta);
    end

    hCfg.updateLog('computeDelta', 'Finished computing delta', 0, 1);
end

%% LOCALFUNCTIONS
function [delta, nNeigh] = computeDeltaSite(siteFeatures, spikeOrder, rhoOrder, n1, n2, distCut2, deltaCK, hCfg)
    %COMPUTEDELTASITE Compute site-wise delta for spike features
    [nC, n12] = size(siteFeatures); % nc is constant with the loop
    dn_max = int32(round((n1 + n2) / hCfg.nClusterIntervals));

    if hCfg.useGPU
        try
            % set every time function is called
            delta = zeros([1, n1], 'single', 'gpuArray');
            nNeigh = zeros([1, n1], 'uint32', 'gpuArray');
            consts = int32([n1, n12, nC, dn_max, hCfg.getOr('fDc_spk', 0)]);
            [delta, nNeigh] = feval(deltaCK, delta, nNeigh, siteFeatures, spikeOrder, rhoOrder, consts, distCut2);

            return;
        catch ME
            warning('CUDA kernel failed: %s', ME.message);
        end
    end

    rhoOrderN1 = rhoOrder(1:n1)';
    spikeOrderN1 = spikeOrder(1:n1)';

    % we can quickly run out of space here, ensure this doesn't happen
    availMem = 2^33; % 8 GiB
    if jrclust.utils.typeBytes(class(spikeOrder))*n1*n12 > availMem
        stepSize = floor(availMem/n12/jrclust.utils.typeBytes(class(spikeOrder)));
        delta = zeros(1, n1, 'single');
        nNeigh = zeros(1, n1, 'uint32');

        for iChunk = 1:stepSize:n1
            iRange = iChunk:min(iChunk+stepSize-1, n1);

            % find spikes with a larger rho and are nearby enough in time
            isDenserChunk = bsxfun(@lt, rhoOrder, rhoOrderN1(iRange));
            nearbyTimeChunk = abs(bsxfun(@minus, spikeOrder, spikeOrderN1(iRange))) <= dn_max;

            distsChunk = pdist2(siteFeatures', siteFeatures(:, iRange)', 'squaredeuclidean');
            distsChunk(~(isDenserChunk & nearbyTimeChunk)) = nan;
            distsChunk = sqrt(distsChunk/distCut2); % normalize

            [delta(iRange), nNeigh(iRange)] = min(distsChunk);
        end
    else
        dists = pdist2(siteFeatures', siteFeatures(:, 1:n1)', 'squaredeuclidean');
        isDenser = bsxfun(@lt, rhoOrder, rhoOrder(1:n1)');
        nearbyInTime = abs(bsxfun(@minus, spikeOrder, spikeOrderN1)) <= dn_max;

        dists(~(isDenser & nearbyInTime)) = nan;
        dists = sqrt(dists/distCut2); % normalize

        [delta, nNeigh] = min(dists);
    end

    % spikes which are maximally dense get SINGLE_INF
    maximallyDense = isnan(delta);
    delta(maximallyDense) = sqrt(3.402E+38/distCut2); % to (more or less) square with CUDA kernel's SINGLE_INF
    nNeigh(maximallyDense) = find(maximallyDense);
end
