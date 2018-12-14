function [delta, nNeigh] = computeDeltaSite(siteFeatures, spikeOrder, rhoOrder, n1, n2, rhoCut, hCfg)
    %COMPUTEDELTASITE Compute site-wise delta for spike features
    [nC, n12] = size(siteFeatures); % nc is constant with the loop
    dn_max = int32(round((n1 + n2) / hCfg.nTime_clu));

    persistent deltaCK;
    chunkSize = 16;
    nC_max = 45;

    if hCfg.fGpu
        try
            if isempty(deltaCK) % create cuda kernel
                ptxFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_delta.ptx');
                cuFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_delta.cu');
                deltaCK = parallel.gpu.CUDAKernel(ptxFile, cuFile);
                deltaCK.ThreadBlockSize = [hCfg.nThreads, 1];
                deltaCK.SharedMemorySize = 4 * chunkSize * (3 + nC_max + 2*hCfg.nThreads); % @TODO: update the size
            end

            % set every time function is called
            deltaCK.GridSize = [ceil(n1 / chunkSize / chunkSize), chunkSize]; %MaxGridSize: [2.1475e+09 65535 65535]
            delta = zeros([1, n1], 'single', 'gpuArray');
            nNeigh = zeros([1, n1], 'uint32', 'gpuArray');
            consts = int32([n1, n12, nC, dn_max, hCfg.getOr('fDc_spk', 0)]);
            [delta, nNeigh] = feval(deltaCK, delta, nNeigh, siteFeatures, spikeOrder, rhoOrder, consts, rhoCut);

            return;
        catch ME
            warning(ME.identifier, 'CUDA kernel failed: %s', ME.message);
            deltaCK = [];
        end
    end

    dists = pdist2(siteFeatures', siteFeatures(:,1:n1)').^2;
    nearby = bsxfun(@lt, rhoOrder, rhoOrder(1:n1)') & abs(bsxfun(@minus, spikeOrder, spikeOrder(1:n1)')) <= dn_max;
    dists(~nearby) = nan;
    [delta, nNeigh] = min(dists);
end