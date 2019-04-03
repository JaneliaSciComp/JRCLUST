function res = computeRho(dRes, res, hCfg)
    %COMPUTERHO Compute rho values for spike features
    res.rhoCutSite = zeros(hCfg.nSites, 1, 'single');
    res.rhoCutGlobal = [];

    % compute rho cutoff globally
    if hCfg.useGlobalDistCut
        res.rhoCutGlobal = estRhoCutGlobal(dRes, hCfg);
    end

    hCfg.updateLog('computeRho', 'Computing rho', 1, 0);

    % create CUDA kernel
    chunkSize = 16;
    nC_max = 45;
    if hCfg.useGPU
        ptxFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_rho.ptx');
        cuFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_rho.cu');
        rhoCK = parallel.gpu.CUDAKernel(ptxFile, cuFile);
        rhoCK.ThreadBlockSize = [hCfg.nThreadsGPU, 1];
        rhoCK.SharedMemorySize = 4 * chunkSize * (2 + nC_max + 2*hCfg.nThreadsGPU);
    else
        rhoCK = [] ;
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
        elseif isfield(dRes, 'spikesBySite3')
            spikeData.spikes3 = dRes.spikesBySite3(dRes.spikesBySite3==iSite);
        end

        [siteFeatures, ~, n1, n2, spikeOrder] = jrclust.features.getSiteFeatures(dRes.spikeFeatures, iSite, spikeData, hCfg);
        if hCfg.useGPU
            rhoCK.GridSize = [ceil(n1/chunkSize^2), chunkSize]; %MaxGridSize: [2.1475e+09 65535 65535]
        end

        if isempty(siteFeatures)
            continue;
        end

        siteFeatures = jrclust.utils.tryGpuArray(siteFeatures, hCfg.useGPU);
        spikeOrder = jrclust.utils.tryGpuArray(spikeOrder, hCfg.useGPU);

        if isempty(res.rhoCutGlobal) % estimate rhoCut in CPU
            siteCut = estRhoCutSite(siteFeatures, spikeOrder, n1, n2, hCfg);
        else
            siteCut = res.rhoCutGlobal.^2;
        end


        siteRho = computeRhoSite(siteFeatures, spikeOrder, n1, n2, siteCut, rhoCK, hCfg);

        res.spikeRho(spikeData.spikes1) = jrclust.utils.tryGather(siteRho);
        res.rhoCutSite(iSite) = jrclust.utils.tryGather(siteCut);
        clear siteFeatures spikeOrder siteRho;
        hCfg.updateLog('rhoSite', sprintf('Site %d: rho cutoff, %0.2f; average rho: %0.5f (%d spikes)', iSite, siteCut, mean(res.spikeRho), n1), 0, 0);
    end

    hCfg.updateLog('computeRho', 'Finished computing rho', 0, 1);
end

%% LOCAL FUNCTIONS
function rho = computeRhoSite(siteFeatures, spikeOrder, n1, n2, rhoCut, rhoCK, hCfg)
    %COMPUTERHOSITE Compute site-wise rho for spike features
    [nC, n12] = size(siteFeatures); % nc is constant with the loop
    dn_max = int32(round((n1 + n2) / hCfg.nClusterIntervals));

    rhoCut = single(rhoCut);
    if hCfg.useGPU
        try
            rho = zeros(1, n1, 'single', 'gpuArray');
            consts = int32([n1, n12, nC, dn_max, hCfg.getOr('fDc_spk', 0)]);
            rho = feval(rhoCK, rho, siteFeatures, spikeOrder, consts, rhoCut);

            return;
        catch ME
            warning('CUDA kernel failed: %s', ME.message);
        end
    end

    % retry in CPU
    [siteFeatures, spikeOrder] = jrclust.utils.tryGather(siteFeatures, spikeOrder);
    spikeOrderN1 = spikeOrder(1:n1)';

    % we can quickly run out of space here, ensure this doesn't happen
    availMem = 2^33; % 8 GiB
    if jrclust.utils.typeBytes(class(spikeOrder))*n1*n12 > availMem
        stepSize = floor(availMem/n12/jrclust.utils.typeBytes(class(spikeOrder)));
        rho = zeros(1, n1, 'single');

        for iChunk = 1:stepSize:n1
            iRange = iChunk:min(iChunk+stepSize-1, n1);
            nearbyTimeChunk = abs(bsxfun(@minus, spikeOrder, spikeOrderN1(iRange))) <= dn_max;
            nearbySpaceChunk = pdist2(siteFeatures', siteFeatures(:, iRange)', 'squaredeuclidean') <= rhoCut;
            rho(iRange) = sum(nearbyTimeChunk & nearbySpaceChunk) ./ sum(nearbyTimeChunk);
        end
    else
        nearbyTime = abs(bsxfun(@minus, spikeOrder, spikeOrderN1)) <= dn_max;
        nearbySpace = (pdist2(siteFeatures', siteFeatures(:, 1:n1)', 'squaredeuclidean') <= rhoCut); % include self
        rho = sum(nearbyTime & nearbySpace);
        rho = single(rho ./ sum(nearbyTime)); % normalize
    end
end

function rhoCut = estRhoCutGlobal(dRes, hCfg, vlRedo_spk)
    %ESTRHOCUTGLOBAL Estimate rho cutoff distance over all sites
    if nargin < 3
        vlRedo_spk = [];
    end

    spikeData = struct('spikeTimes', dRes.spikeTimes, ...
                       'vlRedo_spk', vlRedo_spk);

    hCfg.updateLog('estimateCutDist', 'Estimating cutoff distance', 1, 0);

    % estimate rho cutoff for each site
    siteCuts = nan(1, hCfg.nSites);
    for iSite = 1:hCfg.nSites
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

        [siteFeatures, ~, n1, n2, spikeOrder] = jrclust.features.getSiteFeatures(dRes.spikeFeatures, iSite, spikeData, hCfg);

        if isempty(siteFeatures)
            continue;
        end

        siteCuts(iSite) = estRhoCutSite(siteFeatures, spikeOrder, n1, n2, hCfg);
        hCfg.updateLog('rhoCutSite', sprintf('Estimated rho cutoff for site %d: %0.2f', iSite, res.rhoCutSite(iSite)), 0, 0);
    end

    rhoCut = sqrt(abs(quantile(siteCuts, .5)));

    hCfg.updateLog('estimateCutDist', sprintf('Finished estimating cutoff distance: %0.2f', rhoCut), 0, 1);
end

function rhoCut = estRhoCutSite(siteFeatures, spikeOrder, n1, n2, hCfg)
    %ESTRHOCUTSITE Estimate site-wise rho cutoff
    if hCfg.getOr('fDc_spk', 0)
        rhoCut = (hCfg.distCut/100).^2;
        return;
    end

    nSubsample = 1000;
    [nPrimary, nSecondary] = deal(nSubsample, 4*nSubsample);

    ss = jrclust.utils.subsample(1:n1, nPrimary);
    spikeOrderPrimary = spikeOrder(ss);
    featuresPrimary = siteFeatures(:, ss);

    ss = jrclust.utils.subsample(1:n1, nSecondary);
    spikeOrder = spikeOrder(ss);
    siteFeatures = siteFeatures(:, ss);

    for iRetry = 1:2
        try
            featureDists = pdist2(siteFeatures', featuresPrimary', 'squaredeuclidean');

            fSubset = abs(bsxfun(@minus, spikeOrder, spikeOrderPrimary')) < (n1 + n2) / hCfg.nClusterIntervals;
            featureDists(~fSubset) = nan;
            break;
        catch ME
            siteFeatures = jrclust.utils.tryGather(siteFeatures);
        end
    end

    if hCfg.getOr('fDc_subsample_mode', 0)
        featureDists(featureDists <= 0) = nan;
        rhoCut = quantile(featureDists(~isnan(featureDists)), hCfg.distCut/100);
    else
        featureDists = jrclust.utils.tryGather(featureDists);
        featureDists(featureDists <=0) = nan;
        rhoCut = nanmedian(quantile(featureDists, hCfg.distCut/100));

        if isnan(rhoCut) % featureDists was completely nan
            rhoCut = quantile(featureDists(:), hCfg.distCut/100);
        end
    end
end
