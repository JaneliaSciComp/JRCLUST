function res = computeRho(dRes, res, hCfg)
    %COMPUTERHO Compute rho values for spike features
    res.rhoCutSite = zeros(hCfg.nSites, 1, 'single');
    res.rhoCutGlobal = [];

    % compute rho cutoff globally
    if hCfg.useGlobalDistCut
        res.rhoCutGlobal = estRhoCutGlobal(dRes, hCfg);
    end

    if hCfg.verbose
        fprintf('Computing rho\n\t');
        t1 = tic;
    end

    % create CUDA kernel
    chunkSize = 16;
    nC_max = 45;
    ptxFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_rho.ptx');
    cuFile = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA', 'jrc_cuda_rho.cu');
    rhoCK = parallel.gpu.CUDAKernel(ptxFile, cuFile);
    rhoCK.ThreadBlockSize = [hCfg.nThreadsGPU, 1];
    rhoCK.SharedMemorySize = 4 * chunkSize * (2 + nC_max + 2*hCfg.nThreadsGPU);

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

        if isempty(siteFeatures)
            continue;
        end

        siteFeatures = jrclust.utils.tryGpuArray(siteFeatures);
        spikeOrder = jrclust.utils.tryGpuArray(spikeOrder);

        if isempty(res.rhoCutGlobal) % estimate rhoCut in CPU
            siteCut = estRhoCutSite(siteFeatures, spikeOrder, n1, n2, hCfg);
        else
            siteCut = res.rhoCutGlobal.^2;
        end

        rhoCK.GridSize = [ceil(n1/chunkSize^2), chunkSize]; %MaxGridSize: [2.1475e+09 65535 65535]
        siteRho = computeRhoSite(siteFeatures, spikeOrder, n1, n2, siteCut, rhoCK, hCfg);

        res.spikeRho(spikeData.spikes1) = jrclust.utils.tryGather(siteRho);
        res.rhoCutSite(iSite) = jrclust.utils.tryGather(siteCut);
        clear siteFeatures spikeOrder siteRho;
        if hCfg.verbose
            fprintf('.');
        end
    end

    if hCfg.verbose
        fprintf('\n\ttook %0.1fs\n', toc(t1));
    end
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
    nearby = abs(bsxfun(@minus, spikeOrder, spikeOrderN1)) <= dn_max;
    rho = sum(nearby & pdist2(siteFeatures', siteFeatures(:,1:n1)').^2 < rhoCut) - 1; %do not include self
    rho = single(rho ./ sum(nearby)); % normalize
end

function rhoCut = estRhoCutGlobal(dRes, hCfg, vlRedo_spk)
    %ESTRHOCUTGLOBAL Estimate rho cutoff distance over all sites
    if nargin < 3
        vlRedo_spk = [];
    end

    spikeData = struct('spikeTimes', dRes.spikeTimes, ...
                       'vlRedo_spk', vlRedo_spk);

    if hCfg.verbose
        fprintf('Estimating cutoff distance\n\t');
        t1 = tic;
    end

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
    end

    rhoCut = sqrt(abs(quantile(siteCuts, .5)));

    if hCfg.verbose
        fprintf('\n\ttook %0.1fs\n', toc(t1));
    end
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
            featureDists = pdist2(siteFeatures', featuresPrimary').^2;

            fSubset = abs(bsxfun(@minus, spikeOrder, spikeOrderPrimary')) < (n1 + n2) / hCfg.nClusterIntervals;
            featureDists(~fSubset) = nan;
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

    if hCfg.verbose
        fprintf('.');
    end
end