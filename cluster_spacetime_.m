%--------------------------------------------------------------------------
function S_clu = cluster_spacetime_(S0, P, vlRedo_spk)
    % this clustering is natively suited for 2D electrode arrays and drifting dataset
    % There is no x,y position in the clustering dataset
    global spikeFeatures

    if ~all(size(spikeFeatures) ~= S0.featureDims)
        spikeFeatures = getSpikeFeatures(P);
    end

    if ~isfield(P, 'CHUNK'), P.CHUNK = 16; end
    if ~isfield(P, 'fTwoStep'), P.fTwoStep = 0; end
    if ~isfield(P, 'mrSiteXY'), P.mrSiteXY = []; end
    if ~isfield(P, 'minClusterSize'), P.minClusterSize = []; end
    if nargin<3, vlRedo_spk=[]; end

    % g = gpuDevice();
    t_func = tic;
    nSites = numel(P.chanMap);
    nSpk = numel(S0.spikeTimes);
    vrRho = zeros(nSpk, 1, 'single');
    vrDelta = zeros(nSpk, 1, 'single');
    viNneigh = zeros(nSpk, 1, 'uint32');
    vrDc2_site = zeros(nSites, 1, 'single');
    nTime_clu = getOr(P, 'nTime_clu', 1);
    P.nTime_clu = nTime_clu;
    P.dc_subsample = 1000;
    P.useGPU = 1;

    % clear memory
    cuda_rho_();
    cuda_delta_();
    if getOr(P, 'fDenoise_fet', 0)
        spikeFeatures = denoise_fet_(spikeFeatures, P, vlRedo_spk);
    end

    %-----
    % Calculate dc2 (global)
    if getOr(P, 'fDc_global', 0)
        dc2 = calc_dc2_(S0, P, vlRedo_spk);
    else
        dc2 = [];
    end

    %-----
    % Calculate Rho
    fprintf('Calculating Rho\n\t'); t1=tic;
    for iSite = 1:nSites
        [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = fet12_site_(spikeFeatures, S0, P, iSite, vlRedo_spk);
        if isempty(mrFet12_), continue; end
        [mrFet12_, viiSpk12_ord_] = multifun_(@gpuArray_, mrFet12_, viiSpk12_ord_);
        if isempty(dc2)
            dc2_ = compute_dc2_(mrFet12_, viiSpk12_ord_, n1_, n2_, P); % Compute DC in CPU
        else
            dc2_ = dc2.^2;
        end
        vrRho_ = cuda_rho_(mrFet12_, viiSpk12_ord_, n1_, n2_, dc2_, P);
        viSpk_site_ = S0.cviSpk_site{iSite};
        if ~isempty(vlRedo_spk), viSpk_site_ = viSpk_site_(vlRedo_spk(viSpk_site_)); end
        vrRho(viSpk_site_) = gather_(vrRho_);
        vrDc2_site(iSite) = gather_(dc2_);
        [mrFet12_, viiSpk12_ord_, vrRho_] = deal([]);
        fprintf('.');
    end

    % [vrRho, vrDc2_site] = gather_(vrRho, vrDc2_site);
    fprintf('\n\ttook %0.1fs\n', toc(t1));

    %-----
    % Calculate Delta
    fprintf('Calculating Delta\n\t'); t2=tic;
    for iSite = 1:nSites
        [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = fet12_site_(spikeFeatures, S0, P, iSite, vlRedo_spk);
        if isempty(mrFet12_), continue; end
        viiRho12_ord_ = rankorder_(vrRho(viSpk12_), 'descend');
        [mrFet12_, viiRho12_ord_, viiSpk12_ord_] = ...
        multifun_(@gpuArray_, mrFet12_, viiRho12_ord_, viiSpk12_ord_);
        try
            [vrDelta_, viNneigh_] = cuda_delta_(mrFet12_, viiSpk12_ord_, viiRho12_ord_, n1_, n2_, vrDc2_site(iSite), P);
            [vrDelta_, viNneigh_] = gather_(vrDelta_, viNneigh_);
        catch
            disperr_(sprintf('error at site# %d', iSite));
        end
        viSpk_site_ = S0.cviSpk_site{iSite};
        if ~isempty(vlRedo_spk), viSpk_site_ = viSpk_site_(vlRedo_spk(viSpk_site_)); end
        vrDelta(viSpk_site_) = vrDelta_;
        viNneigh(viSpk_site_) = viSpk12_(viNneigh_);
        [mrFet12_, viiRho12_ord_, viiSpk12_ord_] = deal([]); %vrDelta_, viNneigh_, viSpk12_
        fprintf('.');
    end
    % Deal with nan delta
    viNan_delta = find(isnan(vrDelta));
    if ~isempty(viNan_delta)
        vrDelta(viNan_delta) = max(vrDelta);
    end
    % [vrDelta, viNneigh] = multifun_(@gather_, vrDelta, viNneigh);
    fprintf('\n\ttook %0.1fs\n', toc(t2));

    if ~isempty(vlRedo_spk)
        vrRho = vrRho(vlRedo_spk);
        vrDelta = vrDelta(vlRedo_spk);
        viNneigh = reverse_lookup_(viNneigh(vlRedo_spk), find(vlRedo_spk));
    end

    %-----
    % package
    % if P.useGPU
    %     [vrRho, vrDelta, viNneigh] = multifun_(@gather_, vrRho, vrDelta, viNneigh);
    % end
    t_runtime = toc(t_func);
    trFet_dim = size(spikeFeatures); %[1, size(mrFet1,1), size(mrFet1,2)]; %for postCluster
    [~, ordrho] = sort(vrRho, 'descend');
    S_clu = struct('rho', vrRho, 'delta', vrDelta, 'ordrho', ordrho, 'nneigh', viNneigh, ...
    'P', P, 't_runtime', t_runtime, 'halo', [], 'viiSpk', [], 'trFet_dim', trFet_dim, 'vrDc2_site', vrDc2_site);

    % figure; loglog(vrRho, vrDelta, '.');
end %func
