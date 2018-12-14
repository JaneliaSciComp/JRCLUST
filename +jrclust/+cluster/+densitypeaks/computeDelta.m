function res = computeDelta(dRes, res, hCfg)
    %COMPUTEDELTA Compute delta for spike features
    nSites = numel(hCfg.siteMap);

    if hCfg.verbose
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
        elseif isfield(dRes, 'spikeSites2')
            spikeData.spikes2 = dRes.spikeSites2(dRes.spikeSites2==iSite);
        end

        if isfield(dRes, 'spikesBySite3') && ~isempty(dRes.spikesBySite3)
            spikeData.spikes3 = dRes.spikesBySite3{iSite};
        elseif isfield(dRes, 'spikeSites3')
            spikeData.spikes3 = dRes.spikeSites3(dRes.spikeSites3==iSite);
        end

        [siteFeatures, spikes, n1, n2, spikeOrder] = jrclust.features.getSiteFeatures(dRes.spikeFeatures, iSite, spikeData, hCfg);

        if isempty(siteFeatures)
            continue;
        end

        rhoOrder = jrclust.utils.rankorder(res.spikeRho(spikes), 'descend');
        siteFeatures = jrclust.utils.tryGpuArray(siteFeatures);
        rhoOrder = jrclust.utils.tryGpuArray(rhoOrder);
        spikeOrder = jrclust.utils.tryGpuArray(spikeOrder);

        try
            [siteDelta, siteNN] = jrclust.cluster.densitypeaks.computeDeltaSite(siteFeatures, spikeOrder, rhoOrder, n1, n2, res.rhoCutSite(iSite), hCfg);
            [siteDelta, siteNN] = jrclust.utils.tryGather(siteDelta, siteNN);
        catch ME % can't continue!
            error(ME.identifier, 'Error at site %d: %s', iSite, ME.message);
        end

        % if ~isempty(vlRedo_spk), viSpk_site_ = viSpk_site_(vlRedo_spk(viSpk_site_)); end

        res.spikeDelta(spikeData.spikes1) = siteDelta;
        res.spikeNeigh(spikeData.spikes1) = spikes(siteNN);
        [~, ~, ~] = jrclust.utils.tryGather(siteFeatures, rhoOrder, spikeOrder);
        clear siteFeatures rhoOrder spikeOrder;
        if hCfg.verbose
            fprintf('.');
        end
    end

    % set delta for spikes with no nearest neighbor of higher density to
    % maximum distance, ensure they don't get grouped into another cluster
    nanDelta = find(isnan(res.spikeDelta));
    if ~isempty(nanDelta)
        res.spikeDelta(nanDelta) = max(res.spikeDelta);
    end

    if hCfg.verbose
        fprintf('\n\ttook %0.2fs\n', toc(t2));
    end
end
