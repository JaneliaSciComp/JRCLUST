%--------------------------------------------------------------------------
function [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, sitesOfInterest)
    % show spikes excluding the clusters excluding clu1 and 2
    P = S0.P;
    S_clu = S0.S_clu;
    primaryCluster = S0.primarySelectedCluster;
    secondaryCluster = S0.iCluPaste;

    % select (indices of) subset of spikes on our sites of interest
    siteSpikes = find(ismember(S0.spikeSites, sitesOfInterest));
    siteSpikeTimes = S0.spikeTimes(siteSpikes);

    %time filter
    if ~isfield(P, 'tlim_proj'), P.tlim_proj = []; end
    if ~isempty(P.tlim_proj)
        nlim_proj = round(P.tlim_proj * P.sampleRateHz);
        windowSpikes = find(siteSpikeTimes >= nlim_proj(1) & siteSpikeTimes <= nlim_proj(end));
        siteSpikes = siteSpikes(windowSpikes);
        siteSpikeTimes = siteSpikeTimes(windowSpikes);
    end

    siteClusters = S_clu.spikeClusters(siteSpikes);
    backgroundSpikes = randomSelect_(siteSpikes, P.nShow_proj*2);
    foregroundSpikes = randomSelect_(siteSpikes(siteClusters == primaryCluster), P.nShow_proj);

    if ~isempty(secondaryCluster)
        secondaryForegroundSpikes = randomSelect_(siteSpikes(siteClusters == secondaryCluster), P.nShow_proj);
    else
        [mrMin2, mrMax2] = deal([]);
    end

    switch lower(P.displayFeature)
        case {'pca'} % pca_pv_spk_ (not pca_pv_clu_)
            [mrPv1, mrPv2] = pca_pv_spk_(S_clu.spikesByCluster{primaryCluster}, sitesOfInterest);
            [mrMin0, mrMax0] = pca_pc_spk_(backgroundSpikes, sitesOfInterest, mrPv1, mrPv2); % get all spikes whose center lies in certain range
            [mrMin1, mrMax1] = pca_pc_spk_(foregroundSpikes, sitesOfInterest, mrPv1, mrPv2); % get all spikes whose center lies in certain range

            if ~isempty(secondaryCluster)
                [mrMin2, mrMax2] = pca_pc_spk_(secondaryForegroundSpikes, sitesOfInterest, mrPv1, mrPv2);
            end

        case {'ppca', 'private pca'} % pca_pv_clu_ (not pca_pv_spk_)
            [mrPv1, mrPv2] = pca_pv_clu_(sitesOfInterest, primaryCluster, secondaryCluster);
            [mrMin0, mrMax0] = pca_pc_spk_(backgroundSpikes, sitesOfInterest, mrPv1, mrPv2); % get all spikes whose center lies in certain range
            [mrMin1, mrMax1] = pca_pc_spk_(foregroundSpikes, sitesOfInterest, mrPv1, mrPv2); % get all spikes whose center lies in certain range

            if ~isempty(secondaryCluster)
                [mrMin2, mrMax2] = pca_pc_spk_(secondaryForegroundSpikes, sitesOfInterest, mrPv1, mrPv2);
            end

        otherwise % generic
            [mrMin0, mrMax0] = getFet_spk_(backgroundSpikes, sitesOfInterest, S0); % get all spikes whose center lies in certain range
            [mrMin1, mrMax1] = getFet_spk_(foregroundSpikes, sitesOfInterest, S0); % get all spikes whose center lies in certain range

            if ~isempty(secondaryCluster)
                [mrMin2, mrMax2] = getFet_spk_(secondaryForegroundSpikes, sitesOfInterest, S0);
            end
    end % switch

    [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = multifun_(@(x) abs(x), mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2);
end % func
