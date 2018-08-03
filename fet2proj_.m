%--------------------------------------------------------------------------
function [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, viSites0)
    % show spikes excluding the clusters excluding clu1 and 2
    P = S0.P;
    S_clu = S0.S_clu;
    iClu1 = S0.iCluCopy;
    iClu2 = S0.iCluPaste;

    % select subset of spikes
    viSpk0 = find(ismember(S0.spikeSites, viSites0)); % TODO: rename viSpk0 to spikeIndices
    viTime0 = S0.spikeTimes(viSpk0);

    %time filter
    if ~isfield(P, 'tlim_proj'), P.tlim_proj = []; end
    if ~isempty(P.tlim_proj)
        nlim_proj = round(P.tlim_proj * P.sRateHz);
        viSpk01 = find(viTime0>=nlim_proj(1) & viTime0<=nlim_proj(end));
        viSpk0 = viSpk0(viSpk01);
        viTime0 = viTime0(viSpk01);
    end

    viClu0 = S_clu.spikeClusters(viSpk0);
    viSpk00 = randomSelect_(viSpk0, P.nShow_proj*2); % background spikes
    viSpk01 = randomSelect_(viSpk0(viClu0 == iClu1), P.nShow_proj); % foreground spikes

    if ~isempty(iClu2)
        viSpk02 = randomSelect_(viSpk0(viClu0 == iClu2), P.nShow_proj);
    else
        [mrMin2, mrMax2] = deal([]);
    end

    switch lower(P.displayFeature)
        case {'pca'} % channel by channel pca. do it by channel
            % determine pca vector from cluster 1
            [mrPv1, mrPv2] = pca_pv_spk_(S_clu.cviSpk_clu{iClu1}, viSites0);
            [mrMin0, mrMax0] = pca_pc_spk_(viSpk00, viSites0, mrPv1, mrPv2); %get all spikes whose center lies in certain range
            [mrMin1, mrMax1] = pca_pc_spk_(viSpk01, viSites0, mrPv1, mrPv2); %get all spikes whose center lies in certain range

            if ~isempty(iClu2)
                [mrMin2, mrMax2] = pca_pc_spk_(viSpk02, viSites0, mrPv1, mrPv2);
            end

        case {'ppca', 'private pca'} %channel by channel pca. do it by channel
            % determine pca vector from cluster 1
            [mrPv1, mrPv2] = pca_pv_clu_(viSites0, iClu1, iClu2);
            [mrMin0, mrMax0] = pca_pc_spk_(viSpk00, viSites0, mrPv1, mrPv2); %get all spikes whose center lies in certain range
            [mrMin1, mrMax1] = pca_pc_spk_(viSpk01, viSites0, mrPv1, mrPv2); %get all spikes whose center lies in certain range

            if ~isempty(iClu2)
                [mrMin2, mrMax2] = pca_pc_spk_(viSpk02, viSites0, mrPv1, mrPv2);
            end

        otherwise % generic
            [mrMin0, mrMax0] = getFet_spk_(viSpk00, viSites0, S0); %get all spikes whose center lies in certain range
            [mrMin1, mrMax1] = getFet_spk_(viSpk01, viSites0, S0); %get all spikes whose center lies in certain range

            if ~isempty(iClu2)
                [mrMin2, mrMax2] = getFet_spk_(viSpk02, viSites0, S0);
            end
    end %switch

    [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = multifun_(@(x) abs(x), mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2);
end %func
