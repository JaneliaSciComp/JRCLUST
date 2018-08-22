%--------------------------------------------------------------------------
% 10/22/17 JJJ
function [mrWav_clu1, clusterSites, mrWav_lo_clu1, mrWav_hi_clu1] = clu_wav_(S_clu, tnWav_, iClu, S0)
    if nargin < 4
        S0 = get(0, 'UserData');
    end

    fUseCenterSpk = 0; % set to zero to use all spikes
    nSamples_max = 1000;

    fDrift_merge = getOr(S0.P, 'fDrift_merge', 0);
    [mrWav_clu1, mrWav_lo_clu1, mrWav_hi_clu1] = deal([]);
    centerSite = S_clu.clusterSites(iClu);
    clusterSites = S0.P.miSites(:, centerSite);
    clusterSpikes = S_clu.spikesByCluster{iClu}; %
    spikeSites1 = S0.spikeSites(clusterSpikes);

    if fUseCenterSpk
        vlCentered_spk1 = (centerSite == spikeSites1);
        clusterSpikes = clusterSpikes(vlCentered_spk1);
        spikeSites1 = spikeSites1(vlCentered_spk1);
    end

    if isempty(clusterSpikes)
        return;
    end

    if ~fDrift_merge
        viSpk_clu2 = spk_select_mid_(clusterSpikes, S0.spikeTimes, S0.P);
        mrWav_clu1 = mean(single(tnWav_(:,:,viSpk_clu2)), 3);
        mrWav_clu1 = meanSubt_(mrWav_clu1); %122717 JJJ

        return;
    end

    vrPosY_spk1 = S0.mrPos_spk(clusterSpikes,2); %position based quantile
    vrYLim = quantile(vrPosY_spk1, [0,1,2,3]/3);
    [viSpk_clu_, clusterSites_] = spk_select_pos_(clusterSpikes, vrPosY_spk1, vrYLim(2:3), nSamples_max, spikeSites1);
    mrWav_clu1 = nanmean_int16_(tnWav_(:,:,viSpk_clu_), 3, fUseCenterSpk, centerSite, clusterSites_, S0.P); % * S0.P.uV_per_bit;

    if nargout > 2
        [viSpk_clu_, clusterSites_] = spk_select_pos_(clusterSpikes, vrPosY_spk1, vrYLim(1:2), nSamples_max, spikeSites1);
        mrWav_lo_clu1 = nanmean_int16_(tnWav_(:,:,viSpk_clu_), 3, fUseCenterSpk, centerSite, clusterSites_, S0.P);

        [viSpk_clu_, clusterSites_] = spk_select_pos_(clusterSpikes, vrPosY_spk1, vrYLim(3:4), nSamples_max, spikeSites1);
        mrWav_hi_clu1 = nanmean_int16_(tnWav_(:,:,viSpk_clu_), 3, fUseCenterSpk, centerSite, clusterSites_, S0.P);
    end
end %func
