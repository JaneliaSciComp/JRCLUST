%--------------------------------------------------------------------------
% 10/22/17 JJJ
function [mrWav_clu1, viSite_clu1, mrWav_lo_clu1, mrWav_hi_clu1] = clu_wav_(S_clu, tnWav_, iClu, S0)
    if nargin<4, S0 = get(0, 'UserData'); end
    fUseCenterSpk = 0; % set to zero to use all spikes
    nSamples_max = 1000;

    fDrift_merge = get_set_(S0.P, 'fDrift_merge', 0);
    [mrWav_clu1, viSite_clu1, mrWav_lo_clu1, mrWav_hi_clu1] = deal([]);
    iSite_clu1 = S_clu.viSite_clu(iClu);
    viSite_clu1 = S0.P.miSites(:,iSite_clu1);
    viSpk_clu1 = S_clu.cviSpk_clu{iClu}; %
    viSite_spk1 = S0.viSite_spk(viSpk_clu1);
    vlCentered_spk1 = iSite_clu1 == viSite_spk1;
    if fUseCenterSpk
        viSpk_clu1 = viSpk_clu1(vlCentered_spk1);
        viSite_spk1 = viSite_spk1(vlCentered_spk1);
    end
    if isempty(viSpk_clu1), return; end
    if ~fDrift_merge
        viSpk_clu2 = spk_select_mid_(viSpk_clu1, S0.viTime_spk, S0.P);
        mrWav_clu1 = mean(single(tnWav_(:,:,viSpk_clu2)), 3);
        mrWav_clu1 = meanSubt_(mrWav_clu1); %122717 JJJ
        return;
    end

    vrPosY_spk1 = S0.mrPos_spk(viSpk_clu1,2); %position based quantile
    vrYLim = quantile(vrPosY_spk1, [0,1,2,3]/3);
    [viSpk_clu_, viSite_clu_] = spk_select_pos_(viSpk_clu1, vrPosY_spk1, vrYLim(2:3), nSamples_max, viSite_spk1);
    mrWav_clu1 = nanmean_int16_(tnWav_(:,:,viSpk_clu_), 3, fUseCenterSpk, iSite_clu1, viSite_clu_, S0.P); % * S0.P.uV_per_bit;

    if nargout > 2
        [viSpk_clu_, viSite_clu_] = spk_select_pos_(viSpk_clu1, vrPosY_spk1, vrYLim(1:2), nSamples_max, viSite_spk1);
        mrWav_lo_clu1 = nanmean_int16_(tnWav_(:,:,viSpk_clu_), 3, fUseCenterSpk, iSite_clu1, viSite_clu_, S0.P);

        [viSpk_clu_, viSite_clu_] = spk_select_pos_(viSpk_clu1, vrPosY_spk1, vrYLim(3:4), nSamples_max, viSite_spk1);
        mrWav_hi_clu1 = nanmean_int16_(tnWav_(:,:,viSpk_clu_), 3, fUseCenterSpk, iSite_clu1, viSite_clu_, S0.P);
    end
end %func
