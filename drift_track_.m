%--------------------------------------------------------------------------
function [vrTime_drift, vrPos_drift] = drift_track_(S_clu, vrPosY_spk, P)
    viTime_spk = double(get0_('viTime_spk'));
    vrAmp_spk = get0_('vrAmp_spk');
    nStep = round(get_set_(P, 'tbin_drift', 2) * P.sRateHz);
    nBins = ceil(max(viTime_spk) / nStep);
    vrTime_drift = ((1:nBins) * nStep - round(nStep/2)) / P.sRateHz;
    viiTime_spk_drift = ceil(viTime_spk / nStep); % replaces discritize function
    viClu_use = find(S_clu.vrSnr_clu(:) > quantile(S_clu.vrSnr_clu, .5)); % & S_clu.vnSpk_clu(:) > quantile(S_clu.vnSpk_clu, .5));
    nClu_use = numel(viClu_use);

    % spike based drift estimation
    % pos_lim = quantile(vrPosY_spk, [.25, .75]);
    % viSpk_use = find(vrAmp_spk < median(vrAmp_spk) & vrPosY_spk(:) > pos_lim(1) & vrPosY_spk(:) < pos_lim(2));
    % vrPos_drift = accumarray(viiTime_spk_drift(viSpk_use), double(vrPosY_spk(viSpk_use)),[nBins,1],@mean,nan); %viiTime_clu1, vrPosY_spk(viSpk_clu1);

    % cluster based drift estimation
    mrPos_clu_drift = nan(nBins, nClu_use, 'single');
    for iClu1=1:nClu_use
        iClu = viClu_use(iClu1);
        viSpk_clu1 = S_clu.cviSpk_clu{iClu};
        viiTime_clu1 = viiTime_spk_drift(viSpk_clu1);
        vrPos_clu1 = double(vrPosY_spk(viSpk_clu1));
        mrPos_clu_drift(:,iClu1) = accumarray(viiTime_clu1,vrPos_clu1,[nBins,1],@median,nan); %viiTime_clu1, vrPosY_spk(viSpk_clu1);
    end
    % mrPos_clu_drift = medfilt1(mrPos_clu_drift, 8);
    vrActivity_clu = mean(~isnan(mrPos_clu_drift));
    vrPos_drift = nanmean(mrPos_clu_drift(:,vrActivity_clu>quantile(vrActivity_clu, .9)),2);

    figure; plot(vrTime_drift, vrPos_drift)
end %func
