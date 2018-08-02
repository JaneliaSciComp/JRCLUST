%--------------------------------------------------------------------------
function S_cluInfo = get_cluInfo_(iClu)

    % determine cluster position
    if isempty(iClu), S_cluInfo=[]; return; end
    [S0, P, S_clu] = get0_();

    iSite1 = S_clu.clusterSites(iClu);
    viSite = P.miSites(:, iSite1);

    xyPos = [S_clu.vrPosX_clu(iClu), S_clu.vrPosY_clu(iClu)];
    vcPos = sprintf('Unit %d (x,y):(%0.1f, %0.1f)[pix]', iClu, xyPos/P.um_per_pix);
    if P.fWav_raw_show
        mrWav_clu = S_clu.tmrWav_raw_clu(:,viSite,iClu);
    else
        mrWav_clu = S_clu.tmrWav_clu(:,viSite,iClu);
    end
    trWav = trWav_clu_(iClu, P.nSpk_show * 1);
    if P.fWav_raw_show
        trWav = fft_lowpass_(trWav, get_set_(P, 'fc_spkwav_show', []), P.sRateHz);
    end
    S_cluInfo = makeStruct_(xyPos, iClu, mrWav_clu, viSite, vcPos, trWav);
    try
        S_cluInfo.l_ratio = S_clu.vrLRatio_clu(iClu);
        S_cluInfo.isi_ratio = S_clu.vrIsiRatio_clu(iClu);
        S_cluInfo.iso_dist = S_clu.vrIsoDist_clu(iClu);
        S_cluInfo.snr = S_clu.vrSnr_clu(iClu);
        S_cluInfo.uVmin = S_clu.vrVmin_uv_clu(iClu);
        S_cluInfo.uVpp = S_clu.vrVpp_uv_clu(iClu);
    catch
    end
end %func
