%--------------------------------------------------------------------------
function S_clu = S_clu_update_(S_clu, viClu1, P)
    % update cluster waveform and self correlation score
    % mrWav not needed
    S0 = get(0, 'UserData');

    % find clu center
    for iClu = 1:numel(viClu1)
        iClu1 = viClu1(iClu);
        viSpk_clu1 = find(S_clu.viClu == iClu1);
        S_clu.cviSpk_clu{iClu1} = viSpk_clu1;
        S_clu.viSite_clu(iClu1) = mode(S0.viSite_spk(viSpk_clu1));
        S_clu.vnSpk_clu(iClu1) = numel(viSpk_clu1);
    end

    % update mean waveform
    S_clu = S_clu_wav_(S_clu, viClu1);
    % [~, S_clu.viSite_clu(iClu1)] = min(S_clu.tmrWav_clu(1-P.spkLim(1),:,iClu1));
    % S_clu.viSite_clu(iClu1) = mode(viSite_spk(viSpk_clu1));
    vrSelfCorr_clu = get_diag_(S_clu.mrWavCor);
    S_clu.mrWavCor = S_clu_wavcor_(S_clu, P, viClu1);
    S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, vrSelfCorr_clu);
    for iClu = 1:numel(viClu1)
        iClu1 = viClu1(iClu);
        S_clu.mrWavCor(iClu1,iClu1) = S_clu_self_corr_(S_clu, iClu1, S0);
    end
    S_clu = S_clu_position_(S_clu, viClu1);
    S_clu = S_clu_quality_(S_clu, P, viClu1);

    % [S_clu, S0] = S_clu_commit_(S_clu, 'S_clu_update_');
end %func
