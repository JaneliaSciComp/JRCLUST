%--------------------------------------------------------------------------
function trWav1 = trWav_clu_(iClu1, nSpk_show)
    % Get a subset of spike waveforms of a cluster

    if nargin<2, nSpk_show=inf; end
    S0 = get(0, 'UserData');
    P = S0.P;
    tnWav_ = get_spkwav_(P);
    [viSpk_clu1, viiSpk_clu1] = S0.S_clu.getCenteredSpikes(iClu1);
    viSpk_clu1 = randomSelect_(viSpk_clu1, nSpk_show);
    if P.fWav_raw_show
        trWav1 = raw2uV_(tnWav_(:,:,viSpk_clu1), P);
    else
        trWav1 = tnWav2uV_(tnWav_(:,:,viSpk_clu1), P);
    end
end %func
