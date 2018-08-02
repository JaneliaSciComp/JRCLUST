%--------------------------------------------------------------------------
function export_tnWav_spk_(h,e)
    % Export all spike waveforms from selected cluster

    S0 = get(0, 'UserData');
    P = S0.P;
    tnWav_spk = get_spkwav_(P, 0);
    tnWav_raw = get_spkwav_(P, 1);
    S_clu = S0.S_clu;
    iClu1 = S0.iCluCopy;
    viSpk1 = S_clu.cviSpk_clu{iClu1};
    nSpk1 = numel(viSpk1);
    nSites = numel(P.chanMap);
    dimm_spk1 = size(tnWav_spk);    dimm_spk1(2) = nSites;
    dimm_raw1 = size(tnWav_raw);    dimm_raw1(2) = nSites;
    tnWav_spk1 = zeros(dimm_spk1, 'like', tnWav_spk);
    tnWav_raw1 = zeros(dimm_raw1, 'like', tnWav_raw);
    miSites_spk1 = P.miSites(:, S0.spikeSites);
    for iSpk = 1:nSpk1
        tnWav_spk1(:, miSites_spk1(:,iSpk), iSpk) = tnWav_spk(:,:,iSpk);
        tnWav_raw1(:, miSites_spk1(:,iSpk), iSpk) = tnWav_raw(:,:,iSpk);
    end
    eval(sprintf('tnWav_spk_clu%d = tnWav_spk1;', iClu1));
    eval(sprintf('assignWorkspace_(tnWav_spk_clu%d);', iClu1));
    eval(sprintf('tnWav_raw_clu%d = tnWav_raw1;', iClu1));
    eval(sprintf('assignWorkspace_(tnWav_raw_clu%d);', iClu1));
end %func
