%--------------------------------------------------------------------------
function tnWav_spk1 = mn2tn_wav_spk2_(mnWav1, viSite_spk, viTime_spk, P)
    nSpks = numel(viSite_spk);
    nSites = numel(P.viSite2Chan);
    spkLim_wav = P.spkLim;
    nSites_spk = (P.maxSite * 2) + 1;
    tnWav_spk1 = zeros(diff(spkLim_wav) + 1, nSites_spk, nSpks, 'like', mnWav1);
    for iSite = 1:nSites
        viiSpk11 = find(viSite_spk == iSite);
        if isempty(viiSpk11), continue; end
        viTime_spk11 = viTime_spk(viiSpk11); %already sorted by time
        viSite11 = P.miSites(:,iSite);
        tnWav_spk1(:,:,viiSpk11) = permute(jrclust.utils.tryGather(jrclust.utils.extractWindows(mnWav1, spkLim_wav, viTime_spk11, viSite11)), [1,3,2]); %raw
    end
end %func
