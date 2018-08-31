%--------------------------------------------------------------------------
function spikeWaveforms1 = mn2tn_wav_spk2_(mnWav1, spikeSites, spikeTimes, P)
    nSpks = numel(spikeSites);
    nSites = numel(P.chanMap);
    spkLim_wav = P.spkLim;
    nSites_spk = (P.maxSite * 2) + 1;
    spikeWaveforms1 = zeros(diff(spkLim_wav) + 1, nSites_spk, nSpks, 'like', mnWav1);
    for iSite = 1:nSites
        viiSpk11 = find(spikeSites == iSite);
        if isempty(viiSpk11), continue; end
        spikeTimes11 = spikeTimes(viiSpk11); %already sorted by time
        viSite11 = P.miSites(:,iSite);
        spikeWaveforms1(:,:,viiSpk11) = permute(gather_(mr2tr3_(mnWav1, spkLim_wav, spikeTimes11, viSite11)), [1,3,2]); %raw
    end
end % function
