%--------------------------------------------------------------------------
function [tnWav_raw, tnWav_spk, spikeTimes] = mn2tn_wav_(mnWav_raw, mnWav_spk, spikeSites, spikeTimes, P)
    nSpks = numel(spikeTimes);
    nSites = numel(P.chanMap);
    spkLim_wav = P.spkLim;
    spkLim_raw = P.spkLim_raw;
    nSites_spk = (P.maxSite * 2) + 1;
    tnWav_raw = zeros(diff(spkLim_raw) + 1, nSites_spk, nSpks, 'like', mnWav_raw);
    tnWav_spk = zeros(diff(spkLim_wav) + 1, nSites_spk, nSpks, 'like', mnWav_spk);

    % Realignment parameters
    fRealign_spk = get_set_(P, 'fRealign_spk', 0); %0,1,2
    spikeTimes = gpuArray_(spikeTimes, isGpu_(mnWav_raw));
    spikeSites = gpuArray_(spikeSites, isGpu_(mnWav_raw));
    if isempty(spikeSites)
        tnWav_raw = permute(mr2tr3_(mnWav_raw, spkLim_raw, spikeTimes), [1,3,2]);
        tnWav_spk = permute(mr2tr3_(mnWav_spk, spkLim_wav, spikeTimes), [1,3,2]);
    else
        for iSite = 1:nSites
            viiSpk11 = find(spikeSites == iSite);
            if isempty(viiSpk11), continue; end
            spikeTimes11 = spikeTimes(viiSpk11); %already sorted by time
            viSite11 = P.miSites(:,iSite);
            try
                tnWav_spk1 = mr2tr3_(mnWav_spk, spkLim_wav, spikeTimes11, viSite11);
                if fRealign_spk==1
                    [tnWav_spk1, spikeTimes11] = spkwav_realign_(tnWav_spk1, mnWav_spk, spkLim_wav, spikeTimes11, viSite11, P);
                    spikeTimes(viiSpk11) = spikeTimes11;
                elseif fRealign_spk==2
                    tnWav_spk1 = spkwav_align_(tnWav_spk1, P);
                end
                tnWav_spk(:,:,viiSpk11) = permute(tnWav_spk1, [1,3,2]);
                tnWav_raw(:,:,viiSpk11) = permute(mr2tr3_(mnWav_raw, spkLim_raw, spikeTimes11, viSite11), [1,3,2]); %raw
            catch % GPU failure
                disperr_('mn2tn_wav_: GPU failed');
            end
        end
    end
end %func
