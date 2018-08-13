%--------------------------------------------------------------------------
function [spikeTraces, spikeWaveforms, spikeTimes] = mn2tn_wav_(rawSamples, mnWav_spk, spikeSites, spikeTimes, P)
    nSpikes = numel(spikeTimes);
    nSites = numel(P.chanMap);
    spkLim_wav = P.spkLim;
    spkLim_raw = P.spkLim_raw;
    nSites_spk = (P.maxSite * 2) + 1;
    spikeTraces = zeros(diff(spkLim_raw) + 1, nSites_spk, nSpikes, 'like', rawSamples);
    spikeWaveforms = zeros(diff(spkLim_wav) + 1, nSites_spk, nSpikes, 'like', mnWav_spk);

    % Realignment parameters
    fRealign_spk = get_set_(P, 'fRealign_spk', 0); %0,1,2
    spikeTimes = gpuArray_(spikeTimes, isGpu_(rawSamples));
    spikeSites = gpuArray_(spikeSites, isGpu_(rawSamples));

    if isempty(spikeSites)
        spikeTraces = permute(mr2tr3_(rawSamples, spkLim_raw, spikeTimes), [1,3,2]);
        spikeWaveforms = permute(mr2tr3_(mnWav_spk, spkLim_wav, spikeTimes), [1,3,2]);
    else
        for iSite = 1:nSites
            viiSpk11 = find(spikeSites == iSite);
            if isempty(viiSpk11), continue; end
            spikeTimes11 = spikeTimes(viiSpk11); %already sorted by time
            viSite11 = P.miSites(:,iSite);
            try
                spikeWaveforms1 = mr2tr3_(mnWav_spk, spkLim_wav, spikeTimes11, viSite11);
                if fRealign_spk==1
                    [spikeWaveforms1, spikeTimes11] = spkwav_realign_(spikeWaveforms1, mnWav_spk, spkLim_wav, spikeTimes11, viSite11, P);
                    spikeTimes(viiSpk11) = spikeTimes11;
                elseif fRealign_spk==2
                    spikeWaveforms1 = spkwav_align_(spikeWaveforms1, P);
                end

                spikeWaveforms(:,:,viiSpk11) = permute(spikeWaveforms1, [1,3,2]);
                spikeTraces(:,:,viiSpk11) = permute(mr2tr3_(rawSamples, spkLim_raw, spikeTimes11, viSite11), [1,3,2]); %raw
            catch % GPU failure
                disperr_('mn2tn_wav_: GPU failed');
            end
        end
    end
end %func
