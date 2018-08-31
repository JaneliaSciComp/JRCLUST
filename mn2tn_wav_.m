%--------------------------------------------------------------------------
function [spikeTraces, spikeWaveforms, spikeTimes] = mn2tn_wav_(rawTraces, filteredTraces, spikeSites, spikeTimes, P)
    nSpikes = numel(spikeTimes);
    nSites = numel(P.chanMap);
    waveformLimits = P.spkLim;
    traceLimits = P.spkLim_raw;
    nSites_spk = (P.maxSite * 2) + 1;
    spikeTraces = zeros(diff(traceLimits) + 1, nSites_spk, nSpikes, 'like', rawTraces);
    spikeWaveforms = zeros(diff(waveformLimits) + 1, nSites_spk, nSpikes, 'like', filteredTraces);

    % Realignment parameters
    fRealign_spk = getOr(P, 'fRealign_spk', 0); %0,1,2
    spikeTimes = gpuArray_(spikeTimes, isGpu_(rawTraces));
    spikeSites = gpuArray_(spikeSites, isGpu_(rawTraces));

    if isempty(spikeSites)
        spikeTraces = permute(mr2tr3_(rawTraces, traceLimits, spikeTimes), [1, 3, 2]);
        spikeWaveforms = permute(mr2tr3_(filteredTraces, waveformLimits, spikeTimes), [1, 3, 2]);
    else
        for iSite = 1:nSites
            spikesThisSite = find(spikeSites == iSite);
            if isempty(spikesThisSite)
                continue;
            end

            spikeTimesThisSite = spikeTimes(spikesThisSite); %already sorted by time
            siteNeighborhood = P.miSites(:, iSite);

            try
                spikeWaveforms1 = mr2tr3_(filteredTraces, waveformLimits, spikeTimesThisSite, siteNeighborhood);
                if fRealign_spk == 1
                    [spikeWaveforms1, spikeTimesThisSite] = spkwav_realign_(spikeWaveforms1, filteredTraces, waveformLimits, spikeTimesThisSite, siteNeighborhood, P);
                    spikeTimes(spikesThisSite) = spikeTimesThisSite;
                elseif fRealign_spk == 2
                    spikeWaveforms1 = spkwav_align_(spikeWaveforms1, P);
                end

                spikeWaveforms(:, :, spikesThisSite) = permute(spikeWaveforms1, [1,3,2]);
                spikeTraces(:, :, spikesThisSite) = permute(mr2tr3_(rawTraces, traceLimits, spikeTimesThisSite, siteNeighborhood), [1,3,2]); %raw
            catch % GPU failure
                disperr_('mn2tn_wav_: GPU failed');
            end
        end
    end
end % function
