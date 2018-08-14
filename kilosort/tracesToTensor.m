function spikeTraces = tracesToTensor(traces, spikeSites, spikeTimes, traceLimits, P)
    nSpikes = numel(spikeTimes);
    nSites = numel(P.chanMap);
    nSites_spk = (P.maxSite * 2) + 1;
    spikeTraces = zeros(diff(traceLimits) + 1, nSites_spk, nSpikes, 'like', traces);

    % Realignment parameters
    fRealign_spk = get_set_(P, 'fRealign_spk', 0); %0,1,2
    spikeTimes = gpuArray_(spikeTimes, isGpu_(traces));
    spikeSites = gpuArray_(spikeSites, isGpu_(traces));

    if isempty(spikeSites)
        spikeTraces = permute(mr2tr3_(traces, traceLimits, spikeTimes), [1, 3, 2]);
    else
        for iSite = 1:nSites
            spikesThisSite = find(spikeSites == iSite);
            if isempty(spikesThisSite)
                continue;
            end

            spikeTimesThisSite = spikeTimes(spikesThisSite); %already sorted by time
            siteNeighborhood = P.miSites(:, iSite);

            try
                spikeTraces(:, :, spikesThisSite) = permute(mr2tr3_(traces, traceLimits, spikeTimesThisSite, siteNeighborhood), [1, 3, 2]);
            catch % GPU failure
                disperr_('tracesToTensor: GPU failed');
            end
        end
    end
end
