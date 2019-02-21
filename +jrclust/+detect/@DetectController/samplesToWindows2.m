function spikeWindows2 = samplesToWindows2(obj, samplesIn, spikeSites, spikeTimes)
    %SAMPLESTOWINDOWS2 Get spatiotemporal windows around SECONDARY peaking
    % events as 3D arrays (filtered samples assumed)
    nSpikes = numel(spikeSites);
    nSitesEvt = 1 + obj.hCfg.nSiteDir*2; % includes ref sites

    % nSamples x nSites x nSpikes
    spikeWindows2 = zeros(diff(obj.hCfg.evtWindowSamp) + 1, nSitesEvt, nSpikes, 'like', samplesIn);

    for iSite = 1:obj.hCfg.nSites
        siteSpikes = find(spikeSites == iSite);
        if isempty(siteSpikes)
            continue;
        end

        siteTimes = spikeTimes(siteSpikes); % already sorted by time
        iNeighbors = obj.hCfg.siteNeighbors(:, iSite);
        spikeWindows2(:, :, siteSpikes) = permute(jrclust.utils.tryGather(obj.extractWindows(samplesIn, siteTimes, iNeighbors, 0)), [1,3,2]);
    end
end