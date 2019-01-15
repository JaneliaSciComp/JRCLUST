function spikeData = findPeaks(spikeData, hCfg)
    %FINDPEAKS Detect or import spike peaks
    samplesFilt = spikeData.samplesFilt;
    keepMe = spikeData.keepMe;
    spikeTimes = spikeData.spikeTimes;
    spikeSites = spikeData.spikeSites;
    siteThresh = spikeData.siteThresh;
    nPadPre = spikeData.nPadPre;
    nPadPost = spikeData.nPadPost;

    if isempty(siteThresh)
        try
            siteThresh = jrclust.utils.tryGather(int16(jrclust.utils.estimateRMS(samplesFilt, 1e5)*hCfg.qqFactor));
        catch ME
            hCfg.useGPU = 0;
            siteThresh = int16(jrclust.utils.estimateRMS(jrclust.utils.tryGather(samplesFilt), 1e5)*hCfg.qqFactor);
        end
    end

    if isempty(spikeTimes) || isempty(spikeSites)
        if ~isprop(hCfg, 'nPadPre')
            hCfg.addprop('nPadPre');
        end
        hCfg.nPadPre = nPadPre;

        [spikeTimes, spikeAmps, spikeSites] = jrclust.utils.detectPeaks(samplesFilt, siteThresh, keepMe, hCfg);
    else
        spikeTimes = spikeTimes + nPadPre;
        spikeAmps = samplesFilt(sub2ind(size(samplesFilt), spikeTimes, spikeSites)); % @TODO read spikes at the site and time
    end

    % reject spikes within the overlap region
    if nPadPre > 0 || nPadPost > 0
        bounds = [nPadPre + 1, size(samplesFilt, 1) - nPadPost]; % inclusive
        inBounds = spikeTimes >= bounds(1) & spikeTimes <= bounds(2);

        spikeTimes = spikeTimes(inBounds);
        spikeAmps  = spikeAmps(inBounds);
        spikeSites = spikeSites(inBounds);
    end

    spikeData.siteThresh = siteThresh;
    spikeData.spikeTimes = spikeTimes;
    spikeData.spikeAmps = jrclust.utils.tryGather(spikeAmps);
    spikeData.spikeSites = spikeSites;
end