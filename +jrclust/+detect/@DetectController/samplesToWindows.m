function spikeData = samplesToWindows(obj, spikeData)
    %SAMPLESTOWINDOWS Get spatiotemporal windows around spiking events as 3D arrays
    %   returns spikesRaw nSamples x nSites x nSpikes
    samplesRaw = spikeData.samplesRaw;
    samplesFilt = spikeData.samplesFilt;
    spikeTimes = spikeData.spikeTimes;
    spikeSites = spikeData.spikeSites;

    nSpikes = numel(spikeTimes);
    nSitesEvt = 1 + obj.hCfg.nSiteDir*2; % includes ref sites

    % tensors, nSamples{Raw, Filt} x nSites x nSpikes
    spikesRaw = zeros(diff(obj.hCfg.evtWindowRawSamp) + 1, nSitesEvt, nSpikes, 'like', samplesRaw);
    spikesFilt = zeros(diff(obj.hCfg.evtWindowSamp) + 1, nSitesEvt, nSpikes, 'like', samplesFilt);

    % Realignment parameters
    realignTraces = obj.hCfg.getOr('realignTraces', 0); % 0, 1, 2
    spikeTimes = jrclust.utils.tryGpuArray(spikeTimes, obj.hCfg.useGPU);
    spikeSites = jrclust.utils.tryGpuArray(spikeSites, obj.hCfg.useGPU);

    % extractWindows returns nSamples x nSpikes x nSites
    if isempty(spikeSites)
        spikesRaw = permute(obj.extractWindows(samplesRaw, spikeTimes, [], 1), [1, 3, 2]);
        spikesFilt = permute(obj.extractWindows(samplesFilt, spikeTimes, [], 0), [1, 3, 2]);
    else
        for iSite = 1:obj.hCfg.nSites
            siteSpikes = find(spikeSites == iSite);
            if isempty(siteSpikes)
                continue;
            end

            siteTimes = spikeTimes(siteSpikes); % already sorted by time
            neighbors = obj.hCfg.siteNeighbors(:, iSite);

            try
                spikeWindows = obj.extractWindows(samplesFilt, siteTimes, neighbors, 0);

                if realignTraces == 1
                    [spikeWindows, siteTimes] = obj.CARRealign(spikeWindows, samplesFilt, siteTimes, neighbors);
                    spikeTimes(siteSpikes) = siteTimes;
                elseif realignTraces == 2
                    spikeWindows = interpPeaks(spikeWindows, obj.hCfg);
                end

                spikesFilt(:, :, siteSpikes) = permute(spikeWindows, [1, 3, 2]);
                spikesRaw(:, :, siteSpikes) = permute(obj.extractWindows(samplesRaw, siteTimes, neighbors, 1), [1, 3, 2]);
            catch ME
                obj.errMsg = ME.message;
                obj.isError = 1;
            end
        end
    end

    % get windows around secondary/tertiary peaks as well
    centerSites = spikeSites(:);
    if obj.hCfg.nPeaksFeatures >= 2
        [spikeSites2, spikeSites3] = obj.findSecondaryPeaks(spikesFilt, spikeSites);
        spikesFilt2 = obj.samplesToWindows2(samplesFilt, spikeSites2, spikeTimes);
        centerSites = [spikeSites(:) spikeSites2(:)];
    else
        spikesFilt2 = [];
    end
    if obj.hCfg.nPeaksFeatures == 3
        spikesFilt3 = obj.samplesToWindows2(samplesFilt, spikeSites3, spikeTimes);
        centerSites = [centerSites spikeSites3(:)];
    else
        spikesFilt3 = [];
    end
    if obj.hCfg.getOr('fCancel_overlap', 0)
        try
            [spikesFilt, spikesFilt2] = obj.cancelOverlap(spikesFilt, spikesFilt2, spikeTimes, spikeSites, spikeSites2, siteThresh);
        catch ME
            warning('Cancel overlap failed: %s', ME.message);
        end
    end

    spikeData.spikesRaw = spikesRaw;
    spikeData.spikesFilt = spikesFilt;
    spikeData.spikeTimes = jrclust.utils.tryGather(spikeTimes);
    spikeData.centerSites = jrclust.utils.tryGather(centerSites);
    spikeData.spikesFilt2 = spikesFilt2;
    spikeData.spikesFilt3 = spikesFilt3;
end

%% LOCAL FUNCTIONS
function spikeWindows = interpPeaks(spikeWindows, hCfg)
    %INTERPPEAKS Interpolate waveforms to determine if true peaks are at subpixel locations
    nInterp = 2; % interpolate at midpoint between each sample

    % interpolate on just the primary site
    spikeWindowsInterp = jrclust.utils.interpWindows(spikeWindows(:, :, 1), nInterp);
    [~, interpPeaks] = min(spikeWindowsInterp);

    peakLoc = 1 - hCfg.evtWindowSamp(1)*nInterp; % peak loc of waveform before interpolation
    isRightShifted = find(interpPeaks == peakLoc + 1); % peak is subpixel right of given loc
    isLeftShifted = find(interpPeaks == peakLoc - 1); % peak is subpixel left of given loc

    % replace right-shifted peaks with their interpolated values
    if ~isempty(isRightShifted)
        rightShifted = jrclust.utils.interpWindows(spikeWindows(:, isRightShifted, :), nInterp);
        spikeWindows(1:end-1, isRightShifted, :) = rightShifted(2:nInterp:end, :, :);
    end

    % replace left-shifted peaks with their interpolated values
    if ~isempty(isLeftShifted)
        leftShifted = jrclust.utils.interpWindows(spikeWindows(:, isLeftShifted, :), nInterp);
        spikeWindows(2:end, isLeftShifted, :) = leftShifted(2:nInterp:end, :, :);
    end
end
