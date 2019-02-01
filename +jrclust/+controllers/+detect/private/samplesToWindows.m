function [spikesRaw, spikesFilt, spikeTimes] = samplesToWindows(samplesRaw, samplesFilt, spikeSites, spikeTimes, hCfg)
    %SAMPLESTOWINDOWS Get spatiotemporal windows around spiking events as 3D arrays
    %   returns spikesRaw nSamples x nSites x nSpikes
    nSpikes = numel(spikeTimes);
    nSitesEvt = 1 + hCfg.nSiteDir*2; % includes ref sites

    % tensors, nSamples{Raw, Filt} x nSites x nSpikes
    spikesRaw = zeros(diff(hCfg.evtWindowRawSamp) + 1, nSitesEvt, nSpikes, 'like', samplesRaw);
    spikesFilt = zeros(diff(hCfg.evtWindowSamp) + 1, nSitesEvt, nSpikes, 'like', samplesFilt);

    % Realignment parameters
    realignTraces = hCfg.getOr('realignTraces', 0); % 0, 1, 2
    spikeTimes = jrclust.utils.tryGpuArray(spikeTimes, isa(samplesRaw, 'gpuArray'));
    spikeSites = jrclust.utils.tryGpuArray(spikeSites, isa(samplesRaw, 'gpuArray'));

    % extractWindows returns nSamples x nSpikes x nSites
    if isempty(spikeSites)
        spikesRaw = permute(jrclust.utils.extractWindows(samplesRaw, hCfg.evtWindowRawSamp, spikeTimes), [1, 3, 2]);
        spikesFilt = permute(jrclust.utils.extractWindows(samplesFilt, hCfg.evtWindowSamp, spikeTimes), [1, 3, 2]);
    else
        for iSite = 1:hCfg.nSites
            siteSpikes = find(spikeSites == iSite);
            if isempty(siteSpikes)
                continue;
            end

            siteTimes = spikeTimes(siteSpikes); % already sorted by time
            neighbors = hCfg.siteNeighbors(:, iSite);

            try
                spikeWindows = jrclust.utils.extractWindows(samplesFilt, hCfg.evtWindowSamp, siteTimes, neighbors);

                if realignTraces == 1
                    [spikeWindows, siteTimes] = carRealign(spikeWindows, samplesFilt, siteTimes, neighbors, hCfg);
                    spikeTimes(siteSpikes) = siteTimes;
                elseif realignTraces == 2
                    spikeWindows = interpPeaks(spikeWindows, hCfg);
                end

                spikesFilt(:, :, siteSpikes) = permute(spikeWindows, [1, 3, 2]);
                spikesRaw(:, :, siteSpikes) = permute(jrclust.utils.extractWindows(samplesRaw, hCfg.evtWindowRawSamp, siteTimes, neighbors), [1, 3, 2]);
            catch ME
                obj.errMsg = ME.message;
                obj.isError = 1;
            end
        end
    end
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

function [spikeWindows, spikeTimes] = carRealign(spikeWindows, samplesIn, spikeTimes, neighbors, hCfg)
    %CARREALIGN Realign spike peaks
    if ~strcmpi(hCfg.getOr('vcSpkRef', 'nmean'), 'nmean')
        return;
    end

    % find where true peaks are not in the correct place after applying CAR
    spikeWindowsCAR = jrclust.utils.localCAR(single(spikeWindows), hCfg); % apply LCAR
    [shiftMe, shiftBy] = findShifted(spikeWindowsCAR, hCfg);

    if isempty(shiftMe)
        return;
    end

    % adjust spike times
    shiftedTimes = spikeTimes(shiftMe) - int32(shiftBy(:));
    spikeTimes(shiftMe) = shiftedTimes;

    % extract windows at new shifted times
    spikeWindows(:, shiftMe, :) = jrclust.utils.extractWindows(samplesIn, hCfg.evtWindowSamp, shiftedTimes, neighbors);
end

function [shiftMe, shiftBy] = findShifted(spikeWindows, hCfg)
    %FINDSHIFTED
    %   spikeWindows: nSamples x nSpikes x nSites
    peakLoc = 1 - hCfg.evtWindowSamp(1);

    if hCfg.detectBipolar
        [~, truePeakLoc] = max(abs(spikeWindows(:, :, 1)));
    else
        [~, truePeakLoc] = min(spikeWindows(:, :, 1));
    end

    shiftMe = find(truePeakLoc ~= peakLoc);
    shiftBy = peakLoc - truePeakLoc(shiftMe);

    shiftOkay = (abs(shiftBy) <= 2); % throw out drastic shifts
    shiftMe = shiftMe(shiftOkay);
    shiftBy = shiftBy(shiftOkay);
end
