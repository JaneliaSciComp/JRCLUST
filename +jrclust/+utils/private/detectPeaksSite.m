function [peakLocs, peaks, siteThresh] = detectPeaksSite(samplesIn, siteThresh, hCfg)
    %DETECTPEAKSSITE Detect peaks on a given site, computing threshold if necessary

    % determine threshold if it has still escaped us
    MAX_SAMPLE_QQ = 300000;
    if ~isempty(hCfg.evtManualThreshSamp)
        siteThresh = hCfg.evtManualThreshSamp;
    end
    if siteThresh == 0 % bad site
        [peakLocs, peaks] = deal([]);
        return;
    end
    if isempty(siteThresh)
        siteThresh = hCfg.qqFactor*jrclust.utils.estimateRMS(samplesIn, MAX_SAMPLE_QQ);
    end
    siteThresh = cast(siteThresh, 'like', samplesIn);

    % detect turning point in waveforms exceeding threshold
    peakLocs = findPeaks(samplesIn, siteThresh, hCfg.minNeighborsDetect);

    if hCfg.detectBipolar % detect positive peaks
        peakLocs = [peakLocs; findPeaks(-samplesIn, siteThresh, hCfg.minNeighborsDetect)];
        peakLocs = sort(peakLocs);
    end

    if isempty(peakLocs)
        peakLocs = double([]);
        peaks = int16([]);
    else
        peaks = samplesIn(peakLocs);

        % remove overlarge spikes
        if ~isempty(hCfg.spikeThreshMax)
            threshMax = cast(abs(hCfg.spikeThreshMax)/hCfg.bitScaling, 'like', peaks);
            keepMe = (abs(peaks) < abs(threshMax));
            peakLocs = peakLocs(keepMe);
            peaks = peaks(keepMe);
        end
    end

    peakLocs = jrclust.utils.tryGather(peakLocs);
    peaks = jrclust.utils.tryGather(peaks);
    siteThresh = jrclust.utils.tryGather(siteThresh);
end

%% LOCAL FUNCTIONS
function peaks = findPeaks(samplesIn, thresh, nneighBelow)
    %FINDPEAKS Find samples which exceed threshold
    if isempty(nneighBelow)
        nneighBelow = 1;
    end

    peaks = [];
    if isempty(samplesIn)
        return;
    end

    exceedsThresh = jrclust.utils.tryGather(samplesIn < -abs(thresh));
    peakLocs = find(exceedsThresh);
    if isempty(peakLocs)
        return;
    end

    % got a peak at 1st sample, incomplete waveform (+ indexing error)
    if peakLocs(1) <= 1 && numel(peakLocs) == 1 % only peak found, give up
        return;
    elseif peakLocs(1) <= 1 % discard this one, salvage the others
        peakLocs(1) = [];
    end

    % got a peak at last sample, incomplete waveform (+ indexing error)
    if peakLocs(end) >= numel(samplesIn) && numel(peakLocs) == 1 % only peak found, give up
        return;
    else % discard this one, salvage the others
        peakLocs(end) = [];
    end

    peakCenters = samplesIn(peakLocs);
    % take only "peaks" whose sample neighbors are not larger
    peaks = peakLocs(peakCenters <= samplesIn(peakLocs + 1) & peakCenters <= samplesIn(peakLocs - 1));
    if isempty(peaks)
        return;
    end

    % take only peaks who have one or two sample neighbors exceeding threshold
    if nneighBelow == 1
        peaks = peaks(exceedsThresh(peaks - 1) | exceedsThresh(peaks + 1));
    elseif nneighBelow == 2
        peaks = peaks(exceedsThresh(peaks - 1) & exceedsThresh(peaks + 1));
    end
end
