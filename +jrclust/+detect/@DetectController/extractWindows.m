function [windows, timeRanges] = extractWindows(obj, samplesIn, spikeTimes, spikeSites, fRaw)
    %EXTRACTWINDOWS Isolate windows around spiking events from matrix of samples
    %   returns nSamples x nSpikes x nSites
    if nargin < 4
        spikeSites = [];
    end
    if nargin < 5
        fRaw = 0;
    end

    if isempty(spikeTimes)
        windows = [];
        return;
    end

    [nSamples, nSites] = size(samplesIn);
    if ~isempty(spikeSites)
        nSites = numel(spikeSites);
    end

    spikeTimes = spikeTimes(:)'; % make a row vector

    if fRaw
        winRange = (obj.hCfg.evtWindowRawSamp(1):obj.hCfg.evtWindowRawSamp(end))';
    else
        winRange = (obj.hCfg.evtWindowSamp(1):obj.hCfg.evtWindowSamp(end))';
    end

    timeRanges = bsxfun(@plus, int32(winRange), int32(spikeTimes));
    timeRanges = min(max(timeRanges, 1), nSamples); % ensure we don't go out of bounds
    timeRanges = timeRanges(:);

    if isempty(spikeSites)
        windows = samplesIn(timeRanges, :);
    else
        windows = samplesIn(timeRanges, spikeSites);
    end

    windows = reshape(windows, [numel(winRange), numel(spikeTimes), nSites]);
end
