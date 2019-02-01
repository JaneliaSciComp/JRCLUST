function [windows, timeRanges] = extractWindows(samplesIn, windowLim, spTimes, spSites, fMeanSubt)
    %EXTRACTWINDOWS Isolate windows around spiking events from matrix of samples
    %   returns nSamples x nSpikes x nSites
    if nargin < 4
        spSites = [];
    end
    if nargin < 5
        fMeanSubt = 0;
    end

    if isempty(spTimes)
        windows = [];
        return;
    end

    [nSamples, nSites] = size(samplesIn);
    if ~isempty(spSites)
        nSites = numel(spSites);
    end

    % if iscolumn(spTimes)
    %     spTimes = spTimes';
    % end
    spTimes = reshape(spTimes, 1, numel(spTimes));

    winRange = (windowLim(1):windowLim(end))';
    timeRanges = bsxfun(@plus, int32(winRange), int32(spTimes));
    timeRanges = min(max(timeRanges, 1), nSamples); % ensure we don't go out of bounds
    timeRanges = timeRanges(:);

    if isempty(spSites)
        windows = samplesIn(timeRanges, :);
    else
        windows = samplesIn(timeRanges, spSites);
    end

    windows = reshape(windows, [numel(winRange), numel(spTimes), nSites]);

    if fMeanSubt
        windows = single(windows);
        oldDim = size(windows);
        windows = reshape(windows, size(windows,1), []);
        windows = bsxfun(@minus, windows, mean(windows));
        windows = reshape(windows, oldDim);
    end
end
