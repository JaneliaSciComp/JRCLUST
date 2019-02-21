function [spikeWindows, spikeTimes] = CARRealign(obj, spikeWindows, samplesIn, spikeTimes, neighbors)
    %CARREALIGN Realign spike peaks after applying local CAR
    if ~strcmpi(obj.hCfg.getOr('vcSpkRef', 'nmean'), 'nmean')
        return;
    end

    % find where true peaks are not in the correct place after applying local CAR
    spikeWindowsCAR = jrclust.utils.localCAR(single(spikeWindows), obj.hCfg);
    [shiftMe, shiftBy] = findShifted(spikeWindowsCAR, obj.hCfg);

    if isempty(shiftMe)
        return;
    end

    % adjust spike times
    shiftedTimes = spikeTimes(shiftMe) - int32(shiftBy(:));
    spikeTimes(shiftMe) = shiftedTimes;

    % extract windows at new shifted times
    spikeWindows(:, shiftMe, :) = obj.extractWindows(samplesIn, shiftedTimes, neighbors, 0);
end

%% LOCAL FUNCTIONS
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