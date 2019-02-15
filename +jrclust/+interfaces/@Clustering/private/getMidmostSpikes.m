function subSpikes = getMidmostSpikes(spikes, spikeTimes, nClusterIntervals)
    %GETMIDMOSTSPIKES Get the spikes closest to the center spike
    middlemost = round(numel(spikeTimes)/2); % index of middlemost spike
    nearestToCenter = jrclust.utils.rankorder(abs(spikes - middlemost), 'ascend');
    nSpikesInterval = round(numel(spikes) / nClusterIntervals);
    subSpikes = spikes(nearestToCenter <= nSpikesInterval);
end