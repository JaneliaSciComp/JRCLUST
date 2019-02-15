function [spikesOut, sitesOut] = getSpikesInBounds(spikesIn, sitesIn, yPos, yLims)
    %GETSPIKESINBOUNDS Select spikes by whether their y-positions fall within given limits
    nSamplesMax = 1000;
    inBounds = (yPos >= yLims(1)) & (yPos < yLims(2));

    if ~any(inBounds)
        spikesOut = spikesIn;
        sitesOut = sitesIn;
        return;
    end

    spikesOut = jrclust.utils.subsample(spikesIn(inBounds), nSamplesMax);
    sitesOut = jrclust.utils.subsample(sitesIn(inBounds), nSamplesMax);
end