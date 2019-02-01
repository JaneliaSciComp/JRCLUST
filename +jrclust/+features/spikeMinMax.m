function [spikeMin, spikeMax, vpp] = spikeMinMax(spikeWindows)
    %SPIKEMINMAX Summary of this function goes here
    spikeMin = shiftdim(min(spikeWindows))';
    spikeMax = shiftdim(max(spikeWindows))';

    if nargout == 3
        vpp = spikeMax - spikeMin;
    end

    spikeMin = abs(spikeMin);
    spikeMax = abs(spikeMax);
end

